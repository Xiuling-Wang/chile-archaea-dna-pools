#!/usr/bin/env Rscript

# Supplementary taxonomic composition figure: phylum and genus panels only.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
house_style_run <- identical(Sys.getenv("ROUND48_HOUSESTYLE_OUTPUT"), "1")
theme_path <- file.path(project_root, "scripts", "theme_metagenomics.R")
if (file.exists(theme_path)) {
  source(theme_path)
}

if (!exists("theme_pub")) {
  theme_pub <- function(base_size = 8, base_family = "Helvetica") {
    theme_classic(base_size = base_size, base_family = base_family) +
      theme(
        axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.3, colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.key = element_blank()
      )
  }
}

figure_dir <- file.path(project_root, "figures")
candidate_dir <- file.path(figure_dir, "candidates")
analysis_dir <- if (house_style_run) {
  file.path(tempdir(), "chile_archaea_round48_housestyle", "FigS2")
} else {
  file.path(project_root, "analysis", "reframed_figures")
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

raw_path <- file.path(project_root, "data", "arc_unrarefild", "ASV_Arc_200cm_delete_less200.txt")
retained_samples_path <- file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt")
tax_path <- file.path(project_root, "data", "file_A.txt")
env_path <- file.path(project_root, "data", "env_2024.csv")

missing_inputs <- c(raw_path, retained_samples_path, tax_path, env_path)[!file.exists(c(raw_path, retained_samples_path, tax_path, env_path))]
if (length(missing_inputs) > 0) {
  stop("Missing input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

site_levels <- c("AZ", "SG", "LC", "NB")
site_display <- c(AZ = "AZ", SG = "SG", LC = "LC", NB = "NA")
dna_levels <- c("iDNA", "eDNA")
depth_levels_raw <- c("12_0_5", "11_5_10", "10_10_20", "09_20_40", "08_40_60")
depth_labels <- c(
  "12_0_5" = "0-5",
  "11_5_10" = "5-10",
  "10_10_20" = "10-20",
  "09_20_40" = "20-40",
  "08_40_60" = "40-60"
)
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
other_label <- "Other/unclassified/unassigned"

display_site <- function(x) unname(site_display[as.character(x)])

clean_tax <- function(x) {
  x <- iconv(as.character(x), from = "", to = "UTF-8", sub = "")
  x <- str_replace(x, "^[a-z]__", "")
  x <- str_replace(x, "_(cl|or|fa|ge|sp)$", "")
  x[is.na(x) | x == "" | x == "NA"] <- "Unclassified"
  x
}

is_other_taxon <- function(x) {
  str_detect(str_to_lower(x), "^(unclassified|unassigned|other)$")
}

raw <- read.delim(raw_path, header = TRUE, row.names = 1, check.names = FALSE)
retained_samples <- colnames(read.delim(retained_samples_path, row.names = 1, check.names = FALSE))
tax_source <- read.delim(tax_path, header = TRUE, check.names = FALSE)
env <- read.csv(env_path, header = TRUE, check.names = FALSE) %>%
  filter(sample_id %in% retained_samples, depth < 60, site %in% site_levels, dna_type %in% dna_levels) %>%
  mutate(
    site = factor(site, levels = site_levels),
    site_label = factor(display_site(site), levels = site_display[site_levels]),
    dna_type = factor(dna_type, levels = dna_levels),
    depth_plot = factor(depths_cm_re, levels = rev(depth_levels_raw), labels = rev(depth_labels[depth_levels_raw]))
  )

tax <- tax_source %>%
  filter(asv_id %in% rownames(raw)) %>%
  select(asv_id, all_of(tax_cols))
counts <- raw

matched_samples <- retained_samples[retained_samples %in% env$sample_id & retained_samples %in% colnames(counts)]
if (length(matched_samples) != length(retained_samples)) {
  stop("Expected all 95 retained 0-60 cm samples in the unrarefied table and metadata.")
}
counts <- counts[, matched_samples, drop = FALSE]
counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
tax <- tax %>% filter(asv_id %in% rownames(counts))
env <- env %>% filter(sample_id %in% colnames(counts))

asv_long <- counts %>%
  as.data.frame() %>%
  rownames_to_column("asv_id") %>%
  pivot_longer(-asv_id, names_to = "sample_id", values_to = "reads") %>%
  left_join(tax, by = "asv_id") %>%
  left_join(env %>% select(sample_id, site, site_label, dna_type, depth_plot), by = "sample_id") %>%
  filter(!is.na(site)) %>%
  group_by(sample_id) %>%
  mutate(sample_rel = reads / sum(reads, na.rm = TRUE)) %>%
  ungroup()

make_rank_data <- function(rank, top_n = Inf) {
  base <- asv_long %>%
    mutate(taxon_raw = clean_tax(.data[[rank]]))

  totals <- base %>%
    group_by(taxon_raw) %>%
    summarise(total_abundance = sum(sample_rel, na.rm = TRUE), .groups = "drop") %>%
    mutate(is_other = is_other_taxon(taxon_raw)) %>%
    arrange(desc(total_abundance), taxon_raw)

  keep_taxa <- totals %>%
    filter(!is_other) %>%
    slice_head(n = top_n) %>%
    pull(taxon_raw)

  grouped <- base %>%
    mutate(taxon = if_else(taxon_raw %in% keep_taxa & !is_other_taxon(taxon_raw), taxon_raw, other_label))

  level_order <- grouped %>%
    group_by(taxon) %>%
    summarise(total_abundance = sum(sample_rel, na.rm = TRUE), .groups = "drop") %>%
    filter(taxon != other_label) %>%
    arrange(desc(total_abundance), taxon) %>%
    pull(taxon)
  level_order <- c(level_order, other_label)

  grouped %>%
    group_by(sample_id, site, site_label, dna_type, depth_plot, taxon) %>%
    summarise(prop = sum(sample_rel, na.rm = TRUE), .groups = "drop") %>%
    complete(
      nesting(sample_id, site, site_label, dna_type, depth_plot),
      taxon = level_order,
      fill = list(prop = 0)
    ) %>%
    group_by(site, site_label, dna_type, depth_plot, taxon) %>%
    summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      signed_prop = if_else(dna_type == "iDNA", -prop, prop),
      taxon = factor(taxon, levels = level_order)
    )
}

make_taxon_cols <- function(levels_in) {
  non_other <- setdiff(levels_in, other_label)
  palette <- c(
    "#B2182B", "#2166AC", "#1B9E77", "#E08214", "#7570B3", "#66A61E",
    "#C51B7D", "#8C510A", "#35978F", "#542788", "#D95F02", "#4D4D4D",
    "#A6761D", "#1F78B4", "#B2DF8A", "#FB9A99"
  )
  cols <- setNames(rep(palette, length.out = length(non_other)), non_other)
  c(cols, setNames("#BDBDBD", other_label))
}

make_rank_plot <- function(dat, rank_label) {
  fill_levels <- levels(dat$taxon)
  ggplot(dat, aes(x = signed_prop, y = depth_plot, fill = taxon)) +
    geom_col(width = 0.70, colour = "white", linewidth = 0.05, position = position_stack(reverse = TRUE)) +
    geom_vline(xintercept = 0, colour = "white", linewidth = 0.90) +
    facet_wrap(~ site_label, nrow = 1) +
    scale_fill_manual(values = make_taxon_cols(fill_levels), breaks = fill_levels, drop = FALSE) +
    scale_x_continuous(
      position = "top",
      breaks = seq(-1, 1, by = 0.5),
      labels = function(x) abs(x) * 100,
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(xlim = c(-1, 1), clip = "off") +
    labs(
      title = paste0(rank_label, "-level archaeal composition"),
      subtitle = "iDNA is left of zero and eDNA is right of zero.",
      x = "Relative abundance (%)",
      y = "Depth (cm)",
      fill = rank_label
    ) +
    theme_pub(base_size = 9.4, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 11.0),
      plot.subtitle = element_text(size = 8.2, colour = "grey25", margin = margin(b = 2)),
      strip.text = element_text(face = "bold", size = 9.9),
      axis.text.x = element_text(size = 8.2, angle = 45, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 8.8),
      axis.title = element_text(size = 9.8, face = "bold"),
      panel.grid.major.x = element_line(colour = "grey87", linewidth = 0.20),
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(1.2, "lines"),
      legend.position = "bottom",
      legend.title = element_text(size = 8.9, face = "bold"),
      legend.text = element_text(size = 8.0),
      legend.key.size = unit(0.40, "cm")
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

phylum_data <- make_rank_data("Phylum", top_n = Inf)
genus_data <- make_rank_data("Genus", top_n = 12)

write.csv(phylum_data, file.path(analysis_dir, "SI_taxonomic_composition_phylum_summary.csv"), row.names = FALSE)
write.csv(genus_data, file.path(analysis_dir, "SI_taxonomic_composition_genus_summary.csv"), row.names = FALSE)

p_phylum <- make_rank_plot(phylum_data, "Phylum")
p_genus <- make_rank_plot(genus_data, "Genus")

si_taxonomy <- p_phylum / p_genus +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold", size = 20, family = "Helvetica", colour = "black"
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

out_pdf <- if (house_style_run) {
  file.path(candidate_dir, "FigS2_round48_housestyle.pdf")
} else {
  file.path(figure_dir, "SI_taxonomic_composition.pdf")
}
ggsave(out_pdf, si_taxonomy, width = 11.6, height = 8.1, units = "in", device = cairo_pdf, bg = "white")

message("Wrote: ", out_pdf)
