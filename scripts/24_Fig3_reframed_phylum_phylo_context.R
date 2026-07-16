#!/usr/bin/env Rscript

# Reframed Figure 3: keep the established mirrored phylum-composition grammar,
# then add a phylogenetic panel that still carries site x DNA-pool abundance
# information through bubbles. The tree should not be decorative.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ape)
  library(ggtree)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
theme_path <- file.path(project_root, "scripts", "theme_metagenomics.R")
if (file.exists(theme_path)) {
  source(theme_path)
}

if (!exists("theme_pub")) {
  theme_pub <- function(base_size = 8, base_family = "Helvetica") {
    theme_classic(base_size = base_size, base_family = base_family)
  }
}

analysis_dir <- file.path(project_root, "analysis", "reframed_figures")
figure_dir <- file.path(project_root, "figures")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_analysis_outputs <- identical(Sys.getenv("FIG3_WRITE_ANALYSIS_OUTPUTS"), "1")

asv_path <- file.path(project_root, "data", "arc_unrarefild", "ASV_Arc_200cm_delete_less200.txt")
retained_samples_path <- file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt")
taxa_path <- file.path(project_root, "data", "file_A.txt")
env_path <- file.path(project_root, "data", "env_2024.csv")
tree_path <- file.path(project_root, "data", "tree arch", "Tree archaea 1.tree")
ann_path <- file.path(project_root, "analysis", "tree_audit", "legacy_tree_asv_annotation.csv")

missing_inputs <- c(asv_path, retained_samples_path, taxa_path, env_path, tree_path, ann_path)[
  !file.exists(c(asv_path, retained_samples_path, taxa_path, env_path, tree_path, ann_path))
]
if (length(missing_inputs) > 0) {
  stop("Missing input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

site_levels <- c("AZ", "SG", "LC", "NB")
site_display <- c(AZ = "AZ", SG = "SG", LC = "LC", NB = "NA")
pool_levels <- c("iDNA", "eDNA")
site_pool_levels <- as.vector(t(outer(site_levels, pool_levels, paste, sep = "_")))
site_pool_labels <- paste0(site_display[sub("_.*$", "", site_pool_levels)], " ", sub("^.*_", "", site_pool_levels))
depth_levels_raw <- c("12_0_5", "11_5_10", "10_10_20", "09_20_40", "08_40_60")
depth_labels <- c(
  "12_0_5" = "0-5",
  "11_5_10" = "5-10",
  "10_10_20" = "10-20",
  "09_20_40" = "20-40",
  "08_40_60" = "40-60"
)

composition_class_cols <- c(
  "Nitrososphaeria (AOA-affiliated)" = "#B2182B",
  Thermoplasmata = "#542788",
  "Thermoplasmatota (other/unclassified)" = "#E08214",
  "Other archaea" = "#BDBDBD"
)
composition_class_levels <- names(composition_class_cols)
site_cols <- c(AZ = "#D95F02", SG = "#1B9E77", LC = "#7570B3", NB = "#66A61E")
class_cols <- c(
  Nitrososphaeria = "#B2182B",
  Thermoplasmata = "#542788",
  Thermoplasmatota = "#E08214",
  Thermoprotei = "#1B7837",
  Bathyarchaeia = "#C51B7D",
  Nanoarchaeia = "#4D4D4D",
  unclassified = "#6B8BA4"
)

display_site <- function(x) unname(site_display[as.character(x)])

clean_tax <- function(x) {
  x <- iconv(as.character(x), from = "", to = "UTF-8", sub = "")
  x <- sub("^[a-z]__", "", x)
  x[is.na(x) | x == "" | x == "NA"] <- "Unassigned"
  x <- sub("_(cl|or|fa|ge|sp)$", "", x)
  x
}

clean_tax_lower_unclassified <- function(x) {
  x <- clean_tax(x)
  x[x == "Unassigned"] <- "unclassified"
  x
}

extract_newick <- function(path) {
  x <- readLines(path, warn = FALSE)
  start <- grep("^\\s*\\(", x)[1]
  if (is.na(start)) stop("No Newick tree found in ", path)
  nwk_lines <- x[start:length(x)]
  end <- grep(";\\s*$", nwk_lines)[1]
  if (is.na(end)) stop("No terminating semicolon found in ", path)
  paste(nwk_lines[seq_len(end)], collapse = "")
}

asv_counts <- read.delim(asv_path, row.names = 1, check.names = FALSE)
retained_samples <- colnames(read.delim(retained_samples_path, row.names = 1, check.names = FALSE))
taxa <- read.delim(taxa_path, check.names = FALSE) %>%
  mutate(
    asv_id = as.character(asv_id),
    Phylum_clean = clean_tax(Phylum),
    Class_clean = clean_tax(Class),
    Class_group = case_when(
      Class_clean == "Nitrososphaeria" ~ "Nitrososphaeria (AOA-affiliated)",
      Class_clean == "Thermoplasmata" ~ "Thermoplasmata",
      Phylum_clean == "Thermoplasmatota" & Class_clean != "Thermoplasmata" ~ "Thermoplasmatota (other/unclassified)",
      TRUE ~ "Other archaea"
    )
  ) %>%
  select(asv_id, Phylum_clean, Class_clean, Class_group)
env <- read.csv(env_path, check.names = FALSE) %>%
  filter(sample_id %in% retained_samples, depth < 60, site %in% site_levels, dna_type %in% pool_levels) %>%
  mutate(
    site = factor(site, levels = site_levels),
    site_label = factor(display_site(site), levels = site_display[site_levels]),
    dna_type = factor(dna_type, levels = pool_levels),
    depth_plot = factor(depths_cm_re, levels = rev(depth_levels_raw), labels = rev(depth_labels[depth_levels_raw])),
    site_pool = paste(site, dna_type, sep = "_")
  )

sample_cols <- retained_samples[retained_samples %in% colnames(asv_counts) & retained_samples %in% env$sample_id]
if (length(sample_cols) != length(retained_samples)) {
  stop("Expected all 95 retained 0-60 cm samples in the unrarefied table and metadata.")
}
asv_counts <- asv_counts[, sample_cols, drop = FALSE]
env <- env %>% filter(sample_id %in% sample_cols) %>% arrange(match(sample_id, sample_cols))

asv_long <- asv_counts %>%
  as.data.frame() %>%
  rownames_to_column("asv_id") %>%
  pivot_longer(-asv_id, names_to = "sample_id", values_to = "count") %>%
  left_join(taxa, by = "asv_id") %>%
  left_join(env %>% select(sample_id, site, site_label, dna_type, depth_plot, site_pool), by = "sample_id") %>%
  filter(!is.na(site), !is.na(Class_group)) %>%
  group_by(sample_id) %>%
  mutate(sample_rel = count / sum(count, na.rm = TRUE)) %>%
  ungroup()

sample_class <- asv_long %>%
  group_by(sample_id, site, site_label, dna_type, depth_plot, Class_group) %>%
  summarise(prop = sum(sample_rel, na.rm = TRUE), .groups = "drop") %>%
  complete(
    nesting(sample_id, site, site_label, dna_type, depth_plot),
    Class_group = composition_class_levels,
    fill = list(prop = 0)
  )

class_depth <- sample_class %>%
  group_by(site, site_label, dna_type, depth_plot, Class_group) %>%
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    signed_prop = if_else(dna_type == "iDNA", -prop, prop),
    Class_group = factor(Class_group, levels = composition_class_levels)
  )

# The mirrored direction is labelled directly in the upper-row facets. This is
# more immediate than a sentence-length subtitle and avoids repeating the same
# labels in all four site panels.
pool_side_labels <- tibble(
  site_label = factor(
    rep(c("AZ", "SG"), each = 2),
    levels = unname(site_display[site_levels])
  ),
  x = rep(c(-0.50, 0.50), times = 2),
  y = factor("0-5", levels = levels(class_depth$depth_plot)),
  label = rep(c("iDNA", "eDNA"), times = 2)
)

class_totals <- sample_class %>%
  group_by(Class_group) %>%
  summarise(mean_percent = 100 * mean(prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Class_group = factor(Class_group, levels = composition_class_levels)
  ) %>%
  arrange(Class_group) %>%
  mutate(percent = mean_percent)

message("Panel A class-level composition (mean within-sample relative abundance):")
purrr::pwalk(
  class_totals,
  function(Class_group, mean_percent, percent) {
    message(sprintf("  %s: %.1f%%", as.character(Class_group), percent))
  }
)

if (write_analysis_outputs) {
  write.csv(
    class_depth,
    file.path(analysis_dir, "Fig2_mirrored_class_depth_summary.csv"),
    row.names = FALSE
  )
}

p_composition <- ggplot(class_depth, aes(x = signed_prop, y = depth_plot, fill = Class_group)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.06, position = position_stack(reverse = TRUE)) +
  geom_vline(xintercept = 0, colour = "white", linewidth = 1.0) +
  geom_text(
    data = pool_side_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4.1,
    fontface = "bold",
    family = "Helvetica",
    colour = "grey18",
    vjust = -0.9
  ) +
  facet_wrap(~ site_label, ncol = 2) +
  scale_fill_manual(values = composition_class_cols, breaks = composition_class_levels, drop = FALSE) +
  scale_x_continuous(
    position = "top",
    breaks = seq(-1, 1, by = 0.25),
    labels = function(x) abs(x) * 100,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(
    title = "Archaeal class composition across depth featuring AOA-affiliated lineages",
    subtitle = NULL,
    x = "Relative abundance (%)",
    y = "Depth (cm)",
    fill = "Class"
  ) +
  theme_pub(base_size = 12.0, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 16.0),
    plot.subtitle = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14.0),
    axis.text.x = element_text(size = 12.0, angle = 45, hjust = 0),
    axis.text.y = element_text(size = 12.0),
    axis.title = element_text(size = 14.0, face = "bold"),
    panel.grid.major.x = element_line(colour = "grey87", linewidth = 0.22),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(1.8, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 12.0, face = "bold"),
    legend.text = element_text(size = 11.5),
    legend.key.size = unit(0.56, "cm")
  )

ann <- read.csv(ann_path, check.names = FALSE) %>%
  transmute(
    ASV,
    Phylum = clean_tax_lower_unclassified(Phylum),
    Class = clean_tax_lower_unclassified(Class),
    total_reads_60cm,
    read_share_60cm
  ) %>%
  mutate(
    Class_plot = case_when(
      Class %in% names(class_cols) ~ Class,
      TRUE ~ "unclassified"
    ),
    broad_lineage = case_when(
      Phylum == "Crenarchaeota" ~ "Crenarchaeota",
      Phylum %in% c("Thermoplasmatota", "Nanoarchaeota") ~ "Thermoplasmatota + Nanoarchaeota",
      TRUE ~ "Other"
    )
  )

legacy_tree <- read.tree(text = extract_newick(tree_path))
asv_keep <- intersect(legacy_tree$tip.label, ann$ASV)
if (length(asv_keep) == 0) {
  stop("No ASV tips from annotation table were found in the legacy tree.")
}
asv_tree <- keep.tip(legacy_tree, asv_keep)
tree_base <- ggtree(asv_tree, layout = "rectangular", size = 0.17, colour = "grey67")
tree_data <- tree_base$data
tip_data <- tree_data %>%
  filter(isTip) %>%
  left_join(ann, by = c("label" = "ASV"))

asv_rel <- asv_counts %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  filter(ASV %in% ann$ASV) %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = {
    total_count <- sum(count, na.rm = TRUE)
    if (total_count > 0) count / total_count else rep(0, n())
  }) %>%
  ungroup() %>%
  left_join(env %>% select(sample_id, site, site_label, dna_type, site_pool), by = "sample_id") %>%
  group_by(ASV, site, site_label, dna_type, site_pool) %>%
  summarise(mean_rel = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  complete(
    ASV = ann$ASV,
    site_pool = site_pool_levels,
    fill = list(mean_rel = 0)
  ) %>%
  mutate(
    site = factor(sub("_.*$", "", site_pool), levels = site_levels),
    site_label = factor(display_site(site), levels = site_display[site_levels]),
    dna_type = factor(sub("^.*_", "", site_pool), levels = pool_levels),
    site_pool = factor(site_pool, levels = site_pool_levels),
    site_pool_index = as.numeric(site_pool)
  ) %>%
  left_join(tip_data %>% select(ASV = label, y), by = "ASV") %>%
  filter(!is.na(y))

if (write_analysis_outputs) {
  write.csv(
    asv_rel,
    file.path(analysis_dir, "Fig2_tree_asv_site_pool_mean_relative_abundance.csv"),
    row.names = FALSE
  )
}

lineage_labels <- tip_data %>%
  filter(broad_lineage != "Other") %>%
  group_by(broad_lineage) %>%
  summarise(
    y = median(y, na.rm = TRUE),
    n_tips = n(),
    nitroso = sum(Class_plot == "Nitrososphaeria", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      broad_lineage == "Crenarchaeota" ~ paste0("Crenarchaeota: ", n_tips, " ASVs; ", nitroso, " Nitrososphaeria / AOA-affiliated"),
      broad_lineage == "Thermoplasmatota + Nanoarchaeota" ~ paste0("Thermoplasmatota + Nanoarchaeota: ", n_tips, " ASVs; mostly unresolved"),
      TRUE ~ broad_lineage
    )
  )

max_x <- max(tree_data$x, na.rm = TRUE)
max_y <- max(tree_data$y, na.rm = TRUE)
bubble_step <- max_x * 0.055
bubble_start <- max_x * 1.08
asv_rel <- asv_rel %>%
  mutate(bubble_x = bubble_start + (site_pool_index - 1) * bubble_step)
label_df <- tibble(
  site_pool = factor(site_pool_levels, levels = site_pool_levels),
  site_pool_index = seq_along(site_pool_levels),
  label = site_pool_labels,
  site = factor(sub("_.*$", "", site_pool_levels), levels = site_levels),
  x = bubble_start + (seq_along(site_pool_levels) - 1) * bubble_step,
  y = max_y + 3.0
)
lineage_labels <- lineage_labels %>%
  mutate(
    x = bubble_start + 8.25 * bubble_step,
    label = stringr::str_replace_all(label, "; ", "\n")
  )

p_tree <- tree_base +
  geom_point(
    data = tip_data,
    aes(x = x, y = y, colour = Class_plot),
    inherit.aes = FALSE,
    size = 1.1,
    alpha = 0.92
  ) +
  geom_point(
    data = asv_rel %>% filter(mean_rel > 0),
    aes(x = bubble_x, y = y, size = 100 * mean_rel, fill = site),
    inherit.aes = FALSE,
    shape = 21,
    colour = "grey30",
    stroke = 0.16,
    alpha = 0.86
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label, colour = site),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    size = 2.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_text(
    data = tibble(x = bubble_start + 3.5 * bubble_step, y = max_y + 12.5, label = "ASV mean relative abundance by site x DNA pool"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    fontface = "bold"
  ) +
  geom_segment(
    data = lineage_labels,
    aes(x = max_x * 0.72, xend = x - 0.02, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.24,
    colour = "grey35"
  ) +
  geom_text(
    data = lineage_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.9,
    fontface = "bold",
    colour = "grey15"
  ) +
  scale_colour_manual(
    values = c(class_cols, site_cols),
    breaks = names(class_cols),
    name = "Class",
    drop = FALSE
  ) +
  scale_fill_manual(values = site_cols, name = "Site", labels = site_display[names(site_cols)], drop = FALSE) +
  scale_size_area(
    max_size = 4.2,
    breaks = c(0.1, 0.5, 1, 2),
    limits = c(0, max(2, 100 * quantile(asv_rel$mean_rel, 0.98, na.rm = TRUE))),
    oob = scales::squish,
    name = "Mean relative\nabundance (%)\n(capped)"
  ) +
  coord_cartesian(
    xlim = c(0, bubble_start + 14.5 * bubble_step),
    ylim = c(0, max_y + 18),
    clip = "off"
  ) +
  labs(
    title = "Phylogenetic context of abundant archaeal ASVs with site-pool abundance",
    subtitle = "Bubbles show mean relative abundance from the full 0-60 cm ASV table; blank cells indicate no detected mean abundance."
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 2.0)),
    fill = guide_legend(order = 2, override.aes = list(size = 3.0)),
    size = guide_legend(order = 3)
  ) +
  theme_tree() +
  theme(
    text = element_text(family = "Helvetica", size = 8.4),
    plot.title = element_text(face = "bold", size = 10.2),
    plot.subtitle = element_text(size = 7.5, colour = "grey25", margin = margin(b = 2)),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 7.8, face = "bold"),
    legend.text = element_text(size = 7.0),
    legend.key.size = unit(0.36, "cm"),
    plot.margin = margin(4, 78, 4, 4)
  )

fig3 <- p_composition / p_tree +
  plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

out_analysis <- file.path(analysis_dir, "Fig2_reframed_phylum_phylo_context_round43.pdf")
out_figure <- file.path(figure_dir, "02_Fig2_reframed_phylum_phylo_context_round43.pdf")

if (!identical(Sys.getenv("FIG2_PANEL_A_NO_WRITE"), "1")) {
  ggsave('figures/Fig2_panelA_composition_round43.pdf', p_composition, width = 8.7, height = 5.6, device = cairo_pdf, bg = 'white')
  if (write_analysis_outputs) {
    ggsave(out_analysis, fig3, width = 8.7, height = 10.8, units = "in", device = cairo_pdf, bg = "white")
  }
  ggsave(out_figure, fig3, width = 8.7, height = 10.8, units = "in", device = cairo_pdf, bg = "white")
}

if (!identical(Sys.getenv("FIG2_PANEL_A_NO_WRITE"), "1")) {
  if (write_analysis_outputs) {
    message("Wrote: ", out_analysis)
  }
  message("Wrote: ", out_figure)
}
