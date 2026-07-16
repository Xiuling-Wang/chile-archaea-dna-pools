#!/usr/bin/env Rscript

# Supplementary QC figure: rarefaction curve plus a 500-read alpha-diversity
# sensitivity check for archaeal 0-60 cm samples.

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
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
  file.path(tempdir(), "chile_archaea_round48_housestyle", "FigS1")
} else {
  file.path(project_root, "analysis", "reframed_figures")
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

raw_path <- file.path(project_root, "data", "file_A.txt")
rare202_path <- file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt")
env_path <- file.path(project_root, "data", "env_2024.csv")

missing_inputs <- c(raw_path, rare202_path, env_path)[!file.exists(c(raw_path, rare202_path, env_path))]
if (length(missing_inputs) > 0) {
  stop("Missing input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

site_levels <- c("AZ", "SG", "LC", "NB")
site_display <- c(AZ = "AZ", SG = "SG", LC = "LC", NB = "NA")
dna_levels <- c("iDNA", "eDNA")
metric_levels <- c("Shannon", "Observed_ASVs", "Evenness")
metric_labels <- c(
  Shannon = "Shannon diversity",
  Observed_ASVs = "Observed ASV richness",
  Evenness = "Pielou evenness"
)
site_cols <- c(AZ = "#D95F02", SG = "#1B9E77", LC = "#7570B3", NB = "#66A61E")
pool_cols <- c(iDNA = "dodgerblue2", eDNA = "darkorange")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

display_site <- function(x) unname(site_display[as.character(x)])

raw <- read.delim(raw_path, header = TRUE, row.names = 1, check.names = FALSE)
rare202 <- read.delim(rare202_path, header = TRUE, row.names = 1, check.names = FALSE)
env <- read.csv(env_path, header = TRUE, check.names = FALSE) %>%
  filter(depth < 60, site %in% site_levels, dna_type %in% dna_levels)

archaea_raw <- raw %>%
  filter(Kingdom == "Archaea")
raw_counts <- archaea_raw[, setdiff(colnames(archaea_raw), tax_cols), drop = FALSE]

matched_samples <- intersect(env$sample_id, colnames(raw_counts))
raw_counts <- raw_counts[, matched_samples, drop = FALSE]
raw_counts <- raw_counts[rowSums(raw_counts, na.rm = TRUE) > 0, , drop = FALSE]
raw_depth <- colSums(raw_counts, na.rm = TRUE)

rare202_samples <- intersect(colnames(rare202), colnames(raw_counts))
n_202 <- length(rare202_samples)

otu202 <- t(raw_counts[, rare202_samples, drop = FALSE])
meta202 <- tibble(sample_id = rownames(otu202)) %>%
  left_join(env, by = "sample_id") %>%
  mutate(
    site = factor(site, levels = site_levels),
    site_label = factor(display_site(site), levels = site_display[site_levels]),
    dna_type = factor(dna_type, levels = dna_levels)
  )

set.seed(42)
rarefaction_curve <- rarecurve(otu202, step = 10, sample = 202, tidy = TRUE) %>%
  rename(sample_id = Site, n_reads = Sample, n_asvs = Species) %>%
  left_join(meta202 %>% select(sample_id, site, site_label, dna_type), by = "sample_id")

p_rarefaction <- ggplot(
  rarefaction_curve %>% filter(n_reads <= 1500),
  aes(x = n_reads, y = n_asvs, group = sample_id, colour = site, linetype = dna_type)
) +
  geom_line(alpha = 0.60, linewidth = 0.42) +
  geom_vline(xintercept = 202, colour = "grey30", linetype = "dashed", linewidth = 0.55) +
  annotate(
    "text",
    x = 220,
    y = 1,
    label = "202 reads\nbaseline",
    hjust = 0,
    vjust = 0,
    size = 2.35,
    colour = "grey30",
    family = "Helvetica"
  ) +
  scale_colour_manual(values = site_cols, name = "Site", labels = site_display[names(site_cols)], drop = FALSE) +
  scale_linetype_manual(values = c(iDNA = "solid", eDNA = "dashed"), name = "DNA pool", drop = FALSE) +
  scale_x_continuous(limits = c(0, 1500), labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Archaeal ASV rarefaction curves",
    subtitle = "Curves show the 95 samples retained for the 202-read baseline.",
    x = "Sequencing depth (reads)",
    y = "Observed ASVs"
  ) +
  theme_pub(base_size = 9.8, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 11.2),
    plot.subtitle = element_text(size = 8.4, colour = "grey25", margin = margin(b = 3)),
    axis.text = element_text(size = 8.8),
    axis.title = element_text(size = 10.0, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 9.0, face = "bold"),
    legend.text = element_text(size = 8.3)
  )

samples500 <- names(raw_depth)[raw_depth >= 500]
n_500 <- length(samples500)
otu500 <- t(raw_counts[, samples500, drop = FALSE])

set.seed(202)
otu500_rare <- rrarefy(otu500, sample = 500)

alpha500 <- tibble(
  sample_id = rownames(otu500_rare),
  Shannon = diversity(otu500_rare, index = "shannon"),
  Observed_ASVs = specnumber(otu500_rare)
) %>%
  mutate(Evenness = Shannon / log(Observed_ASVs)) %>%
  left_join(env, by = "sample_id") %>%
  filter(site %in% site_levels, dna_type %in% dna_levels) %>%
  mutate(
    site = factor(site, levels = site_levels),
    site_label = factor(display_site(site), levels = site_display[site_levels]),
    dna_type = factor(dna_type, levels = dna_levels)
  )

alpha500_long <- alpha500 %>%
  select(sample_id, site, site_label, dna_type, all_of(metric_levels)) %>%
  pivot_longer(all_of(metric_levels), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_levels, labels = metric_labels[metric_levels]))

alpha500_medians <- alpha500_long %>%
  group_by(metric, site_label, dna_type) %>%
  summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")

write.csv(
  alpha500,
  file.path(analysis_dir, "SI_QC_alpha_diversity_500rare.csv"),
  row.names = FALSE
)

p_sensitivity <- ggplot(alpha500_long, aes(x = site_label, y = value)) +
  geom_violin(
    aes(fill = dna_type),
    position = position_dodge(width = 0.72),
    width = 0.72,
    alpha = 0.38,
    linewidth = 0.25,
    colour = "grey35",
    trim = TRUE
  ) +
  geom_boxplot(
    aes(fill = dna_type),
    position = position_dodge(width = 0.72),
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.26,
    colour = "grey25",
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = dna_type),
    position = position_jitterdodge(jitter.width = 0.09, dodge.width = 0.72, seed = 11),
    size = 1.05,
    alpha = 0.42,
    stroke = 0,
    show.legend = FALSE
  ) +
  geom_line(
    data = alpha500_medians,
    aes(x = site_label, y = median_value, group = dna_type, colour = dna_type),
    linewidth = 0.72,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = alpha500_medians,
    aes(x = site_label, y = median_value, colour = dna_type),
    size = 2.05,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pool_cols, name = "DNA pool", drop = FALSE) +
  scale_colour_manual(values = pool_cols, name = NULL, drop = FALSE) +
  labs(
    title = "Alpha-diversity sensitivity after rarefaction to 500 reads",
    subtitle = sprintf("Retained samples: 500-read n = %d; 202-read baseline n = %d.", n_500, n_202),
    x = "Site along arid-to-humid gradient",
    y = "Archaeal alpha diversity"
  ) +
  theme_pub(base_size = 9.8, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 11.2),
    plot.subtitle = element_text(size = 8.4, colour = "grey25", margin = margin(b = 3)),
    strip.text = element_text(size = 9.6),
    axis.text = element_text(size = 8.8),
    axis.title = element_text(size = 10.0, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 9.0, face = "bold")
  )

qc_fig <- p_rarefaction / p_sensitivity +
  plot_layout(heights = c(0.95, 1.05), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold", size = 20, family = "Helvetica", colour = "black"
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
  ) &
  theme(legend.position = "bottom")

out_pdf <- if (house_style_run) {
  file.path(candidate_dir, "FigS1_round48_housestyle.pdf")
} else {
  file.path(figure_dir, "SI_QC_rarefaction_plus_500sensitivity.pdf")
}
ggsave(out_pdf, qc_fig, width = 11.6, height = 8.4, units = "in", device = cairo_pdf, bg = "white")

message("Retained samples at 202 reads: ", n_202)
message("Retained samples at 500 reads: ", n_500)
message("Wrote: ", out_pdf)
