#!/usr/bin/env Rscript

# Figure 3 (supporting alpha diversity), Round47 wording refinement of the
# validated Round43 compact-annotation layout.
# Replaces script/45_Fig2_supporting_shannon.R, whose
# bacterial median trends were drawn on an archaeal y-axis (archaeal Shannon
# ~1.4-3.7 vs bacterial Shannon ~6.0) and were therefore clipped off-panel
# while the legend and subtitle still promised them.
#
# Panel A keeps the useful part of script/45: paired iDNA/eDNA archaeal Shannon
#   distributions with BH-adjusted paired Wilcoxon tests per site.
# Panel B restores the cross-domain contrast from script/23 using the ONLY
#   legitimate comparison: site means rescaled to 0-1 within each domain and
#   metric, so the two rarefaction depths never share a raw axis.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
house_style_run <- identical(Sys.getenv("ROUND48_HOUSESTYLE_OUTPUT"), "1")
figure_dir <- file.path(project_root, "figures")
candidate_dir <- file.path(figure_dir, "candidates")
provisional_dir <- file.path(figure_dir, "provisional")
analysis_dir <- if (house_style_run) {
  file.path(tempdir(), "chile_archaea_round48_housestyle", "Fig3")
} else {
  file.path(project_root, "analysis", "reframed_figures")
}

archaea_alpha_path <- file.path(project_root, "analysis", "submission_statistics", "alpha_diversity_60cm.csv")
bacteria_root <- file.path(project_root, "data", "companion_bacteria")
bacteria_alpha_path <- file.path(bacteria_root, "processed", "bac_200_9435rarefild", "alpha_diversity.txt")
bacteria_env_path <- file.path(bacteria_root, "processed", "bac_200_9435rarefild", "rare_env_200.txt")

required <- c(archaea_alpha_path, bacteria_alpha_path, bacteria_env_path)
if (any(!file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = "; "))
for (d in c(figure_dir, candidate_dir, provisional_dir, analysis_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

site_levels <- c("AZ", "SG", "LC", "NB")
site_display <- c(AZ = "AZ", SG = "SG", LC = "LC", NB = "NA")
metric_levels <- c("Shannon", "Observed_ASVs", "Evenness")
metric_labels <- c("Shannon diversity", "Observed ASV richness", "Pielou evenness")
pool_cols <- c(iDNA = "dodgerblue2", eDNA = "darkorange")
domain_cols <- c(Archaea = "#7B3FA0", Bacteria = "#4C9A2A")

lab_site <- function(x) factor(unname(site_display[as.character(x)]), levels = unname(site_display[site_levels]))

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

# ---------------------------------------------------------------- data
archaea_alpha <- read.csv(archaea_alpha_path, check.names = FALSE) %>%
  filter(site %in% site_levels, dna_type %in% c("iDNA", "eDNA")) %>%
  mutate(
    Evenness = Shannon / log(Observed_ASVs),
    domain = "Archaea",
    site = factor(site, levels = site_levels),
    site_label = lab_site(site),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA"))
  )

bacteria_env <- read.delim(bacteria_env_path, check.names = FALSE)
bacteria_alpha <- read.delim(bacteria_alpha_path, check.names = FALSE) %>%
  left_join(bacteria_env %>% select(sample_id, site, dna_type, depth), by = "sample_id") %>%
  # review MC1: restrict bacteria to the exact sample IDs retained in the
  # archaeal dataset, so the cross-domain contrast uses matched samples.
  filter(sample_id %in% archaea_alpha$sample_id,
         depth < 60, site %in% site_levels, dna_type %in% c("iDNA", "eDNA")) %>%
  mutate(
    Observed_ASVs = observed_species,
    Evenness = Pielous,
    domain = "Bacteria",
    site = factor(site, levels = site_levels),
    site_label = lab_site(site),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA"))
  )

# ------------------------------------------------- panel A: paired archaeal Shannon
pair_lookup <- read.csv(file.path(project_root, "data", "env_2024.csv"), check.names = FALSE) %>%
  transmute(sample_id, pair_id = interaction(site, pit, depth_order, drop = TRUE))

paired_tests <- archaea_alpha %>%
  left_join(pair_lookup, by = "sample_id") %>%
  select(site, site_label, dna_type, pair_id, Shannon) %>%
  pivot_wider(names_from = dna_type, values_from = Shannon) %>%
  filter(!is.na(iDNA), !is.na(eDNA)) %>%
  group_by(site, site_label) %>%
  summarise(n_pairs = n(), p = wilcox.test(iDNA, eDNA, paired = TRUE, exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(p_bh = p.adjust(p, method = "BH"), label = if_else(p_bh < 0.05, "*", "ns"))

archaea_medians <- archaea_alpha %>%
  group_by(site_label, dna_type) %>%
  summarise(archaea_median = median(Shannon, na.rm = TRUE), .groups = "drop")

rng <- range(archaea_alpha$Shannon, na.rm = TRUE)
span <- diff(rng)

# No per-site significance brackets: all four paired tests are ns, so brackets
# add no information. Keep only the test, correction, and adjusted p value in
# the panel; matching details remain in the editable figure legend.
sub_a <- if (all(paired_tests$label == "ns")) {
  bquote(paste("Paired Wilcoxon; BH-adjusted ", italic(p), " = ",
               .(sprintf("%.2f", max(paired_tests$p_bh))), " at all sites"))
} else {
  bquote(paste("Paired Wilcoxon; BH-adjusted ", italic(p), " < 0.05 at ",
               .(paste(paired_tests$site_label[paired_tests$label == "*"], collapse = ", "))))
}

pA <- ggplot(archaea_alpha, aes(site_label, Shannon)) +
  geom_violin(aes(fill = dna_type), position = position_dodge(width = 0.72), width = 0.72,
              alpha = 0.36, linewidth = 0.27, colour = "grey35", trim = TRUE) +
  geom_boxplot(aes(fill = dna_type), position = position_dodge(width = 0.72), width = 0.18,
               outlier.shape = NA, alpha = 0.72, linewidth = 0.28, colour = "grey25", show.legend = FALSE) +
  geom_point(aes(colour = dna_type),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.72, seed = 11),
             size = 1.05, alpha = 0.40, stroke = 0, show.legend = FALSE) +
  geom_line(data = archaea_medians, aes(site_label, archaea_median, group = dna_type, colour = dna_type),
            linewidth = 0.75, inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(data = archaea_medians, aes(site_label, archaea_median, colour = dna_type),
             size = 2.2, inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = pool_cols, name = "Archaeal DNA fraction") +
  scale_colour_manual(values = pool_cols, guide = "none") +
  coord_cartesian(ylim = c(rng[1] - 0.04 * span, rng[2] + 0.06 * span)) +
  labs(
    title = "Archaeal diversity peaks at humid Nahuelbuta in both DNA fractions",
    subtitle = sub_a,
    x = NULL, y = "Shannon diversity"
  )

# ------------------------------------------- panel B: cross-domain, within-domain scaling
domain_site_summary <- bind_rows(
  archaea_alpha %>% select(sample_id, domain, site_label, dna_type, all_of(metric_levels)),
  bacteria_alpha %>% select(sample_id, domain, site_label, dna_type, all_of(metric_levels))
) %>%
  filter(dna_type == "iDNA") %>%
  pivot_longer(all_of(metric_levels), names_to = "metric", values_to = "value") %>%
  group_by(domain, site_label, metric) %>%
  summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(domain, metric) %>%
  mutate(scaled_mean = rescale01(mean)) %>%
  ungroup() %>%
  mutate(
    domain = factor(domain, levels = c("Archaea", "Bacteria")),
    metric = factor(metric, levels = metric_levels, labels = metric_labels)
  )

peak_labels <- domain_site_summary %>%
  group_by(domain, metric) %>%
  slice_max(scaled_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = paste0("peak: ", site_label), y = 1.10)

pB <- ggplot(domain_site_summary, aes(site_label, scaled_mean, colour = domain, group = domain, shape = domain)) +
  geom_line(linewidth = 0.72) +
  geom_point(size = 2.1) +
  geom_text(data = peak_labels, aes(site_label, y, label = label, colour = domain),
            inherit.aes = FALSE, size = 2.9, show.legend = FALSE) +
  facet_wrap(~metric, nrow = 1) +
  scale_colour_manual(values = domain_cols, name = NULL) +
  scale_shape_manual(values = c(Archaea = 16, Bacteria = 17), name = NULL) +
  coord_cartesian(ylim = c(-0.06, 1.22)) +
  labs(
    title = "Archaea and bacteria show different site-level diversity trajectories",
    x = "Site along arid-to-humid gradient", y = "Scaled site mean",
    caption = paste0("iDNA only; sample-matched site means scaled to 0-1 within each domain and metric ",
                     "because rarefaction depths differ.")
  )

# ---------------------------------------------------------------- compose
base <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.8, colour = "grey25", margin = margin(b = 4)),
    plot.caption = element_text(size = 8.8, colour = "grey25", hjust = 0,
                                margin = margin(t = 5, b = 0)),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 9.5),
    plot.margin = margin(t = 6, r = 8, b = 4, l = 6)
  )

fig <- (pA & base) / (pB & base) +
  plot_layout(heights = c(1, 0.95)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold", size = 18, family = "Helvetica", colour = "black"
      )
    )
  )

out_pdf <- if (house_style_run) {
  file.path(candidate_dir, "Fig3_round48_housestyle.pdf")
} else {
  file.path(figure_dir, "Fig3_crossdomain_diversity_round47.pdf")
}
ggsave(out_pdf, fig, width = 9.0, height = 8.2, units = "in", device = grDevices::cairo_pdf, bg = "white")
if (!house_style_run) {
  file.copy(out_pdf, file.path(provisional_dir, basename(out_pdf)), overwrite = TRUE)
}
write_csv(paired_tests, file.path(analysis_dir, "Fig3_shannon_paired_site_tests_BH_round47.csv"))

# guard: the bug that broke script/45 was a cross-domain raw-axis overlay
stopifnot(all(domain_site_summary$scaled_mean >= 0, domain_site_summary$scaled_mean <= 1))
message("Wrote: ", out_pdf)
message("Peak sites: ", paste(sprintf("%s/%s=%s", peak_labels$domain, peak_labels$metric, peak_labels$site_label), collapse = "  "))
