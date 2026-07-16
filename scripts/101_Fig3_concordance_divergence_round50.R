#!/usr/bin/env Rscript

# Round50 Figure 3 candidate: the central concordance--divergence evidence chain.
#
# This script does not introduce a new statistical analysis. It reconstructs the
# primary two-dimensional symmetric Procrustes rotation used by the locked
# Round48 paired-pool audit, verifies it against the stored canonical result,
# and combines it with the stored within-site concordance and NA subsoil read-
# share outputs. Previous figures and manuscripts are never overwritten.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(vegan)
})

N_PERMUTATIONS <- 9999L
CONCORDANCE_SEED <- 20260605L
SITE_LEVELS <- c("AZ", "SG", "LC", "NA")
DNA_LEVELS <- c("iDNA", "eDNA")
DNA_COLOURS <- c(iDNA = "#2474B5", eDNA = "#E88919")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  project_root <- dirname(dirname(script_path))
} else {
  project_root <- normalizePath(".", mustWork = TRUE)
}

analysis_dir <- file.path(project_root, "analysis", "round50")
candidate_dir <- file.path(project_root, "figures", "candidates")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  rarefied_counts = file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt"),
  pair_map = file.path(
    project_root, "analysis", "paired_pool_similarity", "round48_complete_pair_map.csv"
  ),
  overall_concordance = file.path(
    project_root, "analysis", "paired_pool_similarity", "round48_paired_concordance_overall.csv"
  ),
  site_concordance = file.path(
    project_root, "analysis", "reproducibility", "round48_control_sensitivity_concordance_by_site.csv"
  ),
  na_subsoil = file.path(
    project_root, "analysis", "reproducibility", "round48_control_sensitivity_NA_subsoil_by_profile.csv"
  )
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

read_count_table <- function(path) {
  x <- as.matrix(read.delim(
    path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))
  storage.mode(x) <- "numeric"
  if (anyNA(x) || any(x < 0)) stop("Invalid count table: ", path)
  x
}

# ---------------------------------------------------------------------------
# Panel A: recreate the exact locked symmetric Procrustes rotation.
# ---------------------------------------------------------------------------
counts <- read_count_table(input_paths[["rarefied_counts"]])
pair_map <- read_csv(
  input_paths[["pair_map"]],
  na = character(),
  show_col_types = FALSE
) %>%
  mutate(site_display = factor(site_display, levels = SITE_LEVELS))

stopifnot(
  nrow(pair_map) == 39L,
  !anyDuplicated(pair_map$pair_id),
  all(pair_map$iDNA_sample %in% colnames(counts)),
  all(pair_map$eDNA_sample %in% colnames(counts)),
  all(colSums(counts) == 202)
)

paired_ids <- c(pair_map$iDNA_sample, pair_map$eDNA_sample)
paired_counts <- counts[, paired_ids, drop = FALSE]
paired_counts <- paired_counts[rowSums(paired_counts) > 0, , drop = FALSE]

i_mat <- t(paired_counts[, pair_map$iDNA_sample, drop = FALSE])
e_mat <- t(paired_counts[, pair_map$eDNA_sample, drop = FALSE])
rownames(i_mat) <- pair_map$pair_id
rownames(e_mat) <- pair_map$pair_id

i_hellinger <- decostand(i_mat, method = "hellinger")
e_hellinger <- decostand(e_mat, method = "hellinger")
i_pcoa <- wcmdscale(vegdist(i_hellinger, method = "bray"), k = 2, eig = TRUE)
e_pcoa <- wcmdscale(vegdist(e_hellinger, method = "bray"), k = 2, eig = TRUE)

set.seed(CONCORDANCE_SEED)
protest_fit <- protest(
  i_pcoa$points,
  e_pcoa$points,
  permutations = N_PERMUTATIONS
)

overall <- read_csv(input_paths[["overall_concordance"]], show_col_types = FALSE) %>%
  filter(scenario == "primary_unfiltered", method == "PROTEST Procrustes")
stopifnot(
  nrow(overall) == 1L,
  abs(unname(protest_fit$t0) - overall$statistic) < 1e-12,
  protest_fit$signif == overall$p_value,
  overall$n_complete_pairs == 39L,
  overall$n_asvs_analyzed == nrow(paired_counts)
)

i_scores <- as.data.frame(protest_fit$X) %>%
  tibble::rownames_to_column("pair_id") %>%
  setNames(c("pair_id", "axis1", "axis2")) %>%
  mutate(dna_pool = "iDNA")
e_scores <- as.data.frame(protest_fit$Yrot) %>%
  tibble::rownames_to_column("pair_id") %>%
  setNames(c("pair_id", "axis1", "axis2")) %>%
  mutate(dna_pool = "eDNA")

ordination_scores <- bind_rows(i_scores, e_scores) %>%
  left_join(
    pair_map %>% select(pair_id, site_internal, site_display, pit, depth_order, depth_cm),
    by = "pair_id"
  ) %>%
  mutate(
    dna_pool = factor(dna_pool, levels = DNA_LEVELS),
    site_display = factor(site_display, levels = SITE_LEVELS)
  )
stopifnot(nrow(ordination_scores) == 78L, !anyNA(ordination_scores$site_display))

pair_segments <- i_scores %>%
  select(pair_id, i_axis1 = axis1, i_axis2 = axis2) %>%
  left_join(
    e_scores %>% select(pair_id, e_axis1 = axis1, e_axis2 = axis2),
    by = "pair_id"
  ) %>%
  left_join(pair_map %>% select(pair_id, site_display), by = "pair_id")

site_centroids <- pair_segments %>%
  transmute(
    site_display,
    midpoint_axis1 = (i_axis1 + e_axis1) / 2,
    midpoint_axis2 = (i_axis2 + e_axis2) / 2
  ) %>%
  group_by(site_display) %>%
  summarise(
    axis1 = mean(midpoint_axis1),
    axis2 = mean(midpoint_axis2),
    .groups = "drop"
  )

p_a <- ggplot() +
  geom_segment(
    data = pair_segments,
    aes(x = i_axis1, y = i_axis2, xend = e_axis1, yend = e_axis2),
    linewidth = 0.42,
    colour = "grey62",
    alpha = 0.72,
    lineend = "round"
  ) +
  geom_point(
    data = ordination_scores,
    aes(axis1, axis2, fill = dna_pool),
    shape = 21,
    size = 2.25,
    stroke = 0.38,
    colour = "grey20",
    alpha = 0.92
  ) +
  geom_label_repel(
    data = site_centroids,
    aes(axis1, axis2, label = site_display),
    seed = 20260716,
    size = 3.25,
    family = "Helvetica",
    fontface = "bold",
    colour = "grey20",
    fill = scales::alpha("white", 0.88),
    label.size = 0.2,
    label.padding = grid::unit(0.12, "lines"),
    box.padding = 0.35,
    point.padding = 0.7,
    min.segment.length = 0,
    segment.colour = "grey60",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = DNA_COLOURS,
    breaks = DNA_LEVELS,
    labels = c("iDNA (intact-cell-enriched)", "eDNA (extracellular)"),
    name = "DNA pool"
  ) +
  coord_equal() +
  labs(
    title = "Paired iDNA and eDNA communities align",
    subtitle = sprintf(
      "PROTEST r = %.3f, p = %.4f; 39 matched horizons",
      overall$statistic,
      overall$p_value
    ),
    x = "Procrustes axis 1",
    y = "Procrustes axis 2"
  )

# ---------------------------------------------------------------------------
# Panel B: site-level primary concordance and conservative control sensitivity.
# ---------------------------------------------------------------------------
site_raw <- read_csv(
  input_paths[["site_concordance"]],
  na = character(),
  show_col_types = FALSE
) %>%
  filter(scenario %in% c("primary_unfiltered", "control_filtered")) %>%
  mutate(
    site_display = factor(site_display, levels = rev(SITE_LEVELS)),
    scenario = factor(
      scenario,
      levels = c("control_filtered", "primary_unfiltered")
    )
  )

stopifnot(
  nrow(site_raw) == 8L,
  all(table(site_raw$site_display) == 2L),
  all(site_raw$n_complete_pairs %in% c(8L, 9L, 11L))
)

site_connections <- site_raw %>%
  select(site_display, scenario, protest_r) %>%
  pivot_wider(names_from = scenario, values_from = protest_r)

site_n <- site_raw %>%
  filter(scenario == "primary_unfiltered") %>%
  select(site_display, n_complete_pairs)

p_b <- ggplot() +
  geom_segment(
    data = site_connections,
    aes(
      x = control_filtered,
      xend = primary_unfiltered,
      y = site_display,
      yend = site_display
    ),
    linewidth = 0.75,
    colour = "grey73",
    lineend = "round"
  ) +
  geom_point(
    data = filter(site_raw, scenario == "control_filtered"),
    aes(protest_r, site_display),
    shape = 21,
    size = 3.2,
    stroke = 0.65,
    fill = "white",
    colour = "#6E7781"
  ) +
  geom_point(
    data = filter(site_raw, scenario == "primary_unfiltered"),
    aes(protest_r, site_display),
    shape = 21,
    size = 3.3,
    stroke = 0.5,
    fill = "#2F5D7C",
    colour = "white"
  ) +
  geom_text(
    data = site_n,
    aes(x = 1.008, y = site_display, label = paste0("n = ", n_complete_pairs)),
    hjust = 0,
    size = 3.05,
    family = "Helvetica",
    colour = "grey30"
  ) +
  scale_x_continuous(
    limits = c(0.80, 1.045),
    breaks = c(0.80, 0.85, 0.90, 0.95, 1.00),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "High concordance persists within sites",
    subtitle = "Filled = primary; open = control-filter sensitivity",
    x = "Within-site Procrustes r",
    y = NULL
  )

# ---------------------------------------------------------------------------
# Panel C: the local quantitative divergence, shown as true matched profiles.
# ---------------------------------------------------------------------------
na_subsoil <- read_csv(
  input_paths[["na_subsoil"]],
  na = character(),
  show_col_types = FALSE
) %>%
  filter(
    scenario == "primary_unfiltered",
    site == "NA",
    depth_cm %in% c("20-40", "40-60")
  ) %>%
  mutate(
    depth_cm = factor(depth_cm, levels = c("20-40", "40-60")),
    dna_pool = factor(dna_pool, levels = DNA_LEVELS)
  )

stopifnot(
  nrow(na_subsoil) == 12L,
  n_distinct(na_subsoil$pair_id) == 6L,
  all(table(na_subsoil$pair_id) == 2L)
)

paired_direction <- na_subsoil %>%
  select(pair_id, pit, depth_cm, dna_pool, archaeal_percent) %>%
  pivot_wider(names_from = dna_pool, values_from = archaeal_percent) %>%
  mutate(eDNA_higher = eDNA > iDNA)
stopifnot(all(paired_direction$eDNA_higher))

na_means <- na_subsoil %>%
  group_by(depth_cm, dna_pool) %>%
  summarise(mean_percent = mean(archaeal_percent), .groups = "drop")

mean_labels <- na_means %>%
  pivot_wider(names_from = dna_pool, values_from = mean_percent) %>%
  mutate(
    label = sprintf("mean %.1f%%  →  %.1f%%", iDNA, eDNA),
    x = 1.5,
    y = 17.0
  )

p_c <- ggplot(
  na_subsoil,
  aes(dna_pool, archaeal_percent, group = pair_id)
  ) +
  geom_line(
    colour = "grey68",
    linewidth = 0.7,
    alpha = 0.9,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = dna_pool),
    shape = 21,
    size = 2.65,
    stroke = 0.42,
    colour = "grey20",
    alpha = 0.96,
    show.legend = FALSE
  ) +
  geom_line(
    data = na_means,
    aes(dna_pool, mean_percent, group = depth_cm),
    inherit.aes = FALSE,
    colour = "grey15",
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_point(
    data = na_means,
    aes(dna_pool, mean_percent, fill = dna_pool),
    inherit.aes = FALSE,
    shape = 23,
    size = 4.0,
    stroke = 0.75,
    colour = "grey10",
    show.legend = FALSE
  ) +
  geom_text(
    data = mean_labels,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    family = "Helvetica",
    fontface = "bold",
    size = 3.25,
    colour = "grey20"
  ) +
  facet_wrap(
    ~depth_cm,
    nrow = 1,
    labeller = as_labeller(c("20-40" = "20–40 cm", "40-60" = "40–60 cm"))
  ) +
  scale_fill_manual(values = DNA_COLOURS, breaks = DNA_LEVELS) +
  scale_y_continuous(
    limits = c(0, 18),
    breaks = c(0, 4, 8, 12, 16),
    labels = scales::label_number(suffix = "%", accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Archaeal read share diverges in humid subsoil",
    subtitle = "Nahuelbuta (NA); grey lines link matched profile–depth pairs",
    x = "Operational DNA pool",
    y = "Archaeal reads\n(% of total prokaryotic reads)"
  )

# ---------------------------------------------------------------------------
# House-style composition: native panel geometry, optical gutters, one legend.
# ---------------------------------------------------------------------------
base_theme <- theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 12.6, hjust = 0),
    plot.subtitle = element_text(
      size = 9.2,
      colour = "grey28",
      margin = margin(t = 1, b = 5)
    ),
    axis.title = element_text(size = 10.8),
    axis.text = element_text(size = 9.7, colour = "grey20"),
    axis.line = element_line(linewidth = 0.45, colour = "grey20"),
    axis.ticks = element_line(linewidth = 0.4, colour = "grey30"),
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey70",
      linewidth = 0.35
    ),
    strip.text = element_text(face = "bold", size = 10.1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 9.8),
    legend.text = element_text(size = 9.3),
    legend.key.width = grid::unit(1.05, "cm"),
    plot.margin = margin(t = 7, r = 8, b = 6, l = 7)
  )

p_a <- p_a + base_theme +
  theme(
    panel.border = element_rect(fill = NA, colour = "grey25", linewidth = 0.45),
    axis.line = element_blank()
  )

p_b <- p_b + base_theme +
  theme(
    panel.grid.major.x = element_line(colour = "grey89", linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 7, r = 15, b = 6, l = 12)
  )

p_c <- p_c + base_theme +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    panel.spacing.x = grid::unit(12, "pt"),
    plot.margin = margin(t = 12, r = 8, b = 6, l = 7)
  )

top_row <- p_a + p_b + plot_layout(widths = c(1.16, 0.84))
figure <- top_row / p_c +
  plot_layout(heights = c(1.12, 0.88), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 20,
        family = "Helvetica",
        colour = "black"
      ),
      plot.tag.position = c(0.002, 0.995)
    )
  ) &
  theme(legend.position = "bottom")

output_pdf <- file.path(
  candidate_dir,
  "Fig3_concordance_divergence_round50_candidate_v2.pdf"
)
ggsave(
  output_pdf,
  figure,
  width = 10.0,
  height = 8.25,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

write_csv(
  ordination_scores %>%
    arrange(site_display, pit, depth_order, dna_pool),
  file.path(analysis_dir, "Fig3_round50_procrustes_coordinates.csv")
)
write_csv(
  site_raw %>% arrange(site_display, scenario),
  file.path(analysis_dir, "Fig3_round50_site_concordance.csv")
)
write_csv(
  na_subsoil %>% arrange(depth_cm, pit, dna_pool),
  file.path(analysis_dir, "Fig3_round50_NA_subsoil_paired_read_share.csv")
)
write_csv(
  paired_direction %>% arrange(depth_cm, pit),
  file.path(analysis_dir, "Fig3_round50_NA_subsoil_pair_directions.csv")
)
write_csv(
  tibble(
    input = names(input_paths),
    project_relative_path = sub(
      paste0("^", normalizePath(project_root, winslash = "/"), "/"),
      "",
      normalizePath(input_paths, winslash = "/")
    ),
    bytes = file.info(input_paths)$size,
    md5 = unname(tools::md5sum(input_paths))
  ),
  file.path(analysis_dir, "Fig3_round50_input_provenance.csv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(analysis_dir, "Fig3_round50_R_session_info.txt")
)

message("Wrote: ", output_pdf)
message(
  "Verified locked PROTEST r = ", sprintf("%.6f", protest_fit$t0),
  "; all six NA subsoil pairs have eDNA > iDNA."
)
