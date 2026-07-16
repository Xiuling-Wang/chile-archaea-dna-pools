suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
})

try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

source(file.path("scripts", "00_config.R"))

report_dir <- ensure_dir(file.path(REPORT_DIR, "paired_pool_similarity"))
supp_fig_dir <- ensure_dir(file.path(PAPER_FIG_DIR, "Supplementary_Figures"))

asv_path <- file.path(DATA_DIR, "arc_60cm_202rare", "asv_202.txt")
env_path <- file.path(DATA_DIR, "env_2024.csv")

asv_raw <- read_tsv_rownames(asv_path)
env_raw <- read.csv(env_path, header = TRUE, check.names = FALSE)

env <- env_raw %>%
  filter(depth < 60, sample_id %in% colnames(asv_raw)) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    site_label = factor(if_else(site == "NB", "NA", as.character(site)),
                        levels = SITE_LABELS),
    dna_type = factor(dna_type, levels = DNA_LEVELS),
    depth_order = as.numeric(depth_order),
    depth_label_clean = factor(
      str_replace_all(depths_cm, "_", "-"),
      levels = DEPTH_60_LABELS
    ),
    pair_id = paste(site, pit, depths_cm, sep = "_")
  )

asv <- asv_raw[, env$sample_id, drop = FALSE]
asv <- asv[rowSums(asv) > 0, colSums(asv) > 0, drop = FALSE]
env <- env %>%
  filter(sample_id %in% colnames(asv)) %>%
  arrange(match(sample_id, colnames(asv)))
asv <- asv[, env$sample_id, drop = FALSE]
stopifnot(identical(colnames(asv), env$sample_id))

pair_table <- env %>%
  count(pair_id, site, site_label, pit, depth_order, depth_label_clean, dna_type) %>%
  pivot_wider(names_from = dna_type, values_from = n, values_fill = 0)

complete_pairs <- pair_table %>%
  filter(iDNA == 1, eDNA == 1) %>%
  select(pair_id, site, site_label, pit, depth_order, depth_label_clean)

paired_env <- env %>%
  semi_join(complete_pairs, by = "pair_id") %>%
  arrange(pair_id, dna_type)

paired_asv <- asv[, paired_env$sample_id, drop = FALSE]
paired_asv <- paired_asv[rowSums(paired_asv) > 0, , drop = FALSE]

sample_matrix <- t(paired_asv)
sample_matrix_pa <- decostand(sample_matrix, method = "pa")

calc_pair_metrics <- function(pair) {
  ids <- paired_env %>%
    filter(pair_id == pair) %>%
    arrange(dna_type) %>%
    pull(sample_id)

  mat <- sample_matrix[ids, , drop = FALSE]
  mat_pa <- sample_matrix_pa[ids, , drop = FALSE]
  rel <- decostand(mat, method = "total")

  bray <- as.numeric(vegdist(mat, method = "bray"))
  jaccard <- as.numeric(vegdist(mat_pa, method = "jaccard", binary = TRUE))
  shared_asv <- sum(mat_pa[1, ] > 0 & mat_pa[2, ] > 0)
  union_asv <- sum(mat_pa[1, ] > 0 | mat_pa[2, ] > 0)
  spearman <- suppressWarnings(cor(rel[1, ], rel[2, ], method = "spearman"))
  cosine <- sum(rel[1, ] * rel[2, ]) /
    sqrt(sum(rel[1, ]^2) * sum(rel[2, ]^2))

  tibble(
    pair_id = pair,
    iDNA_sample = ids[paired_env$dna_type[match(ids, paired_env$sample_id)] == "iDNA"],
    eDNA_sample = ids[paired_env$dna_type[match(ids, paired_env$sample_id)] == "eDNA"],
    bray_dissimilarity = bray,
    bray_similarity = 1 - bray,
    jaccard_dissimilarity = jaccard,
    jaccard_similarity = 1 - jaccard,
    shared_asv = shared_asv,
    union_asv = union_asv,
    shared_asv_fraction = shared_asv / union_asv,
    spearman_similarity = spearman,
    cosine_similarity = cosine
  )
}

paired_metrics <- map_dfr(complete_pairs$pair_id, calc_pair_metrics) %>%
  left_join(complete_pairs, by = "pair_id") %>%
  arrange(site, pit, depth_order)

write_csv(paired_metrics, file.path(report_dir, "paired_iDNA_eDNA_similarity_by_sample.csv"))

site_depth_summary <- paired_metrics %>%
  group_by(site, site_label, depth_order, depth_label_clean) %>%
  summarise(
    n_pairs = n(),
    bray_mean = mean(bray_dissimilarity),
    bray_sd = sd(bray_dissimilarity),
    jaccard_mean = mean(jaccard_dissimilarity),
    jaccard_sd = sd(jaccard_dissimilarity),
    shared_asv_fraction_mean = mean(shared_asv_fraction),
    spearman_mean = mean(spearman_similarity, na.rm = TRUE),
    cosine_mean = mean(cosine_similarity, na.rm = TRUE),
    .groups = "drop"
  )

site_summary <- paired_metrics %>%
  group_by(site, site_label) %>%
  summarise(
    n_pairs = n(),
    bray_mean = mean(bray_dissimilarity),
    bray_sd = sd(bray_dissimilarity),
    bray_median = median(bray_dissimilarity),
    jaccard_mean = mean(jaccard_dissimilarity),
    jaccard_sd = sd(jaccard_dissimilarity),
    jaccard_median = median(jaccard_dissimilarity),
    shared_asv_fraction_mean = mean(shared_asv_fraction),
    spearman_mean = mean(spearman_similarity, na.rm = TRUE),
    cosine_mean = mean(cosine_similarity, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(site_depth_summary, file.path(report_dir, "paired_similarity_summary_by_site_depth.csv"))
write_csv(site_summary, file.path(report_dir, "paired_similarity_summary_by_site.csv"))

same_site_depth_unmatched <- paired_env %>%
  select(pair_id, sample_id, site, site_label, depth_order, depth_label_clean, dna_type) %>%
  inner_join(
    paired_env %>%
      select(pair_id, sample_id, site, site_label, depth_order, depth_label_clean, dna_type),
    by = c("site", "site_label", "depth_order", "depth_label_clean"),
    suffix = c("_i", "_e"),
    relationship = "many-to-many"
  ) %>%
  filter(dna_type_i == "iDNA", dna_type_e == "eDNA", pair_id_i != pair_id_e)

calc_unmatched <- function(i_sample, e_sample) {
  mat <- sample_matrix[c(i_sample, e_sample), , drop = FALSE]
  mat_pa <- sample_matrix_pa[c(i_sample, e_sample), , drop = FALSE]
  tibble(
    bray_dissimilarity = as.numeric(vegdist(mat, method = "bray")),
    jaccard_dissimilarity = as.numeric(vegdist(mat_pa, method = "jaccard", binary = TRUE))
  )
}

unmatched_metrics <- same_site_depth_unmatched %>%
  mutate(row_id = row_number()) %>%
  group_by(row_id) %>%
  group_modify(~ calc_unmatched(.x$sample_id_i, .x$sample_id_e)) %>%
  ungroup() %>%
  bind_cols(same_site_depth_unmatched) %>%
  select(
    site, site_label, depth_order, depth_label_clean,
    iDNA_sample = sample_id_i,
    eDNA_sample = sample_id_e,
    iDNA_pair = pair_id_i,
    eDNA_pair = pair_id_e,
    bray_dissimilarity,
    jaccard_dissimilarity
  )

write_csv(unmatched_metrics, file.path(report_dir, "unmatched_same_site_depth_similarity.csv"))

paired_vs_unmatched_tests <- tibble(
  metric = c("Bray-Curtis dissimilarity", "Jaccard dissimilarity"),
  paired_mean = c(mean(paired_metrics$bray_dissimilarity), mean(paired_metrics$jaccard_dissimilarity)),
  unmatched_mean = c(mean(unmatched_metrics$bray_dissimilarity), mean(unmatched_metrics$jaccard_dissimilarity)),
  paired_median = c(median(paired_metrics$bray_dissimilarity), median(paired_metrics$jaccard_dissimilarity)),
  unmatched_median = c(median(unmatched_metrics$bray_dissimilarity), median(unmatched_metrics$jaccard_dissimilarity)),
  wilcox_p = c(
    wilcox.test(paired_metrics$bray_dissimilarity, unmatched_metrics$bray_dissimilarity)$p.value,
    wilcox.test(paired_metrics$jaccard_dissimilarity, unmatched_metrics$jaccard_dissimilarity)$p.value
  )
)

write_csv(paired_vs_unmatched_tests, file.path(report_dir, "paired_vs_unmatched_same_site_depth_tests.csv"))

dist_full <- as.matrix(vegdist(sample_matrix, method = "bray"))
nearest_neighbor <- map_dfr(complete_pairs$pair_id, function(pair) {
  ids <- paired_env %>%
    filter(pair_id == pair) %>%
    arrange(dna_type) %>%
    pull(sample_id)
  map_dfr(ids, function(sample_id_now) {
    distances <- dist_full[sample_id_now, ]
    distances <- distances[names(distances) != sample_id_now]
    nearest_sample <- names(which.min(distances))
    tibble(
      pair_id = pair,
      sample_id = sample_id_now,
      dna_type = paired_env$dna_type[match(sample_id_now, paired_env$sample_id)],
      paired_counterpart = setdiff(ids, sample_id_now),
      nearest_sample = nearest_sample,
      nearest_distance = unname(min(distances)),
      counterpart_distance = unname(dist_full[sample_id_now, setdiff(ids, sample_id_now)]),
      nearest_is_counterpart = nearest_sample == setdiff(ids, sample_id_now)
    )
  })
}) %>%
  left_join(complete_pairs, by = "pair_id")

nearest_neighbor_summary <- nearest_neighbor %>%
  group_by(site, site_label) %>%
  summarise(
    n_samples = n(),
    n_pairs = n_distinct(pair_id),
    nearest_counterpart_fraction = mean(nearest_is_counterpart),
    median_counterpart_distance = median(counterpart_distance),
    median_nearest_distance = median(nearest_distance),
    .groups = "drop"
  )

write_csv(nearest_neighbor, file.path(report_dir, "nearest_neighbor_pairing_by_sample.csv"))
write_csv(nearest_neighbor_summary, file.path(report_dir, "nearest_neighbor_pairing_summary.csv"))

paired_env_for_permanova <- paired_env %>%
  mutate(pair_id = factor(pair_id), dna_type = factor(dna_type, levels = DNA_LEVELS))
paired_asv_hell <- decostand(t(paired_asv), method = "hellinger")

set.seed(20260605)
paired_permanova <- adonis2(
  paired_asv_hell ~ dna_type,
  data = paired_env_for_permanova,
  method = "bray",
  permutations = 9999,
  strata = paired_env_for_permanova$pair_id
)
capture.output(paired_permanova, file = file.path(report_dir, "paired_PERMANOVA_pool_effect_strata_pair.txt"))

make_pool_matrix <- function(pool) {
  ids <- paired_env %>%
    filter(dna_type == pool) %>%
    arrange(pair_id) %>%
    pull(sample_id)
  mat <- t(paired_asv[, ids, drop = FALSE])
  rownames(mat) <- paired_env %>%
    filter(dna_type == pool) %>%
    arrange(pair_id) %>%
    pull(pair_id)
  mat[rowSums(mat) > 0, , drop = FALSE]
}

idna_mat <- make_pool_matrix("iDNA")
edna_mat <- make_pool_matrix("eDNA")
common_pair_ids <- intersect(rownames(idna_mat), rownames(edna_mat))
idna_mat <- idna_mat[common_pair_ids, , drop = FALSE]
edna_mat <- edna_mat[common_pair_ids, , drop = FALSE]

idna_hell <- decostand(idna_mat, method = "hellinger")
edna_hell <- decostand(edna_mat, method = "hellinger")

idna_dist <- vegdist(idna_hell, method = "bray")
edna_dist <- vegdist(edna_hell, method = "bray")

idna_pcoa <- wcmdscale(idna_dist, k = 2, eig = TRUE)
edna_pcoa <- wcmdscale(edna_dist, k = 2, eig = TRUE)

set.seed(20260605)
protest_all <- protest(idna_pcoa$points, edna_pcoa$points, permutations = 9999)
set.seed(20260605)
mantel_all <- mantel(idna_dist, edna_dist, method = "spearman", permutations = 9999)

procrustes_summary <- tibble(
  analysis = c("PROTEST Procrustes", "Mantel"),
  statistic = c(protest_all$t0, mantel_all$statistic),
  p_value = c(protest_all$signif, mantel_all$signif),
  n_pairs = length(common_pair_ids)
)

site_procrustes <- map_dfr(SITE_LEVELS, function(s) {
  site_pairs <- complete_pairs %>%
    filter(site == s, pair_id %in% common_pair_ids) %>%
    pull(pair_id)
  if (length(site_pairs) < 5) {
    return(tibble(
      site = s,
      site_label = if_else(s == "NB", "NA", s),
      n_pairs = length(site_pairs),
      protest_r = NA_real_,
      protest_p = NA_real_,
      mantel_r = NA_real_,
      mantel_p = NA_real_
    ))
  }
  id <- idna_hell[site_pairs, , drop = FALSE]
  ed <- edna_hell[site_pairs, , drop = FALSE]
  id_dist <- vegdist(id, method = "bray")
  ed_dist <- vegdist(ed, method = "bray")
  id_pcoa <- wcmdscale(id_dist, k = 2, eig = TRUE)
  ed_pcoa <- wcmdscale(ed_dist, k = 2, eig = TRUE)
  set.seed(20260605)
  pr <- protest(id_pcoa$points, ed_pcoa$points, permutations = 9999)
  set.seed(20260605)
  mt <- mantel(id_dist, ed_dist, method = "spearman", permutations = 9999)
  tibble(
    site = s,
    site_label = if_else(s == "NB", "NA", s),
    n_pairs = length(site_pairs),
    protest_r = pr$t0,
    protest_p = pr$signif,
    mantel_r = mt$statistic,
    mantel_p = mt$signif
  )
})

write_csv(procrustes_summary, file.path(report_dir, "procrustes_mantel_overall.csv"))
write_csv(site_procrustes, file.path(report_dir, "procrustes_mantel_by_site.csv"))

plot_metrics <- paired_metrics %>%
  select(pair_id, site_label, depth_label_clean, bray_dissimilarity,
         jaccard_dissimilarity, shared_asv_fraction, cosine_similarity) %>%
  pivot_longer(
    cols = c(bray_dissimilarity, jaccard_dissimilarity, shared_asv_fraction, cosine_similarity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("bray_dissimilarity", "jaccard_dissimilarity", "shared_asv_fraction", "cosine_similarity"),
      labels = c("Bray-Curtis dissimilarity", "Jaccard dissimilarity",
                 "Shared ASV fraction", "Cosine similarity")
    )
  )

fig_similarity <- ggplot(plot_metrics, aes(depth_label_clean, value, colour = site_label, group = site_label)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), alpha = 0.55, size = 1.8) +
  geom_line(
    data = plot_metrics %>%
      group_by(site_label, depth_label_clean, metric) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop"),
    linewidth = 0.7
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("AZ" = "#d73027", "SG" = "#fc8d59", "LC" = "#1a9850", "NA" = "#4575b4")) +
  labs(x = "Soil depth (cm)", y = "Paired iDNA-eDNA value", colour = "Site") +
  theme_clean(10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey92", colour = "grey35")
  )

comparison_plot_data <- bind_rows(
  paired_metrics %>%
    select(site_label, bray_dissimilarity, jaccard_dissimilarity) %>%
    mutate(comparison = "Matched pair"),
  unmatched_metrics %>%
    select(site_label, bray_dissimilarity, jaccard_dissimilarity) %>%
    mutate(comparison = "Unmatched same site-depth")
) %>%
  pivot_longer(
    cols = c(bray_dissimilarity, jaccard_dissimilarity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    comparison = factor(comparison, levels = c("Matched pair", "Unmatched same site-depth")),
    metric = factor(metric, levels = c("bray_dissimilarity", "jaccard_dissimilarity"),
                    labels = c("Bray-Curtis", "Jaccard"))
  )

fig_comparison <- ggplot(comparison_plot_data, aes(comparison, value, fill = comparison)) +
  geom_boxplot(outlier.shape = NA, width = 0.58, alpha = 0.85) +
  geom_jitter(width = 0.16, alpha = 0.35, size = 1.2) +
  facet_grid(metric ~ site_label) +
  scale_fill_manual(values = c("Matched pair" = "grey45", "Unmatched same site-depth" = "grey80")) +
  labs(x = NULL, y = "iDNA-eDNA dissimilarity") +
  theme_clean(10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", colour = "grey35")
  )

make_tree_plot <- function(output_file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(output_file, width = 11, height = 8.5)
  } else {
    png(output_file, width = 4400, height = 3400, res = 400)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })

  par(mfrow = c(2, 2), mar = c(1.2, 1.2, 3.8, 1), oma = c(0, 0, 2.4, 0))
  for (site_now in SITE_LEVELS) {
    site_pairs <- complete_pairs %>%
      filter(site == site_now) %>%
      pull(pair_id)
    site_env <- paired_env %>%
      filter(pair_id %in% site_pairs) %>%
      mutate(
        depth_for_label = str_replace_all(depths_cm, "_", "-"),
        tree_label = paste0("P", pit, " ", depth_for_label, " ", if_else(dna_type == "iDNA", "i", "e")),
        pool_col = DNA_COLORS[as.character(dna_type)]
      )
    site_mat <- sample_matrix[site_env$sample_id, , drop = FALSE]
    site_mat <- site_mat[, colSums(site_mat) > 0, drop = FALSE]
    site_hell <- decostand(site_mat, method = "hellinger")
    site_tree <- as.phylo(hclust(vegdist(site_hell, method = "bray"), method = "average"))
    original_tips <- site_tree$tip.label
    tip_cols <- site_env$pool_col[match(original_tips, site_env$sample_id)]
    site_tree$tip.label <- site_env$tree_label[match(original_tips, site_env$sample_id)]
    plot(
      site_tree,
      type = "phylogram",
      direction = "rightwards",
      cex = 0.72,
      tip.color = tip_cols,
      label.offset = 0.015,
      no.margin = FALSE
    )
    title(main = if_else(site_now == "NB", "NA", site_now), line = 2.2, cex.main = 1.35, font.main = 2)
    legend(
      "bottomleft",
      legend = c("iDNA", "eDNA"),
      col = DNA_COLORS,
      lty = 1,
      lwd = 2,
      bty = "n",
      cex = 0.8
    )
  }
  mtext("Bray-Curtis clustering of matched iDNA/eDNA archaeal assemblages", outer = TRUE, cex = 1.1, font = 2)
}

make_tree_plot(file.path(report_dir, "paired_iDNA_eDNA_cluster_tree_by_site.pdf"), "pdf")
make_tree_plot(file.path(report_dir, "paired_iDNA_eDNA_cluster_tree_by_site.png"), "png")
make_tree_plot(file.path(supp_fig_dir, "13_FigS13_iDNA_eDNA_cluster_tree_by_site.pdf"), "pdf")
make_tree_plot(file.path(supp_fig_dir, "13_FigS13_iDNA_eDNA_cluster_tree_by_site.png"), "png")

ggsave(file.path(report_dir, "paired_iDNA_eDNA_similarity_depth_profiles.pdf"),
       fig_similarity, width = 8.2, height = 5.8)
ggsave(file.path(report_dir, "paired_iDNA_eDNA_similarity_depth_profiles.png"),
       fig_similarity, width = 8.2, height = 5.8, dpi = 400)
ggsave(file.path(report_dir, "paired_vs_unmatched_same_site_depth.pdf"),
       fig_comparison, width = 8.6, height = 5.6)
ggsave(file.path(report_dir, "paired_vs_unmatched_same_site_depth.png"),
       fig_comparison, width = 8.6, height = 5.6, dpi = 400)

ggsave(file.path(supp_fig_dir, "11_FigS11_iDNA_eDNA_paired_similarity.pdf"),
       fig_similarity, width = 8.2, height = 5.8)
ggsave(file.path(supp_fig_dir, "11_FigS11_iDNA_eDNA_paired_similarity.png"),
       fig_similarity, width = 8.2, height = 5.8, dpi = 400)
ggsave(file.path(supp_fig_dir, "12_FigS12_iDNA_eDNA_matched_vs_unmatched.pdf"),
       fig_comparison, width = 8.6, height = 5.6)
ggsave(file.path(supp_fig_dir, "12_FigS12_iDNA_eDNA_matched_vs_unmatched.png"),
       fig_comparison, width = 8.6, height = 5.6, dpi = 400)

summary_lines <- c(
  "# Paired iDNA-eDNA composition similarity",
  "",
  sprintf("ASV table: %s", asv_path),
  sprintf("Samples retained in 0-60 cm rarefied table: %d", nrow(env)),
  sprintf("Complete iDNA-eDNA pairs at identical site, pit/profile, and depth: %d", nrow(complete_pairs)),
  "",
  "## Paired dissimilarity summary by site",
  capture.output(print(site_summary, n = Inf)),
  "",
  "## Matched pairs vs unmatched samples from the same site-depth",
  capture.output(print(paired_vs_unmatched_tests, n = Inf)),
  "",
  "## Nearest-neighbor counterpart checks",
  capture.output(print(nearest_neighbor_summary, n = Inf)),
  "",
  "## Paired PERMANOVA",
  capture.output(print(paired_permanova)),
  "",
  "## Overall Procrustes / Mantel correspondence between iDNA and eDNA ordination spaces",
  capture.output(print(procrustes_summary, n = Inf)),
  "",
  "## Site-level Procrustes / Mantel",
  capture.output(print(site_procrustes, n = Inf)),
  "",
  "## Interpretation guide",
  "High alpha-diversity similarity means the two pools have comparable diversity levels.",
  "Low Bray-Curtis or Jaccard dissimilarity means the two pools contain similar composition.",
  "A high eDNA read proportion alone is not evidence of rapid turnover; rapid contemporary replenishment would be better supported if matched iDNA-eDNA composition is more similar than unmatched samples from the same site-depth."
)
writeLines(summary_lines, file.path(report_dir, "paired_similarity_summary.md"))

message("Wrote paired pool similarity outputs to: ", report_dir)
