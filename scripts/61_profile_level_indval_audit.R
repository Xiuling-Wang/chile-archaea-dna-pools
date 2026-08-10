#!/usr/bin/env Rscript

# Design-corrected indicator-ASV audit. Site is a profile-level treatment, so
# the independent replicate is site x pit (not a depth sample within a pit).

suppressPackageStartupMessages({
  library(dplyr)
  library(labdsv)
  library(readr)
})

set.seed(20260711)
project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
asv <- read.delim(
  file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt"),
  row.names = 1,
  check.names = FALSE
)
env <- read.csv(file.path(project_root, "data", "env_2024.csv"), check.names = FALSE) |>
  filter(sample_id %in% colnames(asv), dna_type %in% c("iDNA", "eDNA")) |>
  mutate(
    site = recode(site, NB = "NA"),
    profile_id = interaction(site, pit, drop = TRUE)
  )

site_levels <- c("AZ", "SG", "LC", "NA")
out_dir <- file.path(project_root, "analysis", "reframed_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_pool <- function(pool, min_retained_depths = 1L) {
  pool_env <- env |>
    filter(dna_type == pool) |>
    arrange(site, pit, depth_order)
  retained_profiles <- pool_env |>
    count(profile_id, name = "retained_depth_samples") |>
    filter(retained_depth_samples >= min_retained_depths) |>
    pull(profile_id) |>
    as.character()
  pool_env <- pool_env |>
    filter(as.character(profile_id) %in% retained_profiles)
  pool_asv <- asv[, pool_env$sample_id, drop = FALSE]
  pool_asv <- pool_asv[rowSums(pool_asv) > 0, , drop = FALSE]
  rel <- sweep(pool_asv, 2, colSums(pool_asv), "/") * 100
  rel <- rel[rowMeans(rel) > 0.01, , drop = FALSE]

  profiles <- unique(as.character(pool_env$profile_id))
  profile_means <- sapply(
    profiles,
    function(profile) rowMeans(rel[, as.character(pool_env$profile_id) == profile, drop = FALSE])
  )
  profile_means <- as.matrix(profile_means)
  colnames(profile_means) <- profiles
  groups <- factor(sub("\\..*$", "", profiles), levels = site_levels)

  result <- indval(t(profile_means), as.integer(groups), numitr = 9999)
  adjusted <- p.adjust(result$pval, method = "BH")
  selected <- result$indcls > 0.8 & adjusted < 0.05

  indicators <- tibble(
    pool = pool,
    ASV = rownames(profile_means)[selected],
    site = site_levels[result$maxcls[selected]],
    indval = result$indcls[selected],
    p_value = result$pval[selected],
    p_bh = adjusted[selected]
  )
  retention <- pool_env |>
    count(site, pit, profile_id, name = "retained_depth_samples") |>
    mutate(pool = pool, .before = 1)
  list(indicators = indicators, retention = retention, n_tested = nrow(profile_means))
}

results <- lapply(c("iDNA", "eDNA"), run_pool)
sensitivity <- lapply(c("iDNA", "eDNA"), run_pool, min_retained_depths = 3L)
indicators <- bind_rows(lapply(results, `[[`, "indicators"))
retention <- bind_rows(lapply(results, `[[`, "retention"))

write_csv(indicators, file.path(out_dir, "TableS5_indicators_profilelevel_BHFDR.csv"))
write_csv(retention, file.path(out_dir, "TableS5_profile_retention.csv"))

summary_lines <- c(
  "# Profile-level indicator-ASV audit",
  "",
  "Independent replicates are site x pit profiles. Sample-wise relative abundances were averaged across the retained depths within each profile before IndVal.",
  "ASVs with mean relative abundance <= 0.01% were excluded; 9,999 permutations; IndVal > 0.8; BH-adjusted p < 0.05.",
  "",
  sprintf("iDNA: %d tested ASVs; %d retained indicators.", results[[1]]$n_tested, nrow(results[[1]]$indicators)),
  sprintf("eDNA: %d tested ASVs; %d retained indicators.", results[[2]]$n_tested, nrow(results[[2]]$indicators)),
  sprintf("Sensitivity requiring at least three retained depths per profile: iDNA %d indicators; eDNA %d indicators.", nrow(sensitivity[[1]]$indicators), nrow(sensitivity[[2]]$indicators)),
  "",
  "The earlier sample-level result (0 iDNA and 8 eDNA indicators) is not retained because multiple depths from the same profile are not independent site-level replicates."
)
write_lines(summary_lines, file.path(out_dir, "TableS5_profilelevel_audit.md"))

message("Wrote profile-level IndVal audit to: ", out_dir)
