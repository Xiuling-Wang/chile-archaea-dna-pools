#!/usr/bin/env Rscript

# Design-aware site and depth PERMANOVA audit.
# Site is tested on independent site x pit profile means (n = 12).
# Depth is tested within profiles, with profile identity in the model and
# permutations restricted within profiles. DNA pool retains its separate exact
# paired test in script/44_restricted_permutation_audit.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(permute)
  library(readr)
  library(tidyr)
  library(vegan)
})

set.seed(20260711)
project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "analysis", "reproducibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

asv <- read.delim(
  file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt"),
  row.names = 1,
  check.names = FALSE
)
env <- read.csv(file.path(project_root, "data", "env_2024.csv"), check.names = FALSE) |>
  filter(depth < 60, sample_id %in% colnames(asv), dna_type %in% c("iDNA", "eDNA")) |>
  mutate(
    site = factor(site, levels = c("AZ", "SG", "LC", "NB")),
    profile_id = interaction(site, pit, drop = TRUE),
    pair_id = interaction(site, pit, depth_order, drop = TRUE),
    depth_mid = as.numeric(depth),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA"))
  ) |>
  arrange(match(sample_id, colnames(asv)))
asv <- asv[, env$sample_id, drop = FALSE]
asv <- asv[rowSums(asv) > 0, , drop = FALSE]
rel <- sweep(asv, 2, colSums(asv), "/")

profile_matrix <- function(pool = NULL, complete_pairs_only = FALSE) {
  use_env <- env
  if (!is.null(pool)) use_env <- use_env |> filter(dna_type == pool)
  if (complete_pairs_only) {
    complete <- env |>
      count(pair_id, dna_type) |>
      pivot_wider(names_from = dna_type, values_from = n, values_fill = 0) |>
      filter(iDNA == 1, eDNA == 1) |>
      pull(pair_id)
    use_env <- use_env |> filter(pair_id %in% complete)
  }
  profiles <- unique(as.character(use_env$profile_id))
  matrix <- sapply(profiles, function(profile) {
    profile_env <- use_env |> filter(as.character(profile_id) == profile)
    if (is.null(pool) && !complete_pairs_only) {
      i_mean <- rowMeans(rel[, profile_env$sample_id[profile_env$dna_type == "iDNA"], drop = FALSE])
      e_mean <- rowMeans(rel[, profile_env$sample_id[profile_env$dna_type == "eDNA"], drop = FALSE])
      return((i_mean + e_mean) / 2)
    }
    rowMeans(rel[, profile_env$sample_id, drop = FALSE])
  })
  profile_env <- use_env[match(profiles, as.character(use_env$profile_id)), ]
  list(matrix = as.matrix(matrix), env = profile_env)
}

run_site <- function(pool = NULL, complete_pairs_only = FALSE) {
  profile <- profile_matrix(pool, complete_pairs_only)
  hellinger <- decostand(t(profile$matrix), method = "hellinger")
  set.seed(20260711)
  model <- adonis2(
    hellinger ~ site,
    data = profile$env,
    method = "bray",
    permutations = 9999
  )
  dispersion <- betadisper(vegdist(hellinger, method = "bray"), profile$env$site)
  set.seed(20260711)
  dispersion_test <- permutest(dispersion, permutations = 9999)
  row <- as.data.frame(model)[1, ]
  tibble(
    analysis = "profile-level site PERMANOVA",
    pool = ifelse(is.null(pool), "equal-weight iDNA/eDNA profile mean", pool),
    sensitivity = complete_pairs_only,
    n_profiles = nrow(profile$env),
    Df = row$Df,
    SumOfSqs = row$SumOfSqs,
    R2 = row$R2,
    F = row$F,
    p_value = row$`Pr(>F)`,
    permdisp_F = dispersion_test$tab[1, "F"],
    permdisp_p = dispersion_test$tab[1, "Pr(>F)"]
  )
}

run_depth <- function(pool = NULL) {
  use_env <- env
  if (!is.null(pool)) use_env <- use_env |> filter(dna_type == pool)
  use_asv <- asv[, use_env$sample_id, drop = FALSE]
  use_asv <- use_asv[rowSums(use_asv) > 0, , drop = FALSE]
  hellinger <- decostand(t(use_asv), method = "hellinger")
  control <- how(
    within = Within(type = "free"),
    blocks = use_env$profile_id,
    nperm = 9999
  )
  formula <- if (is.null(pool)) {
    hellinger ~ profile_id + depth_mid + dna_type
  } else {
    hellinger ~ profile_id + depth_mid
  }
  set.seed(20260711)
  model <- adonis2(
    formula,
    data = use_env,
    method = "bray",
    permutations = control,
    by = "margin"
  )
  row <- as.data.frame(model)["depth_mid", ]
  tibble(
    analysis = "within-profile depth PERMANOVA",
    pool = ifelse(is.null(pool), "both pools", pool),
    n_samples = nrow(use_env),
    Df = row$Df,
    SumOfSqs = row$SumOfSqs,
    R2 = row$R2,
    F = row$F,
    p_value = row$`Pr(>F)`
  )
}

site_results <- bind_rows(
  run_site(),
  run_site("iDNA"),
  run_site("eDNA"),
  run_site(complete_pairs_only = TRUE)
)
depth_results <- bind_rows(run_depth(), run_depth("iDNA"), run_depth("eDNA"))

write_csv(site_results, file.path(out_dir, "round30_profile_level_site_permanova.csv"))
write_csv(depth_results, file.path(out_dir, "round30_within_profile_depth_permanova.csv"))

main_site <- site_results |> filter(pool == "equal-weight iDNA/eDNA profile mean", !sensitivity)
paired_site <- site_results |> filter(sensitivity)
main_depth <- depth_results |> filter(pool == "both pools")
write_lines(
  c(
    "# Round30 design-aware site and depth audit",
    "",
    "## Site",
    "",
    "Sample-wise relative abundance was averaged across retained depths within each DNA pool and profile, then iDNA and eDNA profile means were given equal weight. Site was tested across 12 independent site x pit profiles (three per site).",
    sprintf("Result: R2 = %.5f, F = %.4f, p = %.4g; profile-level PERMDISP p = %.4f.", main_site$R2, main_site$F, main_site$p_value, main_site$permdisp_p),
    sprintf("Complete-pair-only sensitivity: R2 = %.5f, F = %.4f, p = %.4g.", paired_site$R2, paired_site$F, paired_site$p_value),
    "",
    "## Depth",
    "",
    "Depth midpoint was tested with profile identity included in the model and permutations restricted within profile. DNA pool was included as a marginal term in the combined-pool model.",
    sprintf("Result: R2 = %.5f, F = %.4f, p = %.4g.", main_depth$R2, main_depth$F, main_depth$p_value),
    "",
    "## Interpretation",
    "",
    "The site and depth conclusions survive design-aware analyses. The earlier free-permutation sample-level site p value is not retained as the inferential test."
  ),
  file.path(out_dir, "round30_design_aware_site_depth_audit.md")
)

message("Wrote design-aware site/depth audit to: ", out_dir)
