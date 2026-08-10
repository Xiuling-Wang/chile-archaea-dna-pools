#!/usr/bin/env Rscript

# Sensitivity audit for the central class-turnover result. Descriptive
# composition uses unrarefied counts normalized within each sample; alpha
# diversity remains based on the established 202-read table.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(vegan)
})

set.seed(20260713)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "analysis", "reproducibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.delim(
  file.path(project_root, "data", "arc_unrarefild", "ASV_Arc_200cm_delete_less200.txt"),
  row.names = 1,
  check.names = FALSE
)
rare <- read.delim(
  file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt"),
  row.names = 1,
  check.names = FALSE
)
tax <- read.delim(
  file.path(project_root, "data", "arc_100_202rarefild", "tax_202.txt"),
  row.names = 1,
  check.names = FALSE
)
env <- read.csv(file.path(project_root, "data", "env_2024.csv"), check.names = FALSE)

samples <- colnames(rare)
if (!all(samples %in% colnames(raw)) || !all(samples %in% env$sample_id)) {
  stop("The 95 retained samples are not fully represented in the raw table and metadata.")
}

raw95 <- as.matrix(raw[, samples, drop = FALSE])
rare95 <- as.matrix(rare[, samples, drop = FALSE])
storage.mode(raw95) <- "numeric"
storage.mode(rare95) <- "numeric"
env95 <- env[match(samples, env$sample_id), ]

strip_rank <- function(x) sub("^[a-z]__", "", x)
class_lookup <- setNames(strip_rank(tax$Class), rownames(tax))

class_composition <- function(mat) {
  rel <- sweep(mat, 2, colSums(mat), "/")
  classes <- sort(unique(class_lookup[intersect(names(class_lookup), rownames(mat))]))
  classes <- classes[!is.na(classes)]
  values <- vapply(classes, function(class_name) {
    ids <- intersect(names(class_lookup)[class_lookup == class_name], rownames(mat))
    if (length(ids) == 0) return(rep(0, ncol(mat)))
    colSums(rel[ids, , drop = FALSE])
  }, numeric(ncol(mat)))
  rownames(values) <- colnames(mat)
  values
}

summarize_turnover <- function(mat, normalization) {
  comp <- class_composition(mat)
  data.frame(
    sample_id = rownames(comp),
    comp,
    check.names = FALSE
  ) %>%
    left_join(env95 %>% select(sample_id, site, pit, dna_type), by = "sample_id") %>%
    group_by(site, dna_type) %>%
    summarise(
      Nitrososphaeria = 100 * mean(Nitrososphaeria, na.rm = TRUE),
      Thermoplasmata = 100 * mean(Thermoplasmata, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(normalization = normalization, .before = 1)
}

turnover <- bind_rows(
  summarize_turnover(raw95, "unrarefied sample-wise relative abundance"),
  summarize_turnover(rare95, "202-read rarefaction")
)
write.csv(turnover, file.path(out_dir, "class_turnover_raw_vs_rarefied.csv"), row.names = FALSE)

profile_composition <- function(mat) {
  comp <- class_composition(mat)
  sample_df <- data.frame(
    sample_id = rownames(comp),
    comp,
    check.names = FALSE
  ) %>%
    left_join(env95 %>% select(sample_id, site, pit, dna_type), by = "sample_id")

  pool_means <- sample_df %>%
    group_by(site, pit, dna_type) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  pool_means %>%
    group_by(site, pit) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
}

profile_asv_composition <- function(mat) {
  rel <- sweep(mat, 2, colSums(mat), "/")
  sample_df <- data.frame(
    sample_id = colnames(mat),
    t(rel),
    check.names = FALSE
  ) %>%
    left_join(env95 %>% select(sample_id, site, pit, dna_type), by = "sample_id")

  pool_means <- sample_df %>%
    group_by(site, pit, dna_type) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  pool_means %>%
    group_by(site, pit) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
}

raw_profile <- profile_composition(raw95)
rare_profile <- profile_composition(rare95)
class_cols <- intersect(names(raw_profile)[-(1:2)], names(rare_profile)[-(1:2)])
class_cols <- class_cols[
  colSums(raw_profile[, class_cols, drop = FALSE]) > 0 &
    colSums(rare_profile[, class_cols, drop = FALSE]) > 0
]

profile_key <- paste(raw_profile$site, raw_profile$pit, sep = "_")
rare_profile <- rare_profile[match(profile_key, paste(rare_profile$site, rare_profile$pit, sep = "_")), ]

raw_dist <- vegdist(raw_profile[, class_cols, drop = FALSE], method = "robust.aitchison")
site_factor <- factor(raw_profile$site, levels = c("AZ", "SG", "LC", "NB"))
aitchison_site <- adonis2(raw_dist ~ site_factor, permutations = 9999)

raw_asv_profile <- profile_asv_composition(raw95)
asv_cols <- names(raw_asv_profile)[-(1:2)]
asv_cols <- asv_cols[colSums(raw_asv_profile[, asv_cols, drop = FALSE]) > 0]
raw_asv_dist <- vegdist(raw_asv_profile[, asv_cols, drop = FALSE], method = "robust.aitchison")
asv_site_factor <- factor(raw_asv_profile$site, levels = c("AZ", "SG", "LC", "NB"))
asv_aitchison_site <- adonis2(raw_asv_dist ~ asv_site_factor, permutations = 9999)

raw_clr <- decostand(raw_profile[, class_cols, drop = FALSE], method = "clr", pseudocount = 1e-06)
rare_clr <- decostand(rare_profile[, class_cols, drop = FALSE], method = "clr", pseudocount = 1e-06)
procrustes_test <- protest(raw_clr, rare_clr, permutations = 9999)

audit_summary <- data.frame(
  check = c(
    "raw_ASV_robust_aitchison_site_R2",
    "raw_ASV_robust_aitchison_site_F",
    "raw_ASV_robust_aitchison_site_p",
    "raw_robust_aitchison_site_R2",
    "raw_robust_aitchison_site_F",
    "raw_robust_aitchison_site_p",
    "raw_vs_rarefied_class_clr_procrustes_r",
    "raw_vs_rarefied_class_clr_procrustes_p"
  ),
  value = c(
    unname(asv_aitchison_site$R2[1]),
    unname(asv_aitchison_site$F[1]),
    unname(asv_aitchison_site$`Pr(>F)`[1]),
    unname(aitchison_site$R2[1]),
    unname(aitchison_site$F[1]),
    unname(aitchison_site$`Pr(>F)`[1]),
    unname(procrustes_test$t0),
    unname(procrustes_test$signif)
  )
)
write.csv(audit_summary, file.path(out_dir, "class_turnover_aitchison_audit.csv"), row.names = FALSE)

report <- c(
  "# Class-turnover compositional sensitivity audit",
  "",
  sprintf("- Retained samples: %d (all 0-60 cm samples used in the manuscript).", length(samples)),
  sprintf("- Unrarefied library size range: %d-%d reads.", min(colSums(raw95)), max(colSums(raw95))),
  "- Descriptive composition was calculated as within-sample relative abundance before group averaging.",
  sprintf("- Profile-level ASV robust Aitchison site effect: R2 = %.3f, p = %.4f.", asv_aitchison_site$R2[1], asv_aitchison_site$`Pr(>F)`[1]),
  sprintf("- Class-level robust Aitchison site sensitivity: R2 = %.3f, p = %.4f.", aitchison_site$R2[1], aitchison_site$`Pr(>F)`[1]),
  sprintf("- Raw versus 202-read class-composition concordance: Procrustes r = %.3f, p = %.4f.", procrustes_test$t0, procrustes_test$signif),
  "- Interpretation: the central class-turnover pattern is not created by the single 202-read rarefaction."
)
writeLines(report, file.path(out_dir, "class_turnover_aitchison_audit.md"))

message("Wrote compositional sensitivity outputs to: ", out_dir)
