## ============================================================
## Reviewer round-2 transparency tables:
##   1. Good's coverage for retained archaeal samples.
##   2. dbRDA complete-case retention by site, depth, and DNA pool.
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

script_path <- commandArgs(FALSE) |>
  grep(pattern = "^--file=", value = TRUE) |>
  sub(pattern = "^--file=", replacement = "")
if (length(script_path) == 0) {
  script_path <- file.path("script", "36_reviewer2_transparency_tables.R")
}

PROJECT_DIR <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
DATA_DIR <- file.path(PROJECT_DIR, "data")
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis", "reframed_figures")
TABLE_DIR <- file.path(PROJECT_DIR, "01_table")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)

ENV_PATH <- file.path(DATA_DIR, "env_2024.csv")
RAW_ASV_PATH <- file.path(DATA_DIR, "arc_unrarefild", "ASV_Arc_200cm_delete_less200.txt")
RARE_ASV_PATH <- file.path(DATA_DIR, "arc_60cm_202rare", "asv_202.txt")

SITE_LEVELS <- c("AZ", "SG", "LC", "NB")
SITE_LABELS <- c("AZ", "SG", "LC", "NA")
DNA_LEVELS <- c("iDNA", "eDNA")
DEPTH_60_LABELS <- c("0-5", "5-10", "10-20", "20-40", "40-60")
ENV_VARS <- c(
  "pH", "Conductivity", "moisture", "CN",
  "Feo", "Alo", "Mno", "Sio",
  "NH4", "NO3", "Po", "Pi"
)

display_site <- function(x) recode(as.character(x), NB = "NA")

env_raw <- read.csv(ENV_PATH, header = TRUE, check.names = FALSE)
env_60 <- env_raw %>%
  filter(depth < 60) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    site_label = factor(display_site(site), levels = SITE_LABELS),
    dna_type = factor(dna_type, levels = DNA_LEVELS),
    depth_label_display = factor(gsub("_", "-", depths_cm), levels = DEPTH_60_LABELS)
  )

raw_asv <- read.delim(
  RAW_ASV_PATH,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
rare_asv <- read.delim(
  RARE_ASV_PATH,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

retained_samples <- intersect(colnames(rare_asv), env_60$sample_id)
raw_asv_retained <- raw_asv[, retained_samples, drop = FALSE]
rare_asv_retained <- rare_asv[, retained_samples, drop = FALSE]

calc_goods <- function(mat, basis) {
  tibble(
    sample_id = colnames(mat),
    read_depth = colSums(mat),
    observed_asvs = colSums(mat > 0),
    singleton_asvs = colSums(mat == 1),
    goods_coverage = if_else(read_depth > 0, 1 - singleton_asvs / read_depth, NA_real_),
    coverage_basis = basis
  )
}

coverage_by_sample <- bind_rows(
  calc_goods(raw_asv_retained, "pre-rarefaction retained archaeal table"),
  calc_goods(rare_asv_retained, "202-read rarefied archaeal table")
) %>%
  left_join(
    env_60 %>%
      select(sample_id, site_label, dna_type, depth_label_display),
    by = "sample_id"
  ) %>%
  arrange(coverage_basis, site_label, dna_type, depth_label_display, sample_id)

coverage_summary <- coverage_by_sample %>%
  group_by(coverage_basis, site_label, dna_type) %>%
  summarise(
    n_samples = n(),
    median_read_depth = median(read_depth, na.rm = TRUE),
    median_observed_asvs = median(observed_asvs, na.rm = TRUE),
    median_singleton_asvs = median(singleton_asvs, na.rm = TRUE),
    median_goods_coverage = median(goods_coverage, na.rm = TRUE),
    min_goods_coverage = min(goods_coverage, na.rm = TRUE),
    max_goods_coverage = max(goods_coverage, na.rm = TRUE),
    .groups = "drop"
  )

dbRDA_env_60 <- env_60 %>%
  filter(sample_id %in% colnames(rare_asv))

complete_case_retention <- dbRDA_env_60 %>%
  mutate(
    retained_for_dbRDA = if_all(all_of(ENV_VARS), ~ !is.na(.x))
  ) %>%
  group_by(site_label, dna_type, depth_label_display) %>%
  summarise(
    n_total = n(),
    n_dbRDA_complete = sum(retained_for_dbRDA),
    n_dropped = n_total - n_dbRDA_complete,
    .groups = "drop"
  ) %>%
  arrange(site_label, dna_type, depth_label_display)

complete_case_summary <- complete_case_retention %>%
  group_by(site_label, dna_type) %>%
  summarise(
    n_total = sum(n_total),
    n_dbRDA_complete = sum(n_dbRDA_complete),
    n_dropped = sum(n_dropped),
    retention_fraction = n_dbRDA_complete / n_total,
    .groups = "drop"
  )

write.csv(
  coverage_by_sample,
  file.path(ANALYSIS_DIR, "reviewer2_goods_coverage_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  coverage_summary,
  file.path(ANALYSIS_DIR, "reviewer2_goods_coverage_summary.csv"),
  row.names = FALSE
)
write.csv(
  complete_case_retention,
  file.path(ANALYSIS_DIR, "reviewer2_dbRDA_complete_case_retention_site_depth.csv"),
  row.names = FALSE
)
write.csv(
  complete_case_summary,
  file.path(ANALYSIS_DIR, "reviewer2_dbRDA_complete_case_retention_summary.csv"),
  row.names = FALSE
)

message("Wrote Good's coverage tables to: ", ANALYSIS_DIR)
message("Wrote dbRDA retention tables to: ", ANALYSIS_DIR)
