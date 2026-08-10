#!/usr/bin/env Rscript

# Reproduce the reporting quantities requested during the Round39 peer review.
# This script is descriptive only: it does not remove ASVs or alter any analysis.

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

raw_path <- file.path(project_root, "data", "file_A.txt")
fig1_path <- file.path(
  project_root,
  "analysis",
  "reframed_figures",
  "Fig1_archaeal_read_fraction_source.csv"
)
rarefied_path <- file.path(
  project_root,
  "data",
  "arc_60cm_202rare",
  "asv_202.txt"
)
out_dir <- file.path(project_root, "analysis", "reproducibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.delim(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
fig1 <- read.csv(
  fig1_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = ""
)
rarefied <- read.delim(
  rarefied_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxonomy_cols <- c("asv_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
library_cols <- setdiff(names(raw), taxonomy_cols)
control_cols <- library_cols[grepl("(^NTC_|_NC_|^NC_)", library_cols)]
biological_cols <- setdiff(library_cols, control_cols)

stopifnot(
  length(control_cols) == 27L,
  nrow(fig1) == 95L,
  length(unique(fig1$sample_id)) == 95L,
  all(fig1$sample_id %in% biological_cols),
  nrow(rarefied) == 658L,
  ncol(rarefied) == 95L,
  all(colSums(rarefied) == 202L)
)

archaea <- raw[raw$Kingdom == "Archaea", , drop = FALSE]

control_summary <- data.frame(
  control_id = control_cols,
  control_type = ifelse(grepl("^NTC_", control_cols), "no-template control", "negative control"),
  total_prokaryotic_reads = as.integer(colSums(raw[, control_cols, drop = FALSE])),
  archaeal_reads = as.integer(colSums(archaea[, control_cols, drop = FALSE])),
  stringsAsFactors = FALSE
)
control_summary$archaeal_percent <- ifelse(
  control_summary$total_prokaryotic_reads > 0,
  100 * control_summary$archaeal_reads / control_summary$total_prokaryotic_reads,
  NA_real_
)

write.csv(
  control_summary,
  file.path(out_dir, "round40_negative_control_library_audit.csv"),
  row.names = FALSE,
  na = ""
)

control_counts <- archaea[, control_cols, drop = FALSE]
biological_counts <- archaea[, biological_cols, drop = FALSE]
retained_counts <- archaea[, fig1$sample_id, drop = FALSE]
na_subsoil_ids <- fig1$sample_id[
  fig1$site == "NA" &
    fig1$dna_type == "eDNA" &
    fig1$depths_cm_re %in% c("20-40", "40-60")
]
na_subsoil_counts <- archaea[, na_subsoil_ids, drop = FALSE]

asv_audit <- data.frame(
  asv_id = archaea$asv_id,
  phylum = archaea$Phylum,
  class = archaea$Class,
  order = archaea$Order,
  family = archaea$Family,
  genus = archaea$Genus,
  control_reads_total = as.integer(rowSums(control_counts)),
  controls_detected_n = as.integer(rowSums(control_counts > 0)),
  control_prevalence_percent = 100 * rowSums(control_counts > 0) / length(control_cols),
  maximum_reads_in_one_control = as.integer(apply(control_counts, 1, max)),
  all_biological_reads_total = as.integer(rowSums(biological_counts)),
  all_biological_libraries_detected_n = as.integer(rowSums(biological_counts > 0)),
  all_biological_prevalence_percent = 100 * rowSums(biological_counts > 0) / length(biological_cols),
  retained_95_reads_total = as.integer(rowSums(retained_counts)),
  retained_95_samples_detected_n = as.integer(rowSums(retained_counts > 0)),
  retained_95_prevalence_percent = 100 * rowSums(retained_counts > 0) / ncol(retained_counts),
  NA_eDNA_20_60cm_reads_total = as.integer(rowSums(na_subsoil_counts)),
  NA_eDNA_20_60cm_samples_detected_n = as.integer(rowSums(na_subsoil_counts > 0)),
  stringsAsFactors = FALSE
)
asv_audit <- asv_audit[asv_audit$control_reads_total > 0, , drop = FALSE]
asv_audit <- asv_audit[
  order(
    asv_audit$control_reads_total,
    asv_audit$controls_detected_n,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  asv_audit,
  file.path(out_dir, "round40_negative_control_asv_audit.csv"),
  row.names = FALSE,
  na = ""
)

subsoil <- fig1[
  fig1$site %in% c("LC", "NA") &
    fig1$depths_cm_re %in% c("20-40", "40-60"),
  ,
  drop = FALSE
]
subsoil_groups <- split(
  subsoil,
  interaction(subsoil$site, subsoil$dna_type, subsoil$depths_cm_re, drop = TRUE)
)
subsoil_summary <- do.call(
  rbind,
  lapply(subsoil_groups, function(group) {
    data.frame(
      site = group$site[1],
      dna_pool = group$dna_type[1],
      depth_cm = group$depths_cm_re[1],
      n_profiles = nrow(group),
      mean_archaeal_percent = mean(group$arcpercent),
      sd_archaeal_percent = stats::sd(group$arcpercent),
      minimum_archaeal_percent = min(group$arcpercent),
      maximum_archaeal_percent = max(group$arcpercent),
      minimum_archaeal_reads = min(group$archaeal_reads),
      maximum_archaeal_reads = max(group$archaeal_reads),
      stringsAsFactors = FALSE
    )
  })
)
subsoil_summary <- subsoil_summary[
  order(subsoil_summary$site, subsoil_summary$depth_cm, subsoil_summary$dna_pool),
  ,
  drop = FALSE
]
row.names(subsoil_summary) <- NULL

write.csv(
  subsoil_summary,
  file.path(out_dir, "round40_fig1_subsoil_profile_summary.csv"),
  row.names = FALSE
)

pooled_counts <- rowSums(rarefied)
per_sample_rare_fraction <- vapply(
  rarefied,
  function(values) {
    observed <- sum(values > 0)
    if (observed == 0) return(NA_real_)
    (sum(values == 1) + sum(values == 2)) / observed
  },
  numeric(1)
)

rarefaction_summary <- data.frame(
  metric = c(
    "total_ASVs",
    "pooled_singleton_ASVs",
    "pooled_singleton_percent",
    "pooled_doubleton_ASVs",
    "pooled_doubleton_percent",
    "pooled_singleton_plus_doubleton_ASVs",
    "pooled_singleton_plus_doubleton_percent",
    "median_per_sample_singleton_plus_doubleton_fraction",
    "minimum_per_sample_singleton_plus_doubleton_fraction",
    "maximum_per_sample_singleton_plus_doubleton_fraction"
  ),
  value = c(
    nrow(rarefied),
    sum(pooled_counts == 1),
    100 * mean(pooled_counts == 1),
    sum(pooled_counts == 2),
    100 * mean(pooled_counts == 2),
    sum(pooled_counts <= 2),
    100 * mean(pooled_counts <= 2),
    median(per_sample_rare_fraction, na.rm = TRUE),
    min(per_sample_rare_fraction, na.rm = TRUE),
    max(per_sample_rare_fraction, na.rm = TRUE)
  ),
  definition = c(
    "ASVs retained in the pooled 95-sample 202-read table",
    "ASVs with one read after pooling all 95 rarefied samples",
    "Pooled singleton ASVs divided by 658 ASVs",
    "ASVs with two reads after pooling all 95 rarefied samples",
    "Pooled doubleton ASVs divided by 658 ASVs",
    "ASVs with one or two reads after pooling all 95 rarefied samples",
    "Pooled singleton plus doubleton ASVs divided by 658 ASVs",
    "Median within-sample fraction of observed ASVs represented by one or two reads",
    "Minimum within-sample fraction of observed ASVs represented by one or two reads",
    "Maximum within-sample fraction of observed ASVs represented by one or two reads"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  rarefaction_summary,
  file.path(out_dir, "round40_202read_rare_asv_audit.csv"),
  row.names = FALSE
)

message(
  "Negative controls: n=", nrow(control_summary),
  "; archaeal reads median=", median(control_summary$archaeal_reads),
  ", max=", max(control_summary$archaeal_reads),
  ", total=", sum(control_summary$archaeal_reads)
)
message(
  "Control-detected archaeal ASVs: ", nrow(asv_audit),
  "; no contaminant classification was applied."
)
message(
  "Pooled 202-read singleton+doubleton ASVs: ",
  sum(pooled_counts <= 2), " / ", nrow(rarefied),
  " (", sprintf("%.1f", 100 * mean(pooled_counts <= 2)), "%)."
)
