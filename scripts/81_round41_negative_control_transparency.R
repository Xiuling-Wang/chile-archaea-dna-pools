#!/usr/bin/env Rscript

# Build the descriptive negative-control audit used in Round41.
# This script classifies no ASV as a contaminant and changes no analysis table.

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

raw_path <- file.path(project_root, "data", "file_A.txt")
pre_rarefaction_path <- file.path(
  project_root,
  "data",
  "arc_unrarefild",
  "ASV_Arc_200cm_delete_less200.txt"
)
fig1_path <- file.path(
  project_root,
  "analysis",
  "reframed_figures",
  "Fig1_archaeal_read_fraction_source.csv"
)
out_dir <- file.path(project_root, "analysis", "reproducibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.delim(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
pre_rarefaction <- read.delim(
  pre_rarefaction_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fig1 <- read.csv(
  fig1_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = ""
)

taxonomy_cols <- c("asv_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
library_cols <- setdiff(names(raw), taxonomy_cols)

# Negative extraction/process controls contain NC in the identifier; NTC denotes
# no-template PCR controls. Positive controls are excluded from both the negative-
# control and biological-library summaries.
negative_control_cols <- library_cols[
  grepl("(^NTC_|_NC_|^NC_)", library_cols)
]
positive_control_cols <- library_cols[grepl("^PC_", library_cols)]
biological_cols <- setdiff(
  library_cols,
  c(negative_control_cols, positive_control_cols)
)

stopifnot(
  length(negative_control_cols) == 27L,
  length(positive_control_cols) == 6L,
  length(biological_cols) == 320L,
  nrow(fig1) == 95L,
  length(unique(fig1$sample_id)) == 95L,
  all(fig1$sample_id %in% biological_cols)
)

archaea <- raw[raw$Kingdom == "Archaea", , drop = FALSE]

shared_pre_rarefaction_cols <- intersect(
  colnames(pre_rarefaction),
  biological_cols
)
raw_shared_counts <- as.matrix(
  archaea[
    match(rownames(pre_rarefaction), archaea$asv_id),
    shared_pre_rarefaction_cols,
    drop = FALSE
  ]
)
storage.mode(raw_shared_counts) <- "numeric"
pre_rarefaction_shared_counts <- as.matrix(
  pre_rarefaction[, shared_pre_rarefaction_cols, drop = FALSE]
)
storage.mode(pre_rarefaction_shared_counts) <- "numeric"

stopifnot(
  nrow(archaea) == 1166L,
  nrow(pre_rarefaction) == 1166L,
  length(shared_pre_rarefaction_cols) == 214L,
  all(colnames(pre_rarefaction) %in% biological_cols),
  identical(archaea$asv_id, rownames(pre_rarefaction)),
  identical(unname(raw_shared_counts), unname(pre_rarefaction_shared_counts))
)

control_summary <- data.frame(
  control_id = negative_control_cols,
  control_type = ifelse(
    grepl("^NTC_", negative_control_cols),
    "no-template PCR control",
    "negative extraction/process control"
  ),
  total_prokaryotic_reads = as.integer(
    colSums(raw[, negative_control_cols, drop = FALSE])
  ),
  archaeal_reads = as.integer(
    colSums(archaea[, negative_control_cols, drop = FALSE])
  ),
  stringsAsFactors = FALSE
)
control_summary$archaeal_percent <- ifelse(
  control_summary$total_prokaryotic_reads > 0,
  100 * control_summary$archaeal_reads / control_summary$total_prokaryotic_reads,
  NA_real_
)

write.csv(
  control_summary,
  file.path(out_dir, "round41_negative_control_library_audit.csv"),
  row.names = FALSE,
  na = ""
)

control_counts <- archaea[, negative_control_cols, drop = FALSE]
biological_counts <- archaea[, biological_cols, drop = FALSE]
retained_counts <- archaea[, fig1$sample_id, drop = FALSE]
na_subsoil_ids <- fig1$sample_id[
  fig1$site == "NA" &
    fig1$dna_type == "eDNA" &
    fig1$depths_cm_re %in% c("20-40", "40-60")
]
stopifnot(length(na_subsoil_ids) == 6L)
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
  control_prevalence_percent = 100 * rowSums(control_counts > 0) /
    length(negative_control_cols),
  maximum_reads_in_one_control = as.integer(apply(control_counts, 1, max)),
  retained_95_reads_total = as.integer(rowSums(retained_counts)),
  retained_95_samples_detected_n = as.integer(rowSums(retained_counts > 0)),
  retained_95_prevalence_percent = 100 * rowSums(retained_counts > 0) /
    ncol(retained_counts),
  all_biological_reads_total = as.integer(rowSums(biological_counts)),
  all_biological_libraries_detected_n = as.integer(rowSums(biological_counts > 0)),
  all_biological_prevalence_percent = 100 * rowSums(biological_counts > 0) /
    length(biological_cols),
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
row.names(asv_audit) <- NULL

write.csv(
  asv_audit,
  file.path(out_dir, "round41_negative_control_asv_audit.csv"),
  row.names = FALSE,
  na = ""
)

na_subsoil_summary <- fig1[
  match(na_subsoil_ids, fig1$sample_id),
  c(
    "sample_id",
    "site",
    "dna_type",
    "depth",
    "depths_cm_re",
    "total_reads",
    "archaeal_reads",
    "arcpercent"
  ),
  drop = FALSE
]
names(na_subsoil_summary) <- c(
  "sample_id",
  "site",
  "dna_pool",
  "depth_midpoint_cm",
  "depth_interval_cm",
  "total_prokaryotic_reads",
  "archaeal_reads",
  "archaeal_percent"
)
na_subsoil_summary <- na_subsoil_summary[
  order(na_subsoil_summary$depth_midpoint_cm, na_subsoil_summary$sample_id),
  ,
  drop = FALSE
]
row.names(na_subsoil_summary) <- NULL

write.csv(
  na_subsoil_summary,
  file.path(out_dir, "round41_NA_eDNA_20_60cm_library_audit.csv"),
  row.names = FALSE,
  na = ""
)

provenance_audit <- data.frame(
  metric = c(
    "Raw archaeal ASVs",
    "Pre-rarefaction archaeal ASVs",
    "Raw and pre-rarefaction ASV sets identical",
    "Raw and pre-rarefaction ASV order identical",
    "Shared biological libraries compared",
    "Counts identical across shared biological libraries",
    "Control-detected ASVs retained in pre-rarefaction table"
  ),
  value = c(
    nrow(archaea),
    nrow(pre_rarefaction),
    setequal(archaea$asv_id, rownames(pre_rarefaction)),
    identical(archaea$asv_id, rownames(pre_rarefaction)),
    length(shared_pre_rarefaction_cols),
    identical(unname(raw_shared_counts), unname(pre_rarefaction_shared_counts)),
    paste0(
      sum(asv_audit$asv_id %in% rownames(pre_rarefaction)),
      "/",
      nrow(asv_audit)
    )
  ),
  interpretation = c(
    "Archaeal rows in data/file_A.txt.",
    "Rows in data/arc_unrarefild/ASV_Arc_200cm_delete_less200.txt.",
    "No archaeal ASV is absent from either archived table.",
    "The archived ASV row order is unchanged.",
    "All pre-rarefaction columns are biological libraries.",
    "No count change is present in the archived shared columns.",
    "Every control-detected ASV remains in the archived pre-rarefaction table."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  provenance_audit,
  file.path(out_dir, "round41_raw_to_unrarefied_provenance_audit.csv"),
  row.names = FALSE,
  na = ""
)

control_archaeal_reads <- control_summary$archaeal_reads
retained_archaeal_reads <- colSums(retained_counts)

stopifnot(
  median(control_archaeal_reads) == 14,
  min(control_archaeal_reads) == 0,
  max(control_archaeal_reads) == 323,
  sum(control_archaeal_reads) == 1662,
  median(retained_archaeal_reads) == 1513,
  min(retained_archaeal_reads) == 226,
  max(retained_archaeal_reads) == 22364,
  nrow(asv_audit) == 140L,
  min(na_subsoil_summary$archaeal_reads) == 3840,
  max(na_subsoil_summary$archaeal_reads) == 22364
)

summary_lines <- c(
  "Round41 descriptive negative-control audit",
  "",
  paste0(
    "Negative and no-template controls: n = 27; archaeal reads median = 14; ",
    "range = 0-323; total = 1,662."
  ),
  paste0(
    "Retained biological samples: n = 95; archaeal reads median = 1,513; ",
    "range = 226-22,364."
  ),
  paste0(
    "NA eDNA libraries at 20-60 cm: n = 6; archaeal-read range = ",
    "3,840-22,364."
  ),
  paste0(
    "Archaeal ASVs detected at least once in a negative or no-template ",
    "control: 140."
  ),
  paste0(
    "Archived raw-to-pre-rarefaction comparison: 1,166/1,166 archaeal ASVs ",
    "retained; counts identical across 214 shared biological libraries."
  ),
  "",
  paste0(
    "This audit is descriptive only. It applies no contaminant threshold, ",
    "classifies no ASV as a contaminant, and subtracts no reads or ASVs."
  )
)
writeLines(
  summary_lines,
  file.path(out_dir, "round41_negative_control_audit_summary.txt"),
  useBytes = TRUE
)

message(paste(summary_lines[nzchar(summary_lines)], collapse = "\n"))
