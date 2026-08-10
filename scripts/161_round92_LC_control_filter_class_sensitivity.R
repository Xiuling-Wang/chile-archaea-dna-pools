#!/usr/bin/env Rscript

# Round92 targeted sensitivity for the La Campana (LC) class-composition result.
#
# The rule is inherited unchanged from the Round48 control sensitivity:
# remove every archaeal ASV with at least one read in any negative
# extraction/process or no-template PCR control. This is a deliberately
# conservative sensitivity analysis, not a contaminant classification.
# Counts are not subtracted, thresholds are not tuned, and the primary
# unfiltered composition remains the primary analysis.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  project_root <- dirname(dirname(script_path))
} else {
  project_root <- normalizePath(".", mustWork = TRUE)
}

raw_path <- file.path(project_root, "data", "file_A.txt")
unrarefied_path <- file.path(
  project_root, "data", "arc_unrarefild", "ASV_Arc_200cm_delete_less200.txt"
)
rarefied_path <- file.path(
  project_root, "data", "arc_60cm_202rare", "asv_202.txt"
)
metadata_path <- file.path(project_root, "data", "env_2024.csv")
out_dir <- file.path(project_root, "analysis", "reproducibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.delim(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
unrarefied <- read.delim(
  unrarefied_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rarefied <- read.delim(
  rarefied_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.csv(
  metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = ""
)

taxonomy_cols <- c("asv_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
library_cols <- setdiff(names(raw), taxonomy_cols)
negative_control_cols <- library_cols[
  grepl("(^NTC_|_NC_|^NC_)", library_cols)
]
positive_control_cols <- library_cols[grepl("^PC_", library_cols)]
archaea <- raw[raw$Kingdom == "Archaea", , drop = FALSE]

control_counts <- as.matrix(archaea[, negative_control_cols, drop = FALSE])
storage.mode(control_counts) <- "numeric"
control_detected <- rowSums(control_counts) > 0
control_asv_ids <- archaea$asv_id[control_detected]

retained_ids <- colnames(rarefied)
counts <- as.matrix(unrarefied[, retained_ids, drop = FALSE])
storage.mode(counts) <- "numeric"
filtered_counts <- counts[
  !rownames(counts) %in% control_asv_ids,
  ,
  drop = FALSE
]

meta <- metadata[match(retained_ids, metadata$sample_id), , drop = FALSE]
if (anyNA(meta$sample_id) || !identical(meta$sample_id, retained_ids)) {
  stop("Retained metadata could not be aligned to the 95 unrarefied libraries.")
}
meta$site_display <- ifelse(meta$site == "NB", "NA", meta$site)
meta$depth_interval_cm <- gsub("_", "-", meta$depths_cm, fixed = TRUE)

strip_rank <- function(x) sub("^[a-z]__", "", x)
class_lookup <- setNames(strip_rank(archaea$Class), archaea$asv_id)
nit_ids <- names(class_lookup)[class_lookup == "Nitrososphaeria"]
nit_ids <- intersect(nit_ids, rownames(counts))
nit_control_ids <- intersect(nit_ids, control_asv_ids)

class_percent <- function(matrix_now, taxon_ids) {
  totals <- colSums(matrix_now)
  if (any(totals <= 0)) {
    stop("A retained library has zero archaeal reads after filtering.")
  }
  ids <- intersect(taxon_ids, rownames(matrix_now))
  numerator <- if (length(ids) > 0L) {
    colSums(matrix_now[ids, , drop = FALSE])
  } else {
    setNames(rep(0, ncol(matrix_now)), colnames(matrix_now))
  }
  100 * numerator / totals
}

primary_nit <- class_percent(counts, nit_ids)
filtered_nit <- class_percent(filtered_counts, nit_ids)

sample_results <- data.frame(
  sample_id = retained_ids,
  site_internal = meta$site,
  site = meta$site_display,
  pit = meta$pit,
  depth_interval_cm = meta$depth_interval_cm,
  depth_midpoint_cm = meta$depth,
  dna_pool = meta$dna_type,
  primary_archaeal_reads = as.integer(colSums(counts)),
  filtered_archaeal_reads = as.integer(colSums(filtered_counts)),
  retained_archaeal_reads_percent = 100 * colSums(filtered_counts) / colSums(counts),
  primary_Nitrososphaeria_percent = as.numeric(primary_nit),
  control_filtered_Nitrososphaeria_percent = as.numeric(filtered_nit),
  percentage_point_change_after_filter = as.numeric(filtered_nit - primary_nit),
  stringsAsFactors = FALSE
)

aggregate_mean <- function(df, grouping) {
  formula <- as.formula(
    paste(
      "cbind(primary_Nitrososphaeria_percent,",
      "control_filtered_Nitrososphaeria_percent,",
      "retained_archaeal_reads_percent) ~",
      paste(grouping, collapse = " + ")
    )
  )
  aggregate(formula, data = df, FUN = mean)
}

all_site_pool <- aggregate_mean(sample_results, c("site", "dna_pool"))
all_site_pool$n_libraries <- as.integer(
  ave(
    rep(1L, nrow(all_site_pool)),
    interaction(all_site_pool$site, all_site_pool$dna_pool),
    FUN = length
  )
)

pool_gap_table <- do.call(
  rbind,
  lapply(c("AZ", "SG", "LC", "NA"), function(site_now) {
    x <- all_site_pool[all_site_pool$site == site_now, , drop = FALSE]
    i <- x[x$dna_pool == "iDNA", , drop = FALSE]
    e <- x[x$dna_pool == "eDNA", , drop = FALSE]
    data.frame(
      site = site_now,
      iDNA_n = sum(sample_results$site == site_now & sample_results$dna_pool == "iDNA"),
      eDNA_n = sum(sample_results$site == site_now & sample_results$dna_pool == "eDNA"),
      primary_iDNA_Nitrososphaeria_percent = i$primary_Nitrososphaeria_percent,
      primary_eDNA_Nitrososphaeria_percent = e$primary_Nitrososphaeria_percent,
      primary_iDNA_minus_eDNA_percentage_points =
        i$primary_Nitrososphaeria_percent - e$primary_Nitrososphaeria_percent,
      control_filtered_iDNA_Nitrososphaeria_percent =
        i$control_filtered_Nitrososphaeria_percent,
      control_filtered_eDNA_Nitrososphaeria_percent =
        e$control_filtered_Nitrososphaeria_percent,
      control_filtered_iDNA_minus_eDNA_percentage_points =
        i$control_filtered_Nitrososphaeria_percent -
        e$control_filtered_Nitrososphaeria_percent,
      stringsAsFactors = FALSE
    )
  })
)

lc <- sample_results[sample_results$site == "LC", , drop = FALSE]
lc_depth <- aggregate_mean(
  lc,
  c("depth_midpoint_cm", "depth_interval_cm", "dna_pool")
)
lc_depth$n_libraries <- vapply(
  seq_len(nrow(lc_depth)),
  function(i) {
    sum(
      lc$depth_midpoint_cm == lc_depth$depth_midpoint_cm[i] &
        lc$dna_pool == lc_depth$dna_pool[i]
    )
  },
  integer(1)
)
lc_depth <- lc_depth[
  order(lc_depth$depth_midpoint_cm, match(lc_depth$dna_pool, c("iDNA", "eDNA"))),
  ,
  drop = FALSE
]

lc_pair <- merge(
  lc[lc$dna_pool == "iDNA", ],
  lc[lc$dna_pool == "eDNA", ],
  by = c("site_internal", "site", "pit", "depth_interval_cm", "depth_midpoint_cm"),
  suffixes = c("_iDNA", "_eDNA")
)
lc_pair$primary_iDNA_minus_eDNA_percentage_points <-
  lc_pair$primary_Nitrososphaeria_percent_iDNA -
  lc_pair$primary_Nitrososphaeria_percent_eDNA
lc_pair$control_filtered_iDNA_minus_eDNA_percentage_points <-
  lc_pair$control_filtered_Nitrososphaeria_percent_iDNA -
  lc_pair$control_filtered_Nitrososphaeria_percent_eDNA

control_taxon_audit <- data.frame(
  metric = c(
    "negative_or_no_template_controls_n",
    "positive_control_libraries_excluded_n",
    "control_detected_archaeal_ASVs_n",
    "Nitrososphaeria_ASVs_in_unrarefied_table_n",
    "control_detected_Nitrososphaeria_ASVs_n",
    "control_reads_from_Nitrososphaeria_ASVs",
    "retained_95_reads_from_control_detected_Nitrososphaeria_ASVs"
  ),
  value = c(
    length(negative_control_cols),
    length(positive_control_cols),
    length(control_asv_ids),
    length(nit_ids),
    length(nit_control_ids),
    sum(control_counts[match(nit_control_ids, archaea$asv_id), , drop = FALSE]),
    sum(counts[intersect(nit_control_ids, rownames(counts)), , drop = FALSE])
  ),
  interpretation = c(
    "Negative extraction/process and no-template PCR controls used by the fixed sensitivity rule.",
    "PC_1 and PC_2 libraries across three sequencing-library batches; excluded from all biological summaries.",
    "Every ASV with at least one read in a negative or no-template control is removed in the sensitivity analysis.",
    "ASVs assigned to the class Nitrososphaeria in the 95-library unrarefied table.",
    "Nitrososphaeria ASVs removed by the conservative sensitivity rule.",
    "Reads from the removed Nitrososphaeria ASVs across the 27 negative/no-template controls.",
    "Reads from the removed Nitrososphaeria ASVs across the 95 retained biological libraries."
  ),
  stringsAsFactors = FALSE
)

positive_control_audit <- data.frame(
  positive_control_id = positive_control_cols,
  positive_control_label = sub("_lib[0-9]+$", "", positive_control_cols),
  library_batch = sub("^.*_(lib[0-9]+)$", "\\1", positive_control_cols),
  total_prokaryotic_reads = as.integer(
    colSums(raw[, positive_control_cols, drop = FALSE])
  ),
  archaeal_reads = as.integer(
    colSums(archaea[, positive_control_cols, drop = FALSE])
  ),
  archaeal_percent = 100 * as.integer(
    colSums(archaea[, positive_control_cols, drop = FALSE])
  ) / as.integer(colSums(raw[, positive_control_cols, drop = FALSE])),
  archived_expected_composition_available = FALSE,
  used_for_contaminant_classification = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  length(retained_ids) == 95L,
  length(negative_control_cols) == 27L,
  length(positive_control_cols) == 6L,
  length(control_asv_ids) == 140L,
  all(positive_control_cols == c(
    "PC_1_lib1", "PC_2_lib1", "PC_1_lib2",
    "PC_2_lib2", "PC_1_lib3", "PC_2_lib3"
  )),
  nrow(lc) > 0L,
  nrow(lc_pair) > 0L
)

lc_gap <- pool_gap_table[pool_gap_table$site == "LC", , drop = FALSE]
stopifnot(
  abs(lc_gap$primary_iDNA_Nitrososphaeria_percent - 87.9) < 0.15,
  abs(lc_gap$primary_eDNA_Nitrososphaeria_percent - 53.1) < 0.15,
  abs(lc_gap$primary_iDNA_minus_eDNA_percentage_points - 34.9) < 0.15
)

write.csv(
  sample_results,
  file.path(out_dir, "round92_LC_control_filter_class_by_sample.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  pool_gap_table,
  file.path(out_dir, "round92_LC_control_filter_class_all_sites.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  lc_depth,
  file.path(out_dir, "round92_LC_control_filter_class_by_depth.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  lc_pair,
  file.path(out_dir, "round92_LC_control_filter_class_paired_LC.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  control_taxon_audit,
  file.path(out_dir, "round92_LC_control_filter_taxon_audit.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  positive_control_audit,
  file.path(out_dir, "round92_positive_control_library_audit.csv"),
  row.names = FALSE,
  na = ""
)

filtered_largest <- pool_gap_table$site[
  which.max(abs(pool_gap_table$control_filtered_iDNA_minus_eDNA_percentage_points))
]
filtered_lc_positive_pairs <- sum(
  lc_pair$control_filtered_iDNA_minus_eDNA_percentage_points > 0
)
summary_lines <- c(
  "Round92 targeted LC class-composition control-filter sensitivity",
  "",
  sprintf(
    paste0(
      "Primary LC Nitrososphaeria means: iDNA %.3f%%; eDNA %.3f%%; ",
      "iDNA-eDNA gap %.3f percentage points."
    ),
    lc_gap$primary_iDNA_Nitrososphaeria_percent,
    lc_gap$primary_eDNA_Nitrososphaeria_percent,
    lc_gap$primary_iDNA_minus_eDNA_percentage_points
  ),
  sprintf(
    paste0(
      "Control-filtered LC Nitrososphaeria means: iDNA %.3f%%; eDNA %.3f%%; ",
      "iDNA-eDNA gap %.3f percentage points."
    ),
    lc_gap$control_filtered_iDNA_Nitrososphaeria_percent,
    lc_gap$control_filtered_eDNA_Nitrososphaeria_percent,
    lc_gap$control_filtered_iDNA_minus_eDNA_percentage_points
  ),
  sprintf(
    "LC remained the largest absolute site-level pool gap after filtering: %s.",
    ifelse(filtered_largest == "LC", "yes", paste0("no; largest was ", filtered_largest))
  ),
  sprintf(
    "Matched LC horizons retaining iDNA > eDNA after filtering: %d/%d.",
    filtered_lc_positive_pairs,
    nrow(lc_pair)
  ),
  sprintf(
    "Control-detected archaeal ASVs removed: %d, including %d Nitrososphaeria ASVs.",
    length(control_asv_ids),
    length(nit_control_ids)
  ),
  sprintf(
    "Positive-control libraries present and excluded from biological summaries: %d (%s).",
    length(positive_control_cols),
    paste(positive_control_cols, collapse = ", ")
  ),
  "",
  paste(
    "Interpretation boundary: this fixed worst-case removal is a sensitivity analysis,",
    "not a contaminant classification, and the positive-control identities do not",
    "establish their expected composition or acceptance criteria."
  )
)
writeLines(
  summary_lines,
  file.path(out_dir, "round92_LC_control_filter_class_summary.txt"),
  useBytes = TRUE
)

message(paste(summary_lines[nzchar(summary_lines)], collapse = "\n"))
