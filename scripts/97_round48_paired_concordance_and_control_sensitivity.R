#!/usr/bin/env Rscript

# Round48 canonical paired-pool analysis and control-informed sensitivity audit.
#
# Prespecified control rule (fixed before this script was run): remove every
# archaeal ASV with at least one read in any negative extraction/process or
# no-template PCR control. This is a deliberately conservative sensitivity
# filter, not a contaminant classification and not a replacement primary
# pipeline. No abundance subtraction or tuned threshold is used.

suppressPackageStartupMessages({
  library(dplyr)
  library(permute)
  library(purrr)
  library(readr)
  library(tidyr)
  library(vegan)
})

N_PERMUTATIONS <- 9999L
CONCORDANCE_SEED <- 20260605L
PAIRED_PERMANOVA_SEED <- 20260710L
SITE_TEST_SEED <- 20260713L
SITE_LEVELS <- c("AZ", "SG", "LC", "NB")
DNA_LEVELS <- c("iDNA", "eDNA")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  project_root <- dirname(dirname(script_path))
} else {
  project_root <- normalizePath(".", mustWork = TRUE)
}

paired_out_dir <- file.path(project_root, "analysis", "paired_pool_similarity")
sensitivity_out_dir <- file.path(project_root, "analysis", "reproducibility")
dir.create(paired_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sensitivity_out_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  rarefied_202 = file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt"),
  unrarefied = file.path(
    project_root, "data", "arc_unrarefild", "ASV_Arc_200cm_delete_less200.txt"
  ),
  taxonomy = file.path(project_root, "data", "arc_100_202rarefild", "tax_202.txt"),
  metadata = file.path(project_root, "data", "env_2024.csv"),
  raw_all_libraries = file.path(project_root, "data", "file_A.txt"),
  fig1_read_fractions = file.path(
    project_root, "analysis", "reframed_figures",
    "Fig1_archaeal_read_fraction_source.csv"
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
  if (anyNA(x) || any(x < 0)) {
    stop("Count table contains missing or negative values: ", path)
  }
  x
}

display_site <- function(x) ifelse(as.character(x) == "NB", "NA", as.character(x))

rare_counts <- read_count_table(input_paths[["rarefied_202"]])
unrarefied_counts <- read_count_table(input_paths[["unrarefied"]])
taxonomy <- read.delim(
  input_paths[["taxonomy"]],
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.csv(
  input_paths[["metadata"]],
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = ""
)
raw_all <- read.delim(
  input_paths[["raw_all_libraries"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fig1_reads <- read.csv(
  input_paths[["fig1_read_fractions"]],
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = ""
)

if (!all(colSums(rare_counts) == 202)) {
  stop("The canonical paired input is not uniformly rarefied to 202 reads.")
}
if (!all(colnames(rare_counts) %in% colnames(unrarefied_counts))) {
  stop("Not every rarefied sample is represented in the unrarefied table.")
}
if (!all(colnames(rare_counts) %in% metadata$sample_id)) {
  stop("Not every rarefied sample is represented in metadata.")
}

retained_sample_ids <- colnames(rare_counts)
unrarefied_95 <- unrarefied_counts[, retained_sample_ids, drop = FALSE]
metadata_95 <- metadata[match(retained_sample_ids, metadata$sample_id), , drop = FALSE]
if (!identical(metadata_95$sample_id, retained_sample_ids)) {
  stop("Metadata could not be aligned to the canonical 95 retained libraries.")
}

build_pair_map <- function(sample_ids, metadata_all) {
  meta <- metadata_all[match(sample_ids, metadata_all$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id)) stop("Pair metadata contains unmatched sample IDs.")

  meta <- meta %>%
    filter(depth < 60, dna_type %in% DNA_LEVELS) %>%
    mutate(
      site_internal = as.character(site),
      site_display = display_site(site),
      dna_type = factor(dna_type, levels = DNA_LEVELS),
      depth_order = as.numeric(depth_order),
      depth_cm = gsub("_", "-", depths_cm, fixed = TRUE),
      pair_id = paste(site_internal, pit, depths_cm, sep = "_")
    )

  pair_counts <- meta %>%
    count(
      pair_id, site_internal, site_display, pit, depth_order, depth_cm,
      dna_type, name = "n"
    ) %>%
    pivot_wider(names_from = dna_type, values_from = n, values_fill = 0)

  complete_ids <- pair_counts %>%
    filter(iDNA == 1L, eDNA == 1L) %>%
    pull(pair_id)

  paired_meta <- meta %>%
    filter(pair_id %in% complete_ids) %>%
    arrange(pair_id, dna_type)

  if (!all(table(paired_meta$pair_id) == 2L)) {
    stop("Every retained pair must contain exactly one iDNA and one eDNA library.")
  }

  pair_map <- paired_meta %>%
    select(
      pair_id, site_internal, site_display, pit, depth_order, depth_cm,
      dna_type, sample_id
    ) %>%
    pivot_wider(names_from = dna_type, values_from = sample_id, names_glue = "{dna_type}_sample") %>%
    arrange(factor(site_internal, levels = SITE_LEVELS), pit, depth_order)

  list(pair_map = pair_map, paired_meta = paired_meta)
}

run_concordance <- function(counts, metadata_all, scenario) {
  paired <- build_pair_map(colnames(counts), metadata_all)
  pair_map <- paired$pair_map

  paired_ids <- c(pair_map$iDNA_sample, pair_map$eDNA_sample)
  paired_counts <- counts[, paired_ids, drop = FALSE]
  paired_counts <- paired_counts[rowSums(paired_counts) > 0, , drop = FALSE]

  i_mat <- t(paired_counts[, pair_map$iDNA_sample, drop = FALSE])
  e_mat <- t(paired_counts[, pair_map$eDNA_sample, drop = FALSE])
  rownames(i_mat) <- pair_map$pair_id
  rownames(e_mat) <- pair_map$pair_id
  if (any(rowSums(i_mat) == 0) || any(rowSums(e_mat) == 0)) {
    stop("A retained paired library has zero reads in scenario: ", scenario)
  }

  i_hellinger <- decostand(i_mat, method = "hellinger")
  e_hellinger <- decostand(e_mat, method = "hellinger")
  i_dist <- vegdist(i_hellinger, method = "bray")
  e_dist <- vegdist(e_hellinger, method = "bray")
  i_pcoa <- wcmdscale(i_dist, k = 2, eig = TRUE)
  e_pcoa <- wcmdscale(e_dist, k = 2, eig = TRUE)

  set.seed(CONCORDANCE_SEED)
  protest_all <- protest(i_pcoa$points, e_pcoa$points, permutations = N_PERMUTATIONS)
  set.seed(CONCORDANCE_SEED)
  mantel_all <- mantel(
    i_dist, e_dist, method = "spearman", permutations = N_PERMUTATIONS
  )

  overall <- tibble(
    scenario = scenario,
    method = c("PROTEST Procrustes", "Mantel Spearman"),
    statistic = c(unname(protest_all$t0), unname(mantel_all$statistic)),
    p_value = c(protest_all$signif, mantel_all$signif),
    n_complete_pairs = nrow(pair_map),
    n_asvs_analyzed = nrow(paired_counts),
    transformation = "Hellinger",
    dissimilarity = "Bray-Curtis",
    ordination_dimensions = c(2L, NA_integer_),
    permutations = N_PERMUTATIONS,
    seed = CONCORDANCE_SEED
  )

  by_site <- map_dfr(SITE_LEVELS, function(site_now) {
    ids <- pair_map$pair_id[pair_map$site_internal == site_now]
    i_site <- i_hellinger[ids, , drop = FALSE]
    e_site <- e_hellinger[ids, , drop = FALSE]
    i_site_dist <- vegdist(i_site, method = "bray")
    e_site_dist <- vegdist(e_site, method = "bray")
    i_site_pcoa <- wcmdscale(i_site_dist, k = 2, eig = TRUE)
    e_site_pcoa <- wcmdscale(e_site_dist, k = 2, eig = TRUE)

    set.seed(CONCORDANCE_SEED)
    protest_site <- protest(
      i_site_pcoa$points, e_site_pcoa$points, permutations = N_PERMUTATIONS
    )
    set.seed(CONCORDANCE_SEED)
    mantel_site <- mantel(
      i_site_dist, e_site_dist,
      method = "spearman", permutations = N_PERMUTATIONS
    )

    tibble(
      scenario = scenario,
      site_internal = site_now,
      site_display = display_site(site_now),
      n_complete_pairs = length(ids),
      n_asvs_analyzed = nrow(paired_counts),
      protest_r = unname(protest_site$t0),
      protest_p = protest_site$signif,
      mantel_r = unname(mantel_site$statistic),
      mantel_p = mantel_site$signif,
      transformation = "Hellinger",
      dissimilarity = "Bray-Curtis",
      ordination_dimensions = 2L,
      permutations = N_PERMUTATIONS,
      seed = CONCORDANCE_SEED
    )
  })

  list(overall = overall, by_site = by_site, pair_map = pair_map)
}

run_paired_permanova <- function(counts, metadata_all) {
  paired <- build_pair_map(colnames(counts), metadata_all)
  paired_meta <- paired$paired_meta %>%
    mutate(
      pair_id = droplevels(factor(pair_id)),
      dna_type = factor(dna_type, levels = DNA_LEVELS)
    ) %>%
    arrange(pair_id, dna_type)

  paired_counts <- counts[, paired_meta$sample_id, drop = FALSE]
  paired_counts <- paired_counts[rowSums(paired_counts) > 0, , drop = FALSE]
  paired_hellinger <- decostand(t(paired_counts), method = "hellinger")

  # Estimand: the average global fraction-associated compositional shift after
  # absorbing all between-pair differences with a site x pit x depth pair
  # fixed effect. Under the null, iDNA/eDNA labels may swap only within pairs.
  restricted_control <- how(
    within = Within(type = "free"),
    blocks = paired_meta$pair_id,
    nperm = N_PERMUTATIONS
  )

  set.seed(PAIRED_PERMANOVA_SEED)
  fit <- adonis2(
    paired_hellinger ~ pair_id + dna_type,
    data = paired_meta,
    method = "bray",
    permutations = restricted_control,
    by = "margin"
  )

  fit_df <- as.data.frame(fit) %>% tibble::rownames_to_column("term")
  fraction_row <- fit_df %>% filter(term == "dna_type")
  residual_row <- fit_df %>% filter(term == "Residual")
  if (nrow(fraction_row) != 1L) stop("Could not extract the canonical fraction term.")

  tibble(
    test = "global paired iDNA-eDNA fraction effect",
    estimand = paste(
      "average fraction-associated shift conditional on site x pit x depth pair",
      "identity"
    ),
    model = "Hellinger ASV composition ~ pair_id + dna_type",
    n_complete_pairs = nlevels(paired_meta$pair_id),
    n_libraries = nrow(paired_meta),
    n_asvs_analyzed = nrow(paired_counts),
    numerator_df = fraction_row$Df,
    denominator_df = residual_row$Df,
    sum_of_squares = fraction_row$SumOfSqs,
    F = fraction_row$F,
    R2 = fraction_row$R2,
    p_value = fraction_row$`Pr(>F)`,
    transformation = "Hellinger",
    dissimilarity = "Bray-Curtis",
    permutation_scheme = paste(
      "free iDNA/eDNA label exchange only within site x pit x depth pairs;",
      "pair identity included as a fixed effect"
    ),
    permutations = N_PERMUTATIONS,
    seed = PAIRED_PERMANOVA_SEED
  )
}

profile_asv_composition <- function(counts, metadata_all) {
  meta <- metadata_all[match(colnames(counts), metadata_all$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id)) stop("Profile analysis contains unmatched sample IDs.")
  if (any(colSums(counts) == 0)) stop("Profile analysis contains a zero-sum library.")

  rel <- sweep(counts, 2, colSums(counts), "/")
  sample_df <- data.frame(
    sample_id = colnames(counts),
    t(rel),
    check.names = FALSE
  ) %>%
    left_join(meta %>% select(sample_id, site, pit, dna_type), by = "sample_id")

  pool_means <- sample_df %>%
    group_by(site, pit, dna_type) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  pool_means %>%
    group_by(site, pit) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
}

class_profile_composition <- function(counts, metadata_all, taxonomy_table) {
  strip_rank <- function(x) sub("^[a-z]__", "", x)
  class_lookup <- setNames(strip_rank(taxonomy_table$Class), rownames(taxonomy_table))
  rel <- sweep(counts, 2, colSums(counts), "/")
  classes <- sort(unique(class_lookup[intersect(names(class_lookup), rownames(counts))]))
  classes <- classes[!is.na(classes)]

  values <- vapply(classes, function(class_name) {
    ids <- intersect(names(class_lookup)[class_lookup == class_name], rownames(counts))
    if (length(ids) == 0L) return(rep(0, ncol(counts)))
    colSums(rel[ids, , drop = FALSE])
  }, numeric(ncol(counts)))
  rownames(values) <- colnames(counts)

  meta <- metadata_all[match(rownames(values), metadata_all$sample_id), , drop = FALSE]
  sample_df <- data.frame(
    sample_id = rownames(values),
    values,
    check.names = FALSE
  ) %>%
    left_join(meta %>% select(sample_id, site, pit, dna_type), by = "sample_id")

  pool_means <- sample_df %>%
    group_by(site, pit, dna_type) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  pool_means %>%
    group_by(site, pit) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
}

run_site_test <- function(counts, metadata_all, permutation_matrix, scenario) {
  profile <- profile_asv_composition(counts, metadata_all)
  asv_cols <- setdiff(names(profile), c("site", "pit"))
  asv_cols <- asv_cols[colSums(profile[, asv_cols, drop = FALSE]) > 0]
  profile_dist <- vegdist(
    profile[, asv_cols, drop = FALSE], method = "robust.aitchison"
  )
  site_factor <- factor(profile$site, levels = SITE_LEVELS)
  fit <- adonis2(profile_dist ~ site_factor, permutations = permutation_matrix)

  tibble(
    scenario = scenario,
    n_libraries = ncol(counts),
    n_independent_profiles = nrow(profile),
    n_asvs_analyzed = length(asv_cols),
    R2 = unname(fit$R2[1]),
    F = unname(fit$F[1]),
    p_value = unname(fit$`Pr(>F)`[1]),
    distance = "robust Aitchison",
    aggregation = paste(
      "sample relative abundance -> mean within site x pit x fraction ->",
      "mean across fractions"
    ),
    permutations = N_PERMUTATIONS,
    seed_stream = paste0(
      SITE_TEST_SEED,
      " after archived class-level test; fixed matrix shared across scenarios"
    )
  )
}

# Identify the control-detected ASVs directly from the raw table.
taxonomy_cols <- c("asv_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
library_cols <- setdiff(names(raw_all), taxonomy_cols)
negative_control_cols <- library_cols[grepl("(^NTC_|_NC_|^NC_)", library_cols)]
positive_control_cols <- library_cols[grepl("^PC_", library_cols)]
archaea_raw <- raw_all[raw_all$Kingdom == "Archaea", , drop = FALSE]
control_counts <- as.matrix(archaea_raw[, negative_control_cols, drop = FALSE])
storage.mode(control_counts) <- "numeric"
control_detected <- rowSums(control_counts) > 0
control_asv_ids <- archaea_raw$asv_id[control_detected]

if (length(negative_control_cols) != 27L || length(control_asv_ids) != 140L) {
  stop("Control audit no longer matches 27 controls and 140 detected archaeal ASVs.")
}
if (!all(control_asv_ids %in% rownames(unrarefied_counts))) {
  stop("A control-detected ASV is absent from the unrarefied biological table.")
}

control_filter_asvs <- archaea_raw[control_detected, taxonomy_cols, drop = FALSE] %>%
  mutate(
    control_reads_total = as.integer(rowSums(control_counts[control_detected, , drop = FALSE])),
    controls_detected_n = as.integer(rowSums(control_counts[control_detected, , drop = FALSE] > 0)),
    in_unrarefied_95 = asv_id %in% rownames(unrarefied_95),
    in_rarefied_202 = asv_id %in% rownames(rare_counts)
  ) %>%
  arrange(desc(control_reads_total), desc(controls_detected_n), asv_id)

control_rule <- tibble(
  sensitivity_only = TRUE,
  rule = paste(
    "remove every archaeal ASV with at least one read in any negative",
    "extraction/process or no-template PCR control"
  ),
  justification = paste(
    "transparent worst-case filter because no preserved batch-aware decontam model,",
    "DNA concentration record, or prespecified dominance threshold exists"
  ),
  negative_and_ntc_controls_n = length(negative_control_cols),
  positive_controls_excluded_n = length(positive_control_cols),
  control_detected_archaeal_asvs_n = length(control_asv_ids),
  control_detected_asvs_in_unrarefied_95_n = sum(
    control_asv_ids %in% rownames(unrarefied_95)
  ),
  control_detected_asvs_in_rarefied_202_n = sum(
    control_asv_ids %in% rownames(rare_counts)
  ),
  subtraction = "none",
  tuned_threshold = "none",
  local_percentage_denominator = paste(
    "original total prokaryotic reads; only the archaeal numerator is filtered",
    "(conservative)"
  ),
  primary_pipeline_changed = FALSE
)

write_csv(
  control_filter_asvs,
  file.path(sensitivity_out_dir, "round48_control_filter_ASVs.csv")
)
write_csv(
  control_rule,
  file.path(sensitivity_out_dir, "round48_control_filter_rule.csv")
)

# Canonical paired analysis on the established 202-read table.
canonical_concordance <- run_concordance(
  rare_counts, metadata, "primary_unfiltered"
)
canonical_permanova <- run_paired_permanova(rare_counts, metadata)

write_csv(
  canonical_concordance$pair_map,
  file.path(paired_out_dir, "round48_complete_pair_map.csv")
)
write_csv(
  canonical_concordance$overall,
  file.path(paired_out_dir, "round48_paired_concordance_overall.csv")
)
write_csv(
  canonical_concordance$by_site,
  file.path(paired_out_dir, "round48_paired_concordance_by_site.csv")
)
write_csv(
  canonical_permanova,
  file.path(paired_out_dir, "round48_paired_pool_permanova.csv")
)

# Label every conflicting or invalid legacy paired output without deleting it.
superseded_outputs <- tribble(
  ~legacy_path, ~round48_status, ~replacement, ~reason,
  "analysis/paired_pool_similarity/paired_PERMANOVA_pool_effect_strata_pair.txt",
  "SUPERSEDED_NONCANONICAL_SPECIFICATION",
  "analysis/paired_pool_similarity/round48_paired_pool_permanova.csv",
  paste(
    "restricted labels but no pair fixed effect; its F=0.602 denominator includes",
    "between-pair variation rather than the within-pair residual"
  ),
  "analysis/reproducibility/round25_pair_restricted_pool_permanova.csv",
  "SUPERSEDED_FILE_SPECIFICATION_RETAINED",
  "analysis/paired_pool_similarity/round48_paired_pool_permanova.csv",
  "pair-fixed estimand retained and cleanly recomputed in the Round48 canonical script",
  "analysis/reproducibility/round25_pair_restricted_pool_permanova.txt",
  "SUPERSEDED_FILE_SPECIFICATION_RETAINED",
  "analysis/paired_pool_similarity/round48_paired_pool_permanova.csv",
  "pair-fixed estimand retained and cleanly recomputed in the Round48 canonical script",
  "analysis/paired_pool_similarity/procrustes_mantel_overall.csv",
  "SUPERSEDED_BY_ROUND48_RECOMPUTATION",
  "analysis/paired_pool_similarity/round48_paired_concordance_overall.csv",
  "legacy values retained only for provenance",
  "analysis/paired_pool_similarity/procrustes_mantel_by_site.csv",
  "SUPERSEDED_BY_ROUND48_RECOMPUTATION",
  "analysis/paired_pool_similarity/round48_paired_concordance_by_site.csv",
  "legacy values retained only for provenance",
  "analysis/paired_pool_similarity/paired_vs_unmatched_same_site_depth_tests.csv",
  "SUPERSEDED_EXCLUDED_FROM_INFERENCE",
  "none",
  "unmatched combinations reuse samples and do not provide independent-pair inference",
  "analysis/paired_pool_similarity/unmatched_same_site_depth_similarity.csv",
  "SUPERSEDED_EXCLUDED_FROM_INFERENCE",
  "none",
  "source table for the pseudoreplicated unmatched comparison; no Round48 analogue",
  "analysis/paired_pool_similarity/paired_vs_unmatched_same_site_depth.pdf",
  "SUPERSEDED_EXCLUDED_FROM_INFERENCE",
  "none",
  "visualization of the pseudoreplicated unmatched comparison; no Round48 analogue"
)
write_csv(
  superseded_outputs,
  file.path(paired_out_dir, "round48_superseded_outputs.csv")
)

# Reproduce the archived site-test permutation stream once, then share the
# resulting fixed permutation matrix across all sensitivity scenarios.
original_profile <- profile_asv_composition(unrarefied_95, metadata_95)
class_profile <- class_profile_composition(unrarefied_95, metadata_95, taxonomy)
class_cols <- setdiff(names(class_profile), c("site", "pit"))
class_cols <- class_cols[colSums(class_profile[, class_cols, drop = FALSE]) > 0]
class_dist <- vegdist(
  class_profile[, class_cols, drop = FALSE], method = "robust.aitchison"
)
class_site_factor <- factor(class_profile$site, levels = SITE_LEVELS)
set.seed(SITE_TEST_SEED)
invisible(adonis2(class_dist ~ class_site_factor, permutations = N_PERMUTATIONS))
site_permutation_matrix <- shuffleSet(
  nrow(original_profile),
  nset = N_PERMUTATIONS,
  control = how(nperm = N_PERMUTATIONS)
)

# The smallest retained biological archaeal library is defined from the full
# library counts used for Fig. 1, not from the 202-read table.
smallest_index <- which(fig1_reads$archaeal_reads == min(fig1_reads$archaeal_reads))
if (length(smallest_index) != 1L) stop("The smallest retained library is not unique.")
smallest_sample <- fig1_reads$sample_id[smallest_index]
if (fig1_reads$archaeal_reads[smallest_index] != 226L) {
  stop("The smallest retained biological library is no longer 226 archaeal reads.")
}
smallest_meta <- metadata[match(smallest_sample, metadata$sample_id), , drop = FALSE]

filtered_unrarefied_95 <- unrarefied_95[
  !rownames(unrarefied_95) %in% control_asv_ids, , drop = FALSE
]
filtered_rare_counts <- rare_counts[
  !rownames(rare_counts) %in% control_asv_ids, , drop = FALSE
]
if (any(colSums(filtered_unrarefied_95) == 0) || any(colSums(filtered_rare_counts) == 0)) {
  stop("The control filter creates a zero-sum retained library.")
}

control_filter_library_retention <- metadata_95 %>%
  transmute(
    sample_id,
    site_internal = site,
    site_display = display_site(site),
    pit,
    dna_pool = dna_type,
    depth_cm = gsub("_", "-", depths_cm, fixed = TRUE),
    rarefied_reads_before_filter = as.integer(colSums(rare_counts)),
    rarefied_reads_after_filter = as.integer(colSums(filtered_rare_counts)),
    rarefied_reads_retained_percent = 100 * rarefied_reads_after_filter /
      rarefied_reads_before_filter,
    unrarefied_archaeal_reads_before_filter = as.integer(colSums(unrarefied_95)),
    unrarefied_archaeal_reads_after_filter = as.integer(
      colSums(filtered_unrarefied_95)
    ),
    unrarefied_archaeal_reads_retained_percent =
      100 * unrarefied_archaeal_reads_after_filter /
      unrarefied_archaeal_reads_before_filter
  ) %>%
  arrange(rarefied_reads_after_filter, sample_id)
write_csv(
  control_filter_library_retention,
  file.path(sensitivity_out_dir, "round48_control_filter_library_retention.csv")
)

scenario_counts <- list(
  primary_unfiltered = list(raw = unrarefied_95, rare = rare_counts),
  control_filtered = list(raw = filtered_unrarefied_95, rare = filtered_rare_counts),
  drop_226_library = list(
    raw = unrarefied_95[, colnames(unrarefied_95) != smallest_sample, drop = FALSE],
    rare = rare_counts[, colnames(rare_counts) != smallest_sample, drop = FALSE]
  ),
  control_filtered_plus_drop_226 = list(
    raw = filtered_unrarefied_95[
      , colnames(filtered_unrarefied_95) != smallest_sample, drop = FALSE
    ],
    rare = filtered_rare_counts[
      , colnames(filtered_rare_counts) != smallest_sample, drop = FALSE
    ]
  )
)

site_sensitivity <- imap_dfr(scenario_counts, function(x, scenario_name) {
  run_site_test(x$raw, metadata, site_permutation_matrix, scenario_name)
})

concordance_sensitivity <- imap(scenario_counts, function(x, scenario_name) {
  run_concordance(x$rare, metadata, scenario_name)
})
concordance_sensitivity_overall <- map_dfr(concordance_sensitivity, "overall")
concordance_sensitivity_by_site <- map_dfr(concordance_sensitivity, "by_site")

write_csv(
  site_sensitivity,
  file.path(sensitivity_out_dir, "round48_control_sensitivity_site_test.csv")
)
write_csv(
  concordance_sensitivity_overall,
  file.path(sensitivity_out_dir, "round48_control_sensitivity_concordance_overall.csv")
)
write_csv(
  concordance_sensitivity_by_site,
  file.path(sensitivity_out_dir, "round48_control_sensitivity_concordance_by_site.csv")
)

# Local NA subsoil percentages. The sensitivity numerator excludes the 140
# flagged ASVs; the denominator remains the original total prokaryotic reads.
local_meta <- fig1_reads %>%
  filter(site == "NA", depths_cm_re %in% c("20-40", "40-60")) %>%
  mutate(
    pit = metadata$pit[match(sample_id, metadata$sample_id)],
    pair_id = paste("NB", pit, gsub("-", "_", depths_cm_re, fixed = TRUE), sep = "_")
  )

raw_archaeal_matrix <- as.matrix(archaea_raw[, local_meta$sample_id, drop = FALSE])
storage.mode(raw_archaeal_matrix) <- "numeric"
if (!identical(
  as.numeric(colSums(raw_archaeal_matrix)),
  as.numeric(local_meta$archaeal_reads)
)) {
  stop("Raw archaeal reads do not reproduce the Fig. 1 source table.")
}

make_local_detail <- function(scenario_name, apply_control_filter, drop_sample) {
  asv_keep <- if (apply_control_filter) {
    !archaea_raw$asv_id %in% control_asv_ids
  } else {
    rep(TRUE, nrow(archaea_raw))
  }
  meta_now <- local_meta
  if (drop_sample) meta_now <- meta_now %>% filter(sample_id != smallest_sample)
  counts_now <- as.matrix(archaea_raw[asv_keep, meta_now$sample_id, drop = FALSE])
  storage.mode(counts_now) <- "numeric"
  original_now <- as.matrix(archaea_raw[, meta_now$sample_id, drop = FALSE])
  storage.mode(original_now) <- "numeric"

  meta_now %>%
    transmute(
      scenario = scenario_name,
      sample_id,
      site,
      pit,
      pair_id,
      depth_cm = depths_cm_re,
      dna_pool = dna_type,
      original_total_prokaryotic_reads = total_reads,
      original_archaeal_reads = as.integer(colSums(original_now)),
      removed_control_detected_asv_reads = as.integer(
        colSums(original_now) - colSums(counts_now)
      ),
      retained_archaeal_reads = as.integer(colSums(counts_now)),
      archaeal_percent = 100 * retained_archaeal_reads /
        original_total_prokaryotic_reads
    )
}

local_detail <- bind_rows(
  make_local_detail("primary_unfiltered", FALSE, FALSE),
  make_local_detail("control_filtered", TRUE, FALSE),
  make_local_detail("drop_226_library", FALSE, TRUE),
  make_local_detail("control_filtered_plus_drop_226", TRUE, TRUE)
)

local_pool_summary <- local_detail %>%
  group_by(scenario, depth_cm, dna_pool) %>%
  summarise(
    n_profiles = n(),
    mean_percent = mean(archaeal_percent),
    sd_percent = sd(archaeal_percent),
    minimum_percent = min(archaeal_percent),
    maximum_percent = max(archaeal_percent),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = dna_pool,
    values_from = c(n_profiles, mean_percent, sd_percent, minimum_percent, maximum_percent),
    names_glue = "{dna_pool}_{.value}"
  )

local_pair_summary <- local_detail %>%
  select(scenario, depth_cm, pit, dna_pool, archaeal_percent) %>%
  pivot_wider(names_from = dna_pool, values_from = archaeal_percent) %>%
  group_by(scenario, depth_cm) %>%
  summarise(
    n_complete_pairs = sum(!is.na(eDNA) & !is.na(iDNA)),
    eDNA_higher_pairs_n = sum(eDNA > iDNA, na.rm = TRUE),
    mean_paired_difference_percentage_points = mean(eDNA - iDNA, na.rm = TRUE),
    .groups = "drop"
  )

local_summary <- local_pool_summary %>%
  left_join(local_pair_summary, by = c("scenario", "depth_cm")) %>%
  mutate(
    eDNA_to_iDNA_mean_ratio = eDNA_mean_percent / iDNA_mean_percent,
    control_rule_percentage_denominator = "original total prokaryotic reads"
  ) %>%
  arrange(scenario, factor(depth_cm, levels = c("20-40", "40-60")))

write_csv(
  local_detail,
  file.path(sensitivity_out_dir, "round48_control_sensitivity_NA_subsoil_by_profile.csv")
)
write_csv(
  local_summary,
  file.path(sensitivity_out_dir, "round48_control_sensitivity_NA_subsoil_summary.csv")
)

smallest_library_audit <- tibble(
  sample_id = smallest_sample,
  site_internal = smallest_meta$site,
  site_display = display_site(smallest_meta$site),
  pit = smallest_meta$pit,
  dna_pool = smallest_meta$dna_type,
  depth_cm = gsub("_", "-", smallest_meta$depths_cm, fixed = TRUE),
  original_archaeal_reads = fig1_reads$archaeal_reads[smallest_index],
  original_total_prokaryotic_reads = fig1_reads$total_reads[smallest_index],
  original_archaeal_percent = fig1_reads$arcpercent[smallest_index],
  control_filtered_archaeal_reads = sum(
    archaea_raw[!archaea_raw$asv_id %in% control_asv_ids, smallest_sample]
  ),
  affected_pair_id = paste(
    smallest_meta$site, smallest_meta$pit, smallest_meta$depths_cm, sep = "_"
  ),
  present_in_rarefied_202_table = smallest_sample %in% colnames(rare_counts),
  member_of_complete_pair_before_drop = paste(
    smallest_meta$site, smallest_meta$pit, smallest_meta$depths_cm, sep = "_"
  ) %in% canonical_concordance$pair_map$pair_id,
  paired_concordance_pairs_after_drop = nrow(
    concordance_sensitivity$drop_226_library$pair_map
  ),
  NA_subsoil_read_share_affected = smallest_meta$depths_cm %in% c("20_40", "40_60")
)
write_csv(
  smallest_library_audit,
  file.path(sensitivity_out_dir, "round48_smallest_226_library_audit.csv")
)

# Prespecified survival calls.
get_overall <- function(scenario_name) {
  concordance_sensitivity_overall %>% filter(scenario == scenario_name)
}
get_by_site <- function(scenario_name) {
  concordance_sensitivity_by_site %>% filter(scenario == scenario_name)
}
concordance_survives <- function(scenario_name) {
  overall <- get_overall(scenario_name)
  sites <- get_by_site(scenario_name)
  all(overall$statistic >= 0.70) &&
    all(overall$p_value <= 0.05) &&
    all(sites$protest_r >= 0.70)
}
site_survives <- function(scenario_name) {
  site_sensitivity$p_value[site_sensitivity$scenario == scenario_name] < 0.05
}
local_survives <- function(scenario_name) {
  values <- local_summary %>% filter(scenario == scenario_name)
  nrow(values) == 2L && all(values$eDNA_mean_percent > values$iDNA_mean_percent)
}

format_site <- function(scenario_name) {
  x <- site_sensitivity %>% filter(scenario == scenario_name)
  sprintf("R2=%.4f; F=%.3f; p=%.4g", x$R2, x$F, x$p_value)
}
format_concordance <- function(scenario_name) {
  x <- get_overall(scenario_name)
  pr <- x %>% filter(method == "PROTEST Procrustes")
  mt <- x %>% filter(method == "Mantel Spearman")
  sprintf(
    "PROTEST r=%.4f, p=%.4g; Mantel r=%.4f, p=%.4g; n=%d",
    pr$statistic, pr$p_value, mt$statistic, mt$p_value, pr$n_complete_pairs
  )
}
format_local <- function(scenario_name) {
  x <- local_summary %>%
    filter(scenario == scenario_name) %>%
    arrange(factor(depth_cm, levels = c("20-40", "40-60")))
  paste(sprintf(
    "%s cm eDNA %.2f%% vs iDNA %.2f%% (%d/%d paired directions)",
    x$depth_cm, x$eDNA_mean_percent, x$iDNA_mean_percent,
    x$eDNA_higher_pairs_n, x$n_complete_pairs
  ), collapse = "; ")
}

claim_survival <- tibble(
  claim_id = c(
    "A_ASV_profile_site_structure",
    "B_paired_cross_pool_concordance",
    "C_NA_subsoil_archaeal_read_share"
  ),
  survival_rule = c(
    "profile-level ASV robust-Aitchison site p < 0.05",
    paste(
      "overall PROTEST and Mantel p <= 0.05 with r >= 0.70, and",
      "within-site PROTEST r >= 0.70 at all four sites"
    ),
    "mean eDNA archaeal percentage remains greater than iDNA at both depths"
  ),
  original = c(
    format_site("primary_unfiltered"),
    format_concordance("primary_unfiltered"),
    format_local("primary_unfiltered")
  ),
  control_filtered = c(
    format_site("control_filtered"),
    format_concordance("control_filtered"),
    format_local("control_filtered")
  ),
  control_filtered_status = c(
    ifelse(site_survives("control_filtered"), "SURVIVED", "DID_NOT_SURVIVE"),
    ifelse(concordance_survives("control_filtered"), "SURVIVED", "DID_NOT_SURVIVE"),
    ifelse(local_survives("control_filtered"), "SURVIVED", "DID_NOT_SURVIVE")
  ),
  drop_226_library = c(
    format_site("drop_226_library"),
    format_concordance("drop_226_library"),
    format_local("drop_226_library")
  ),
  drop_226_status = c(
    ifelse(site_survives("drop_226_library"), "SURVIVED", "DID_NOT_SURVIVE"),
    ifelse(concordance_survives("drop_226_library"), "SURVIVED", "DID_NOT_SURVIVE"),
    ifelse(local_survives("drop_226_library"), "SURVIVED", "DID_NOT_SURVIVE")
  ),
  combined_control_plus_drop = c(
    format_site("control_filtered_plus_drop_226"),
    format_concordance("control_filtered_plus_drop_226"),
    format_local("control_filtered_plus_drop_226")
  )
)
write_csv(
  claim_survival,
  file.path(sensitivity_out_dir, "round48_core_claim_survival.csv")
)

input_provenance <- tibble(
  input = names(input_paths),
  project_relative_path = sub(
    paste0("^", project_root, "/"), "", unname(input_paths)
  ),
  bytes = as.numeric(file.info(input_paths)$size),
  md5 = unname(tools::md5sum(input_paths))
)
write_csv(
  input_provenance,
  file.path(sensitivity_out_dir, "round48_input_provenance.csv")
)
capture.output(
  sessionInfo(),
  file = file.path(sensitivity_out_dir, "round48_R_session_info.txt")
)

# Hard validation against the reported canonical values.
canonical_pr <- canonical_concordance$overall %>%
  filter(method == "PROTEST Procrustes")
canonical_mt <- canonical_concordance$overall %>%
  filter(method == "Mantel Spearman")
expected_site_r <- c(AZ = 0.9394950397622431, SG = 0.9945611394726042,
                     LC = 0.9251237635137265, NB = 0.9700384622455281)
if (abs(canonical_pr$statistic - 0.984639884367172) > 1e-12 ||
    abs(canonical_mt$statistic - 0.7803787720920279) > 1e-12 ||
    canonical_pr$n_complete_pairs != 39L ||
    any(abs(
      setNames(canonical_concordance$by_site$protest_r,
               canonical_concordance$by_site$site_internal)[names(expected_site_r)] -
        expected_site_r
    ) > 1e-12)) {
  stop("Canonical concordance no longer reproduces the independently verified values.")
}
primary_site <- site_sensitivity %>% filter(scenario == "primary_unfiltered")
if (abs(primary_site$R2 - 0.738023409902886) > 1e-12 ||
    abs(primary_site$F - 7.51235988558422) > 1e-10 ||
    primary_site$p_value != 0.0014) {
  stop("The original robust-Aitchison site result was not reproduced exactly.")
}

cat(sprintf(
  "Confirmed concordance: PROTEST r=%.4f, p=%.4g, n=%d; Mantel r=%.4f, p=%.4g.\n",
  canonical_pr$statistic, canonical_pr$p_value,
  canonical_pr$n_complete_pairs, canonical_mt$statistic, canonical_mt$p_value
))
cat(sprintf(
  "Within-site PROTEST: %s.\n",
  paste(sprintf(
    "%s %.3f",
    canonical_concordance$by_site$site_display,
    canonical_concordance$by_site$protest_r
  ), collapse = " / ")
))
walk(seq_len(nrow(claim_survival)), function(i) {
  cat(sprintf(
    "%s: %s under control filter (%s).\n",
    claim_survival$claim_id[i],
    claim_survival$control_filtered_status[i],
    claim_survival$control_filtered[i]
  ))
})
cat(sprintf(
  "Dropping the 226-read library (%s): site %s; concordance %s; NA subsoil %s.\n",
  smallest_sample,
  claim_survival$drop_226_status[1],
  claim_survival$drop_226_status[2],
  claim_survival$drop_226_status[3]
))
