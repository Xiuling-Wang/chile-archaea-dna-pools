#!/usr/bin/env Rscript

# Reproducibility audit for the only exchangeability restriction supported by
# the retained 202-read table: iDNA/eDNA labels may be swapped within each
# complete site x profile x depth pair. Site-level profile swaps are not valid
# because retained profile sizes are unequal.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(vegan)
  library(permute)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
asv_path <- file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt")
env_path <- file.path(project_root, "data", "env_2024.csv")
out_dir <- file.path(project_root, "analysis", "reproducibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(asv_path) || !file.exists(env_path)) {
  stop("Required input missing: ", paste(c(asv_path, env_path)[!file.exists(c(asv_path, env_path))], collapse = "; "))
}

asv_raw <- read.delim(asv_path, header = TRUE, row.names = 1, check.names = FALSE)
env <- read.csv(env_path, check.names = FALSE) %>%
  filter(depth < 60, sample_id %in% colnames(asv_raw), dna_type %in% c("iDNA", "eDNA")) %>%
  mutate(
    profile_id = interaction(site, pit, drop = TRUE),
    pair_id = interaction(site, pit, depth_order, drop = TRUE),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA"))
  ) %>%
  arrange(match(sample_id, colnames(asv_raw)))

asv <- asv_raw[, env$sample_id, drop = FALSE]
asv <- asv[rowSums(asv) > 0, colSums(asv) > 0, drop = FALSE]
env <- env %>% filter(sample_id %in% colnames(asv)) %>% arrange(match(sample_id, colnames(asv)))
asv <- asv[, env$sample_id, drop = FALSE]
stopifnot(identical(colnames(asv), env$sample_id))

profile_sizes <- env %>%
  count(profile_id, site, pit, name = "n_samples") %>%
  arrange(site, pit)
write_csv(profile_sizes, file.path(out_dir, "round25_profile_sizes.csv"))

pair_counts <- env %>%
  count(pair_id, site, pit, depth_order, dna_type, name = "n") %>%
  pivot_wider(names_from = dna_type, values_from = n, values_fill = 0) %>%
  mutate(complete_pair = iDNA == 1 & eDNA == 1)
write_csv(pair_counts, file.path(out_dir, "round25_pair_retention.csv"))

complete_pair_ids <- pair_counts %>% filter(complete_pair) %>% pull(pair_id)
paired_env <- env %>%
  filter(pair_id %in% complete_pair_ids) %>%
  mutate(pair_id = droplevels(factor(as.character(pair_id)))) %>%
  arrange(pair_id, dna_type)

if (nrow(paired_env) != 2L * length(complete_pair_ids)) {
  stop("Complete pair retention is inconsistent with two samples per pair.")
}
if (!all(table(paired_env$pair_id) == 2L)) {
  stop("Each restricted permutation block must contain exactly one iDNA and one eDNA sample.")
}

paired_asv <- asv[, paired_env$sample_id, drop = FALSE]
paired_asv <- paired_asv[rowSums(paired_asv) > 0, , drop = FALSE]
paired_hellinger <- decostand(t(paired_asv), method = "hellinger")

# With two observations per block, free within-block permutation is a label swap
# or identity operation. It is valid for the paired DNA-pool contrast only.
restricted_pool_control <- how(
  within = Within(type = "free"),
  blocks = paired_env$pair_id,
  nperm = 9999
)

set.seed(20260710)
# Pair identity is included as a blocking term so the residual denominator also
# represents the paired design (review MC4). Without it, community ~ dna_type
# gives an uninterpretable F < 1; the R2 and restricted p value are unchanged.
paired_pool_permanova <- adonis2(
  paired_hellinger ~ pair_id + dna_type,
  data = paired_env,
  method = "bray",
  permutations = restricted_pool_control,
  by = "margin"
)

pool_result <- as.data.frame(paired_pool_permanova) %>%
  tibble::rownames_to_column("term") %>%
  filter(term == "dna_type") %>%
  transmute(
    test = "paired iDNA-eDNA pool contrast",
    n_complete_pairs = length(complete_pair_ids),
    n_samples = nrow(paired_env),
    Df = Df,
    SumOfSqs = SumOfSqs,
    transformation = "Hellinger",
    dissimilarity = "Bray-Curtis",
    permutation_scheme = "free iDNA/eDNA label exchange within site x pit x depth pair",
    permutations = 9999L,
    F = F,
    R2 = R2,
    p_value = `Pr(>F)`
  )
write_csv(pool_result, file.path(out_dir, "round25_pair_restricted_pool_permanova.csv"))
capture.output(
  paired_pool_permanova,
  file = file.path(out_dir, "round25_pair_restricted_pool_permanova.txt")
)

readr::write_lines(
  c(
    "# Round25 restricted-permutation design audit",
    "",
    "## Scope",
    "",
    sprintf("The 202-read archaeal table retained %d samples and %d ASVs.", nrow(env), nrow(asv)),
    sprintf("There were %d complete iDNA/eDNA pairs (n = %d samples) defined by site x pit x depth.", length(complete_pair_ids), nrow(paired_env)),
    "",
    "## Valid restricted test",
    "",
    "The paired DNA-pool contrast was tested by permuting iDNA/eDNA labels freely within each complete site x pit x depth pair (9999 permutations; Hellinger transformation and Bray-Curtis dissimilarity). Pair identity was included as a blocking term (community ~ pair_id + dna_type, marginal test) so the residual denominator represents the paired design.",
    sprintf("Result: F = %.4f, R2 = %.5f, p = %.4g.", pool_result$F, pool_result$R2, pool_result$p_value),
    "",
    "## Why no profile-restricted site or dbRDA p-values are reported",
    "",
    "The retained profile blocks are unbalanced. Whole-profile exchange therefore violates the balanced-plot requirement in permute and cannot supply a valid restricted test of the site effect. Site and depth effect sizes in the manuscript remain the standard marginal PERMANOVA descriptions, not profile-restricted p-values. The dbRDA and hierarchical-partitioning panels are presented as exploratory associations/model contributions; restricted dbRDA significance claims are not retained.",
    "",
    "## Outputs",
    "",
    "- round25_profile_sizes.csv",
    "- round25_pair_retention.csv",
    "- round25_pair_restricted_pool_permanova.csv",
    "- round25_pair_restricted_pool_permanova.txt"
  ),
  file.path(out_dir, "round25_design_audit.md")
)

message("Wrote restricted-permutation audit to: ", out_dir)
