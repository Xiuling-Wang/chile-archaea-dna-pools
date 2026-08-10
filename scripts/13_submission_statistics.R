suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

script_path <- commandArgs(FALSE) |>
  grep(pattern = "^--file=", value = TRUE) |>
  sub(pattern = "^--file=", replacement = "")
if (length(script_path) == 0) {
  script_path <- file.path("script", "13_submission_statistics.R")
}

project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
data_dir <- file.path(project_root, "data")
report_dir <- file.path(project_root, "analysis", "submission_statistics")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

asv_path <- file.path(data_dir, "arc_60cm_202rare/asv_202.txt")
tax_path <- file.path(data_dir, "arc_100_202rarefild/tax_202.txt")
env_path <- file.path(data_dir, "env_2024.csv")

asv_raw <- read.table(asv_path, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
tax <- read.table(tax_path, header = TRUE, check.names = FALSE, sep = "\t")
env <- read.csv(env_path, header = TRUE, check.names = FALSE)

env_60 <- env %>%
  filter(depth < 60, sample_id %in% colnames(asv_raw)) %>%
  mutate(
    site = factor(site, levels = c("AZ", "SG", "LC", "NB")),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA")),
    depth_order = as.numeric(depth_order),
    depth_mid = as.numeric(depth)
  ) %>%
  arrange(match(sample_id, colnames(asv_raw)))

asv <- asv_raw[, env_60$sample_id, drop = FALSE]
asv <- asv[rowSums(asv) > 0, colSums(asv) > 0, drop = FALSE]
env_60 <- env_60 %>% filter(sample_id %in% colnames(asv))
asv <- asv[, env_60$sample_id, drop = FALSE]
stopifnot(identical(colnames(asv), env_60$sample_id))

asv_t <- t(asv) %>% as.data.frame()
asv_hell <- decostand(asv_t, method = "hellinger")
permutation_n <- 9999
permutation_seed <- 20260522

sample_summary <- env_60 %>%
  count(site, dna_type, name = "n_samples") %>%
  arrange(site, dna_type)
write.csv(sample_summary, file.path(report_dir, "sample_summary.csv"), row.names = FALSE)

alpha <- tibble(
  sample_id = rownames(asv_t),
  Shannon = diversity(asv_t, index = "shannon"),
  Simpson = diversity(asv_t, index = "simpson"),
  Chao1 = estimateR(asv_t)["S.chao1", ] %>% as.numeric(),
  Observed_ASVs = specnumber(asv_t)
) %>%
  left_join(env_60 %>% select(sample_id, site, dna_type, depth_order, depth), by = "sample_id")
write.csv(alpha, file.path(report_dir, "alpha_diversity_60cm.csv"), row.names = FALSE)

alpha_summary <- alpha %>%
  group_by(site, dna_type) %>%
  summarise(
    n = n(),
    Shannon_mean = mean(Shannon),
    Shannon_sd = sd(Shannon),
    Chao1_mean = mean(Chao1),
    Chao1_sd = sd(Chao1),
    Simpson_mean = mean(Simpson),
    Simpson_sd = sd(Simpson),
    .groups = "drop"
  )
write.csv(alpha_summary, file.path(report_dir, "alpha_summary_by_site_dna.csv"), row.names = FALSE)

set.seed(permutation_seed)
permanova_all <- adonis2(
  asv_hell ~ site + depth_mid + dna_type,
  data = env_60,
  permutations = permutation_n,
  method = "bray",
  by = "margin",
  na.action = na.fail
)
capture.output(permanova_all, file = file.path(report_dir, "PERMANOVA_all_site_depth_dna.txt"))

run_pool_permanova <- function(pool) {
  env_pool <- env_60 %>% filter(dna_type == pool)
  asv_pool <- asv[, env_pool$sample_id, drop = FALSE]
  asv_pool <- asv_pool[rowSums(asv_pool) > 0, , drop = FALSE]
  asv_pool_hell <- decostand(t(asv_pool), method = "hellinger")
  set.seed(permutation_seed)
  adonis2(
    asv_pool_hell ~ site + depth_mid,
    data = env_pool,
    permutations = permutation_n,
    method = "bray",
    by = "margin",
    na.action = na.fail
  )
}
permanova_i <- run_pool_permanova("iDNA")
permanova_e <- run_pool_permanova("eDNA")
capture.output(permanova_i, file = file.path(report_dir, "PERMANOVA_iDNA_site_depth.txt"))
capture.output(permanova_e, file = file.path(report_dir, "PERMANOVA_eDNA_site_depth.txt"))

env_vars <- c("pH", "Conductivity", "moisture", "CN", "Feo", "Alo", "Mno", "Sio", "NH4", "NO3", "Po", "Pi")
env_db <- env_60 %>% select(all_of(env_vars), site, dna_type, depth_order, depth_mid, sample_id) %>% na.omit()
asv_db <- asv[, env_db$sample_id, drop = FALSE]
asv_db <- asv_db[rowSums(asv_db) > 0, , drop = FALSE]
asv_db_hell <- decostand(t(asv_db), method = "hellinger")

db_formula <- as.formula(paste("asv_db_hell ~", paste(env_vars, collapse = " + ")))
db_all <- dbrda(db_formula, data = env_db, distance = "bray", scale = TRUE)
db_all_r2 <- RsquareAdj(db_all)
capture.output(db_all_r2, file = file.path(report_dir, "dbRDA_all_R2.txt"))

run_pool_dbrda <- function(pool) {
  env_pool <- env_db %>% filter(dna_type == pool)
  asv_pool <- asv[, env_pool$sample_id, drop = FALSE]
  asv_pool <- asv_pool[rowSums(asv_pool) > 0, , drop = FALSE]
  asv_pool_hell <- decostand(t(asv_pool), method = "hellinger")
  pool_formula <- as.formula(paste("asv_pool_hell ~", paste(env_vars, collapse = " + ")))
  mod <- dbrda(pool_formula, data = env_pool, distance = "bray", scale = TRUE)
  list(model = mod, r2 = RsquareAdj(mod))
}
db_i <- run_pool_dbrda("iDNA")
db_e <- run_pool_dbrda("eDNA")
capture.output(db_i$r2, file = file.path(report_dir, "dbRDA_iDNA_R2.txt"))
capture.output(db_e$r2, file = file.path(report_dir, "dbRDA_eDNA_R2.txt"))

tax_clean <- tax %>%
  mutate(asv_id = as.character(asv_id)) %>%
  filter(asv_id %in% rownames(asv))

rel_long <- asv %>%
  as.data.frame() %>%
  rownames_to_column("asv_id") %>%
  pivot_longer(-asv_id, names_to = "sample_id", values_to = "count") %>%
  left_join(tax_clean, by = "asv_id") %>%
  left_join(env_60 %>% select(sample_id, site, dna_type, depth_order, depth), by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(rel_abundance = count / sum(count)) %>%
  ungroup()

phylum_overall <- rel_long %>%
  group_by(Phylum) %>%
  summarise(relative_abundance_percent = 100 * sum(count) / sum(rel_long$count), .groups = "drop") %>%
  arrange(desc(relative_abundance_percent))
write.csv(phylum_overall, file.path(report_dir, "phylum_overall_relative_abundance.csv"), row.names = FALSE)

phylum_by_site_dna <- rel_long %>%
  group_by(site, dna_type, Phylum) %>%
  summarise(count = sum(count), .groups = "drop_last") %>%
  mutate(relative_abundance_percent = 100 * count / sum(count)) %>%
  ungroup() %>%
  arrange(site, dna_type, desc(relative_abundance_percent))
write.csv(phylum_by_site_dna, file.path(report_dir, "phylum_relative_abundance_by_site_dna.csv"), row.names = FALSE)

class_overall <- rel_long %>%
  group_by(Class) %>%
  summarise(relative_abundance_percent = 100 * sum(count) / sum(rel_long$count), .groups = "drop") %>%
  arrange(desc(relative_abundance_percent))
write.csv(class_overall, file.path(report_dir, "class_overall_relative_abundance.csv"), row.names = FALSE)

format_permanova <- function(x) {
  as.data.frame(x) %>%
    rownames_to_column("term") %>%
    filter(!term %in% c("Residual", "Total")) %>%
    transmute(
      term,
      Df,
      SumOfSqs,
      R2,
      F = F,
      p = `Pr(>F)`
    )
}

summary_lines <- c(
  "# Archaea submission statistics",
  "",
  sprintf("Input ASV table: %s", asv_path),
  sprintf("Retained samples: %d; retained ASVs: %d; rarefaction depth: %d reads per sample.", nrow(env_60), nrow(asv), unique(colSums(asv))[1]),
  sprintf("PERMANOVA settings: Bray-Curtis dissimilarity after Hellinger transformation; depth treated as interval midpoint; marginal tests; %d permutations; set.seed(%d).", permutation_n, permutation_seed),
  "",
  "## PERMANOVA: all samples",
  capture.output(print(format_permanova(permanova_all))),
  "",
  "## PERMANOVA: iDNA",
  capture.output(print(format_permanova(permanova_i))),
  "",
  "## PERMANOVA: eDNA",
  capture.output(print(format_permanova(permanova_e))),
  "",
  "## dbRDA R2",
  sprintf("All samples: R2 = %.3f; adjusted R2 = %.3f.", db_all_r2$r.squared, db_all_r2$adj.r.squared),
  sprintf("iDNA: R2 = %.3f; adjusted R2 = %.3f.", db_i$r2$r.squared, db_i$r2$adj.r.squared),
  sprintf("eDNA: R2 = %.3f; adjusted R2 = %.3f.", db_e$r2$r.squared, db_e$r2$adj.r.squared),
  "",
  "## Dominant phyla",
  capture.output(print(phylum_overall, n = Inf)),
  "",
  "## Alpha diversity summary",
  capture.output(print(alpha_summary, n = Inf))
)
writeLines(summary_lines, file.path(report_dir, "submission_statistics_summary.md"))

message("Wrote submission statistics to: ", report_dir)
