#!/usr/bin/env Rscript

# Round55 profile-aware analysis of paired archaeal read-share offsets.
#
# The inferential unit is the soil profile (three profiles per site), not an
# individual sequencing read or horizon library. Horizon-pair offsets are
# retained for transparent descriptive plotting. The read-share response is a
# relative, primer-dependent library percentage and is not an absolute-abundance
# or DNA-mass measurement.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

SITE_LEVELS <- c("AZ", "SG", "LC", "NA")
DEPTH_LEVELS <- c("0-5", "5-10", "10-20", "20-40", "40-60")
SITE_COLOURS <- c(
  AZ = "#C25B3F",
  SG = "#D6A141",
  LC = "#6FA055",
  "NA" = "#3D6F97"
)
PROFILE_SHAPES <- c(`1` = 21, `2` = 22, `3` = 24)
PROFILE_LINETYPES <- c(`1` = "solid", `2` = "22", `3` = "42")
BOOTSTRAP_SEED <- 20260716L
BOOTSTRAP_REPLICATES <- 100000L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  project_root <- dirname(dirname(script_path))
} else {
  project_root <- normalizePath(".", mustWork = TRUE)
}

input_paths <- c(
  pair_map = file.path(
    project_root, "analysis", "paired_pool_similarity",
    "round48_complete_pair_map.csv"
  ),
  read_share = file.path(
    project_root, "analysis", "reframed_figures",
    "Fig1_archaeal_read_fraction_source.csv"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

analysis_dir <- file.path(project_root, "analysis", "round55")
figure_path <- file.path(
  project_root, "figures", "Figure_S7_paired_read_share_offsets_round55.pdf"
)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)

pair_map <- read_csv(
  input_paths[["pair_map"]],
  na = character(),
  show_col_types = FALSE
)
read_share <- read_csv(
  input_paths[["read_share"]],
  na = character(),
  show_col_types = FALSE
)

required_pair_columns <- c(
  "pair_id", "site_display", "pit", "depth_order", "depth_cm",
  "iDNA_sample", "eDNA_sample"
)
required_share_columns <- c("sample_id", "arcpercent")
if (!all(required_pair_columns %in% names(pair_map))) {
  stop("Pair map lacks one or more required columns.")
}
if (!all(required_share_columns %in% names(read_share))) {
  stop("Read-share source lacks one or more required columns.")
}
if (anyDuplicated(pair_map$pair_id) || anyDuplicated(read_share$sample_id)) {
  stop("Pair IDs and read-share sample IDs must be unique.")
}

i_share <- read_share %>%
  transmute(iDNA_sample = sample_id, iDNA_pct = arcpercent)
e_share <- read_share %>%
  transmute(eDNA_sample = sample_id, eDNA_pct = arcpercent)

paired <- pair_map %>%
  select(all_of(required_pair_columns)) %>%
  left_join(i_share, by = "iDNA_sample") %>%
  left_join(e_share, by = "eDNA_sample") %>%
  transmute(
    pair_id,
    site = factor(site_display, levels = SITE_LEVELS),
    pit = as.integer(pit),
    profile = paste(site_display, pit, sep = "_P"),
    depth_order = as.integer(depth_order),
    depth_cm = factor(depth_cm, levels = DEPTH_LEVELS),
    iDNA_pct,
    eDNA_pct,
    diff_pp = eDNA_pct - iDNA_pct
  ) %>%
  arrange(site, pit, depth_order)

if (nrow(paired) != 39L || anyNA(paired)) {
  stop("Expected 39 fully joined horizon pairs with no missing values.")
}
if (!identical(sort(unique(as.character(paired$site))), sort(SITE_LEVELS))) {
  stop("Unexpected site inventory after joining source data.")
}

profile_summary <- paired %>%
  group_by(site, pit, profile) %>%
  summarise(
    n_horizon_pairs = n(),
    profile_mean_pp = mean(diff_pp),
    .groups = "drop"
  ) %>%
  arrange(site, pit)

if (nrow(profile_summary) != 12L || any(profile_summary$n_horizon_pairs < 1L)) {
  stop("Expected 12 non-empty profile means.")
}

paired <- paired %>%
  left_join(
    profile_summary %>% select(site, pit, profile_mean_pp),
    by = c("site", "pit")
  )

site_profile_means <- profile_summary %>%
  group_by(site) %>%
  summarise(
    n_profiles = n(),
    site_equal_profile_mean_pp = mean(profile_mean_pp),
    .groups = "drop"
  )

pooled_depth <- paired %>%
  group_by(depth_order, depth_cm) %>%
  summarise(
    n_horizon_pairs = n(),
    pooled_pair_mean_pp = mean(diff_pp),
    .groups = "drop"
  ) %>%
  arrange(depth_order)

site_depth <- paired %>%
  group_by(site, depth_order, depth_cm) %>%
  summarise(
    n_horizon_pairs = n(),
    site_depth_pair_mean_pp = mean(diff_pp),
    .groups = "drop"
  ) %>%
  arrange(site, depth_order)

# Reconciliation gates from the two independently reviewed source-data audits.
stopifnot(
  sum(paired$diff_pp > 0) == 37L,
  sum(paired$diff_pp < 0) == 2L,
  all(profile_summary$profile_mean_pp > 0),
  all(site_profile_means$n_profiles == 3L)
)
expected_site_means <- c(AZ = 1.18, SG = 2.31, LC = 2.24, "NA" = 6.65)
observed_site_means <- setNames(
  site_profile_means$site_equal_profile_mean_pp,
  as.character(site_profile_means$site)
)
if (any(abs(observed_site_means[names(expected_site_means)] - expected_site_means) > 0.01)) {
  stop("Equal-profile site means do not reconcile with the reviewed values.")
}
expected_pooled_depth <- c(0.749, 1.133, 2.089, 4.540, 6.542)
if (any(abs(pooled_depth$pooled_pair_mean_pp - expected_pooled_depth) > 0.001)) {
  stop("Pooled depth means do not reconcile with the reviewed values.")
}
strictly_increasing_sites <- site_depth %>%
  group_by(site) %>%
  summarise(strictly_increasing = all(diff(site_depth_pair_mean_pp) > 0), .groups = "drop")
if (!identical(
  as.character(strictly_increasing_sites$site[strictly_increasing_sites$strictly_increasing]),
  "NA"
)) {
  stop("Site-specific depth heterogeneity no longer matches the reviewed data.")
}

# Exact one-sided sign test on the 12 independent profile-mean directions.
n_profiles <- nrow(profile_summary)
n_positive_profiles <- sum(profile_summary$profile_mean_pp > 0)
sign_test_result <- binom.test(
  n_positive_profiles,
  n_profiles,
  p = 0.5,
  alternative = "greater"
)

# Stratified profile bootstrap: resample the three profiles within each site,
# then average the 12 resampled profile means with equal profile weight.
set.seed(BOOTSTRAP_SEED)
site_profile_vectors <- split(profile_summary$profile_mean_pp, profile_summary$site)
if (!all(lengths(site_profile_vectors) == 3L)) {
  stop("The stratified bootstrap requires three profiles per site.")
}
bootstrap_means <- vapply(seq_len(BOOTSTRAP_REPLICATES), function(iteration) {
  resampled <- unlist(lapply(
    site_profile_vectors,
    function(x) sample(x, size = length(x), replace = TRUE)
  ), use.names = FALSE)
  mean(resampled)
}, numeric(1))
bootstrap_ci <- unname(quantile(
  bootstrap_means,
  probs = c(0.025, 0.975),
  names = FALSE,
  type = 7
))

profile_effect <- mean(profile_summary$profile_mean_pp)
profile_median <- median(profile_summary$profile_mean_pp)
inference <- tibble(
  independent_unit = "profile",
  n_profiles = n_profiles,
  positive_profile_means = n_positive_profiles,
  negative_profile_means = sum(profile_summary$profile_mean_pp < 0),
  zero_profile_means = sum(profile_summary$profile_mean_pp == 0),
  equal_profile_mean_pp = profile_effect,
  profile_median_pp = profile_median,
  profile_min_pp = min(profile_summary$profile_mean_pp),
  profile_max_pp = max(profile_summary$profile_mean_pp),
  bootstrap_ci_level = 0.95,
  bootstrap_ci_low_pp = bootstrap_ci[1],
  bootstrap_ci_high_pp = bootstrap_ci[2],
  bootstrap_method = "site-stratified profile resampling; equal profile weight; percentile CI",
  bootstrap_replicates = BOOTSTRAP_REPLICATES,
  bootstrap_seed = BOOTSTRAP_SEED,
  sign_test = "exact one-sided binomial sign test against Pr(positive)=0.5",
  sign_test_p_value = unname(sign_test_result$p.value),
  horizon_pairs_total = nrow(paired),
  horizon_pairs_positive = sum(paired$diff_pp > 0)
)

write_csv(
  paired %>%
    transmute(
      pair_id,
      site = as.character(site),
      pit,
      depth_cm = as.character(depth_cm),
      iDNA_pct,
      eDNA_pct,
      diff_pp,
      profile_mean_pp
    ),
  file.path(analysis_dir, "paired_read_share_offsets.csv")
)
write_csv(profile_summary, file.path(analysis_dir, "profile_read_share_offsets.csv"))
write_csv(inference, file.path(analysis_dir, "profile_read_share_inference.csv"))
write_csv(site_profile_means, file.path(analysis_dir, "site_equal_profile_mean_offsets.csv"))
write_csv(pooled_depth, file.path(analysis_dir, "pooled_depth_pair_mean_offsets.csv"))
write_csv(site_depth, file.path(analysis_dir, "site_depth_pair_mean_offsets.csv"))
write_csv(
  strictly_increasing_sites,
  file.path(analysis_dir, "site_depth_monotonicity_audit.csv")
)

# Figure S7: every horizon-pair offset, joined within its profile, plus the
# descriptive pooled depth mean. The profile panels make site heterogeneity
# visible; the pooled panel is explicitly labelled as descriptive.
plot_pairs <- paired %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    pit_factor = factor(pit, levels = 1:3),
    depth_cm = factor(depth_cm, levels = DEPTH_LEVELS)
  )

base_theme <- theme_classic(base_size = 10.5, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 11.2, hjust = 0),
    plot.subtitle = element_text(size = 9.2, colour = "grey30", hjust = 0),
    axis.title = element_text(size = 9.8),
    axis.text = element_text(size = 8.6, colour = "grey20"),
    axis.line = element_line(linewidth = 0.35, colour = "grey25"),
    axis.ticks = element_line(linewidth = 0.3, colour = "grey25"),
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text = element_text(face = "bold", size = 9.3),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 8.8),
    legend.text = element_text(size = 8.5),
    panel.spacing.x = grid::unit(0.65, "lines"),
    plot.margin = margin(4, 6, 3, 4)
  )

p_profiles <- ggplot(
  plot_pairs,
  aes(x = depth_cm, y = diff_pp, group = profile)
) +
  geom_hline(yintercept = 0, linewidth = 0.36, linetype = "dashed", colour = "grey55") +
  geom_line(
    aes(colour = site, linetype = pit_factor),
    linewidth = 0.70,
    alpha = 0.86,
    na.rm = TRUE
  ) +
  geom_point(
    aes(fill = site, shape = pit_factor),
    size = 2.65,
    stroke = 0.55,
    colour = "white"
  ) +
  facet_grid(. ~ site, labeller = as_labeller(c(
    AZ = "AZ · Pan de Azúcar",
    SG = "SG · Santa Gracia",
    LC = "LC · La Campana",
    "NA" = "NA · Nahuelbuta"
  ))) +
  scale_colour_manual(values = SITE_COLOURS, guide = "none") +
  scale_fill_manual(values = SITE_COLOURS, guide = "none") +
  scale_shape_manual(
    values = PROFILE_SHAPES,
    name = "Profile",
    labels = paste("Profile", 1:3)
  ) +
  scale_linetype_manual(
    values = PROFILE_LINETYPES,
    name = "Profile",
    labels = paste("Profile", 1:3)
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Profile trajectories show site-specific depth heterogeneity",
    subtitle = "39 complete horizon pairs across 12 profiles; positive values indicate eDNA > iDNA",
    x = NULL,
    y = "eDNA − iDNA archaeal read share\n(percentage points)"
  ) +
  base_theme

p_profiles <- p_profiles +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.8))

p_pooled <- ggplot(
  pooled_depth,
  aes(x = factor(depth_cm, levels = DEPTH_LEVELS), y = pooled_pair_mean_pp, group = 1)
) +
  geom_hline(yintercept = 0, linewidth = 0.36, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.75, colour = "grey25") +
  geom_point(shape = 21, size = 2.9, stroke = 0.45, fill = "grey25", colour = "white") +
  geom_text(
    aes(y = pooled_pair_mean_pp + 0.48, label = paste0("n = ", n_horizon_pairs)),
    vjust = 0.5,
    size = 2.75,
    family = "Helvetica",
    colour = "grey25"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.18))) +
  labs(
    title = "Descriptive pooled paired mean increases across depth classes",
    subtitle = "The pooled sequence does not imply monotonic change within every site",
    x = "Soil depth (cm)",
    y = "Pooled mean offset\n(percentage points)"
  ) +
  base_theme +
  theme(legend.position = "none")

figure <- (p_profiles / p_pooled) +
  plot_layout(heights = c(2.25, 1)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste0(
      "Read share is a relative, universal-primer library outcome; it is not absolute abundance or DNA mass.\n",
      "Profile-level inference: 12/12 means positive; exact one-sided sign-test p = 1/4096 (",
      formatC(sign_test_result$p.value, format = "f", digits = 6), ")."
    ),
    theme = theme(
      plot.tag = element_text(
        family = "Helvetica", face = "bold", size = 15, colour = "black"
      ),
      plot.caption = element_text(
        family = "Helvetica", size = 8.2, colour = "grey30", hjust = 0,
        margin = margin(t = 5)
      )
    )
  )

ggsave(
  filename = figure_path,
  plot = figure,
  width = 190,
  height = 150,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white"
)

session_path <- file.path(analysis_dir, "R_session_info.txt")
sink(session_path)
cat("Round55 profile-aware read-share analysis\n")
cat("Script: script/104_round55_profile_read_share_offsets.R\n")
cat("Bootstrap seed:", BOOTSTRAP_SEED, "\n")
cat("Bootstrap replicates:", BOOTSTRAP_REPLICATES, "\n\n")
print(sessionInfo())
sink()

cat(sprintf(
  paste0(
    "PASS Round55 profile analysis: %d/%d horizon pairs positive; ",
    "%d/%d profile means positive; equal-profile mean = %.6f pp; ",
    "95%% stratified-bootstrap CI = [%.6f, %.6f] pp; exact one-sided sign p = %.12f.\n"
  ),
  sum(paired$diff_pp > 0), nrow(paired),
  n_positive_profiles, n_profiles,
  profile_effect, bootstrap_ci[1], bootstrap_ci[2], sign_test_result$p.value
))
cat("Wrote:", file.path(analysis_dir, "paired_read_share_offsets.csv"), "\n")
cat("Wrote:", figure_path, "\n")
