#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

script_path <- commandArgs(FALSE) |>
  grep(pattern = "^--file=", value = TRUE) |>
  sub(pattern = "^--file=", replacement = "")
if (length(script_path) == 0) {
  script_path <- file.path("scripts", "68_FigS3_iDNA_eDNA_ASV_overlap.R")
}

project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
house_style_run <- identical(Sys.getenv("ROUND48_HOUSESTYLE_OUTPUT"), "1")
house_style_v2_run <- house_style_run &&
  identical(Sys.getenv("ROUND48_FIGS3_HOUSESTYLE_V2"), "1")
asv_path <- file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt")
env_path <- file.path(project_root, "data", "env_2024.csv")
candidate_dir <- file.path(project_root, "figures", "candidates")
figure_path <- if (house_style_v2_run) {
  file.path(candidate_dir, "FigS3_round48_housestyle_v2.pdf")
} else if (house_style_run) {
  file.path(candidate_dir, "FigS3_round48_housestyle.pdf")
} else {
  file.path(project_root, "figures", "FigS3_iDNA_eDNA_ASV_overlap.pdf")
}
summary_path <- if (house_style_run) {
  file.path(tempdir(), "chile_archaea_round48_housestyle", "FigS3", "overlap_summary.csv")
} else {
  file.path(
    project_root,
    "analysis",
    "reframed_figures",
    "FigS3_iDNA_eDNA_ASV_overlap_summary.csv"
  )
}
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

site_levels <- c("AZ", "SG", "LC", "NA")
depth_levels <- c("0-5", "5-10", "10-20", "20-40", "40-60")
pool_colours <- c("iDNA" = "#2474B5", "eDNA" = "#E88919")

asv <- read.table(
  asv_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)
env <- read.csv(env_path, check.names = FALSE) |>
  filter(sample_id %in% colnames(asv), depth < 60) |>
  mutate(
    site = recode(site, "NB" = "NA"),
    site = factor(site, levels = site_levels),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA")),
    depth_group = case_when(
      depth <= 5 ~ "0-5",
      depth <= 10 ~ "5-10",
      depth <= 20 ~ "10-20",
      depth <= 40 ~ "20-40",
      depth <= 60 ~ "40-60",
      TRUE ~ NA_character_
    ),
    depth_group = factor(depth_group, levels = depth_levels)
  ) |>
  filter(!is.na(site), !is.na(depth_group), !is.na(dna_type))

stopifnot(nrow(env) == 95, setequal(env$sample_id, colnames(asv)))

asv_set <- function(sample_ids) {
  sample_ids <- intersect(sample_ids, colnames(asv))
  if (length(sample_ids) == 0) return(character())
  rownames(asv)[rowSums(asv[, sample_ids, drop = FALSE]) > 0]
}

overlap <- expand_grid(site = site_levels, depth_group = depth_levels) |>
  pmap_dfr(function(site, depth_group) {
    site_value <- site
    depth_value <- depth_group
    group_env <- env |>
      filter(
        as.character(.data$site) == .env$site_value,
        as.character(.data$depth_group) == .env$depth_value
      )
    i_samples <- group_env |> filter(dna_type == "iDNA") |> pull(sample_id)
    e_samples <- group_env |> filter(dna_type == "eDNA") |> pull(sample_id)
    i_set <- asv_set(i_samples)
    e_set <- asv_set(e_samples)
    shared <- length(intersect(i_set, e_set))
    tibble(
      site = factor(site, levels = site_levels),
      depth_group = factor(depth_group, levels = depth_levels),
      n_iDNA = length(i_samples),
      n_eDNA = length(e_samples),
      iDNA_total = length(i_set),
      eDNA_total = length(e_set),
      shared = shared,
      iDNA_unique = length(i_set) - shared,
      eDNA_unique = length(e_set) - shared,
      union = length(union(i_set, e_set))
    )
  }) |>
  mutate(
    across(c(iDNA_unique, shared, eDNA_unique), ~ 100 * .x / union, .names = "pct_{.col}"),
    label_iDNA = sprintf("%.0f%%\n(%d)", pct_iDNA_unique, iDNA_unique),
    label_shared = sprintf("%.0f%%\n(%d)", pct_shared, shared),
    label_eDNA = sprintf("%.0f%%\n(%d)", pct_eDNA_unique, eDNA_unique),
    sample_label_iDNA = paste0("n=", n_iDNA),
    sample_label_eDNA = paste0("n=", n_eDNA)
  )

write_csv(overlap, summary_path)

theta <- seq(0, 2 * pi, length.out = 241)
circle_template <- bind_rows(
  tibble(dna_type = "iDNA", x = -0.38 + 0.77 * cos(theta), y = 0.77 * sin(theta)),
  tibble(dna_type = "eDNA", x = 0.38 + 0.77 * cos(theta), y = 0.77 * sin(theta))
)

circles <- overlap |>
  select(site, depth_group) |>
  crossing(circle_template) |>
  mutate(group = interaction(site, depth_group, dna_type, drop = TRUE))

value_labels <- bind_rows(
  overlap |> transmute(site, depth_group, x = -0.72, y = 0.06, label = label_iDNA),
  overlap |> transmute(site, depth_group, x = 0, y = 0.06, label = label_shared),
  overlap |> transmute(site, depth_group, x = 0.72, y = 0.06, label = label_eDNA)
)

sample_labels <- bind_rows(
  overlap |> transmute(site, depth_group, x = -0.55, y = -0.88, label = sample_label_iDNA, dna_type = "iDNA"),
  overlap |> transmute(site, depth_group, x = 0.55, y = -0.88, label = sample_label_eDNA, dna_type = "eDNA")
)

p <- ggplot() +
  geom_polygon(
    data = circles,
    aes(x = x, y = y, group = group, fill = dna_type, colour = dna_type),
    alpha = 0.26,
    linewidth = 0.55
  ) +
  geom_text(data = value_labels, aes(x = x, y = y, label = label), size = 4.35, lineheight = 0.88) +
  geom_text(
    data = sample_labels,
    aes(x = x, y = y, label = label, colour = dna_type),
    size = 3.75,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_grid(depth_group ~ site, switch = "y") +
  scale_fill_manual(values = pool_colours, name = "DNA pool") +
  scale_colour_manual(values = pool_colours, name = "DNA pool") +
  coord_fixed(xlim = c(-1.22, 1.22), ylim = c(-1.04, 0.92), clip = "off") +
  labs(
    title = "Archaeal ASV overlap between iDNA and eDNA",
    subtitle = "Percent of the pooled site-by-depth ASV union (ASV count); n gives retained samples per DNA pool",
    x = NULL,
    y = "Depth (cm)"
  ) +
  theme_bw(base_size = 13, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey55", linewidth = 0.35),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "grey96", colour = "grey65", linewidth = 0.35),
    strip.text.x = element_text(face = "bold", size = 12.5),
    strip.text.y.left = element_text(face = "bold", size = 11.5, angle = 0),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.8, colour = "grey25", margin = margin(b = 5)),
    axis.title.y = element_text(face = "bold", size = 12.5, margin = margin(r = 7)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10.5),
    legend.key.width = unit(0.65, "cm"),
    panel.spacing = unit(0.15, "lines"),
    plot.margin = margin(5, 7, 5, 5)
  )

figure_width <- 10.0
figure_to_save <- p
central_group_gutter_add_in <- 0
if (house_style_v2_run) {
  # AZ/SG and LC/NA are the two side-by-side site groups. Preserve the narrow
  # within-group facet spacing and add only one modest group-level spacer.
  central_group_gutter_add_in <- 0.15
  figure_to_save <- ggplotGrob(p)
  panel_columns <- sort(unique(
    figure_to_save$layout$l[grepl("^panel", figure_to_save$layout$name)]
  ))
  stopifnot(length(panel_columns) == 4)
  figure_to_save <- gtable::gtable_add_cols(
    figure_to_save,
    grid::unit(central_group_gutter_add_in, "in"),
    pos = panel_columns[2]
  )
  figure_width <- figure_width + central_group_gutter_add_in
}

ggsave(
  figure_path,
  figure_to_save,
  width = figure_width,
  height = 9.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Retained samples: ", nrow(env))
message(sprintf(
  "Fig. S3 central group gutter: retained 0.15-line facet gap plus %.2f in (%.3f%% added canvas width)",
  central_group_gutter_add_in,
  100 * central_group_gutter_add_in / figure_width
))
message("Wrote vector PDF: ", figure_path)
message("Wrote overlap summary: ", summary_path)
