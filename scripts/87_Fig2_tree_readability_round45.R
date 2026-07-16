#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(grid)
  library(gridGraphics)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
house_style_run <- identical(Sys.getenv("ROUND48_HOUSESTYLE_OUTPUT"), "1")
fig2_house_style_v2_run <- house_style_run &&
  identical(Sys.getenv("ROUND48_FIG2_HOUSESTYLE_V2"), "1")
fig2_house_style_v3_run <- house_style_run &&
  identical(Sys.getenv("ROUND48_FIG2_HOUSESTYLE_V3"), "1")
fig2_spine_layout_run <- fig2_house_style_v2_run || fig2_house_style_v3_run
figure_dir <- file.path(project_root, "figures")
candidate_dir <- file.path(figure_dir, "candidates")
provisional_dir <- file.path(figure_dir, "provisional")
dir.create(provisional_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(FIG3_WRITE_ANALYSIS_OUTPUTS = "0")
Sys.setenv(FIG2_PANEL_A_NO_WRITE = "1")
source(file.path(project_root, "scripts", "24_Fig3_reframed_phylum_phylo_context.R"), local = FALSE)
Sys.unsetenv("FIG2_PANEL_A_NO_WRITE")

Sys.setenv(FIG3_TREE_SOURCE_ONLY = "1")
source(file.path(project_root, "scripts", "22_legacy_tree_refmeta_phylum_panels.R"), local = FALSE)
Sys.unsetenv("FIG3_TREE_SOURCE_ONLY")

# Sourced scripts define their own output directories. Re-establish the final
# Figure 2 destinations here so this compositor controls only its own PDF.
figure_dir <- file.path(project_root, "figures")
candidate_dir <- file.path(figure_dir, "candidates")
provisional_dir <- file.path(figure_dir, "provisional")

tree_panel_width <- 19.2
tree_panel_height <- 20.0
# Widen the composite instead of stretching either panel. At a 170-mm final
# width this 21 x 28.32 inch canvas is about 229 mm high, while the tree and
# composition panel retain their native render boxes.
figure_width <- 21.0

# Panel A's established standalone render is 8.7 x 5.6 inches. The previous
# compositor forced it into a 15.36 x 5.32 inch viewport, visibly flattening
# the facets. Scale it uniformly to 12 inches wide and retain its native ratio.
composition_native_width <- 8.7
composition_native_height <- 5.6
composition_width <- 12.0
composition_height <- composition_width * composition_native_height / composition_native_width

# A 0.60-inch empty band is 2.86% of the 21-inch figure width: visible but
# still subordinate to the outer margins, matching the Fig. 1 V10 house rule.
panel_gutter <- 0.60
figure_height <- tree_panel_height + panel_gutter + composition_height

crenarchaeota_asvs <- ann$ASV[ann$phylum_panel == "Crenarchaeota"]
crenarchaeota_reads <- rowSums(
  asv_counts[intersect(rownames(asv_counts), crenarchaeota_asvs), , drop = FALSE],
  na.rm = TRUE
)
main_asvs <- names(sort(crenarchaeota_reads, decreasing = TRUE))[
  seq_len(min(24, length(crenarchaeota_reads)))
]
main_share <- sum(crenarchaeota_reads[main_asvs]) / sum(crenarchaeota_reads)
main_panel_defs <- list(
  Crenarchaeota = list(
    title = "Abundant Crenarchaeota lineages",
    asv_filter = bquote(ASV %in% .(main_asvs))
  )
)

Sys.setenv(FIG3_MAIN_TREE = "1")
Sys.setenv(FIG2_MAIN_TREE_READABILITY_V2 = "1")
# Keep the legacy base-graphics capture off the implicit Rplots.pdf device.
# grid.echo() can otherwise leave a default device open until process exit.
grDevices::pdf(file = NULL, width = tree_panel_width, height = tree_panel_height)
capture_device <- grDevices::dev.cur()
tree_grob <- grid.grabExpr(
  {
    grid.echo(
      function() {
        plot_panel(
          "Crenarchaeota",
          file = NULL,
          defs = main_panel_defs,
          title_cols = phylum_panel_cols
        )
      },
      newpage = FALSE
    )
  },
  width = tree_panel_width,
  height = tree_panel_height,
  wrap = TRUE
)
if (capture_device %in% grDevices::dev.list()) {
  grDevices::dev.off(capture_device)
}
Sys.unsetenv("FIG2_MAIN_TREE_READABILITY_V2")
Sys.unsetenv("FIG3_MAIN_TREE")

composition_x_centered <- (figure_width - composition_width) / (2 * figure_width)
composition_y <- (tree_panel_height + panel_gutter) / figure_height
composition_w <- composition_width / figure_width
composition_h <- composition_height / figure_height
tree_x <- (figure_width - tree_panel_width) / (2 * figure_width)
tree_w <- tree_panel_width / figure_width
tree_h <- tree_panel_height / figure_height

# House-style v2/v3: Panel A and Panel B share the tree panel's established
# left edge. Both native render boxes remain unchanged.
composition_x <- if (fig2_spine_layout_run) tree_x else composition_x_centered
tag_a_x <- if (fig2_spine_layout_run) tree_x else composition_x - 0.030
tag_b_x <- if (fig2_spine_layout_run) tree_x else tree_x + 0.006

# House-style v3: align the top of the B tag to the title row's measured top
# on the actual tree capture device. The v2 position is retained for QA only;
# no guessed vertical offset is reused.
tag_b_y_v2 <- tree_h - 0.006
if (fig2_house_style_v3_run) {
  tree_metrics <- round48_tree_layout_metrics[["Crenarchaeota"]]
  if (is.null(tree_metrics$title_top_ndc) || !is.finite(tree_metrics$title_top_ndc)) {
    stop("Missing rendered title-top measurement for Fig. 2 Panel B")
  }
  tag_b_y <- tree_h * tree_metrics$title_top_ndc
} else {
  tag_b_y <- tag_b_y_v2
}

fig2 <- ggdraw() +
  draw_plot(
    p_composition,
    x = composition_x,
    y = composition_y,
    width = composition_w,
    height = composition_h
  ) +
  draw_grob(tree_grob, x = tree_x, y = 0, width = tree_w, height = tree_h) +
  draw_label(
    "A",
    x = tag_a_x,
    y = 0.994,
    hjust = 0,
    vjust = 1,
    size = 32,
    fontface = "bold",
    fontfamily = "sans",
    colour = "black"
  ) +
  draw_label(
    "B",
    x = tag_b_x,
    y = tag_b_y,
    hjust = 0,
    vjust = 1,
    size = 32,
    fontface = "bold",
    fontfamily = "sans",
    colour = "black"
  )

out_candidate <- if (fig2_house_style_v3_run) {
  file.path(candidate_dir, "Fig2_round48_housestyle_v3.pdf")
} else if (fig2_house_style_v2_run) {
  file.path(candidate_dir, "Fig2_round48_housestyle_v2.pdf")
} else if (house_style_run) {
  file.path(candidate_dir, "Fig2_round48_housestyle.pdf")
} else {
  file.path(figure_dir, "Fig2_V6_composition_plus_tree_readable_round45.pdf")
}
out_release <- file.path(
  figure_dir,
  "Fig2_composition_AOA_tree_round45.pdf"
)
out_provisional <- file.path(
  provisional_dir,
  "Fig2_composition_AOA_tree_round45.pdf"
)

ggsave(
  out_candidate,
  fig2,
  width = figure_width,
  height = figure_height,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)
if (!house_style_run) {
  file.copy(out_candidate, out_release, overwrite = TRUE)
  file.copy(out_candidate, out_provisional, overwrite = TRUE)
}

message(sprintf("Fig. 2 output dimensions: %.1f x %.1f inches", figure_width, figure_height))
message(sprintf(
  "Fig. 2 panel A native ratio retained at %.2f x %.2f inches; gutter = %.2f inches (%.3f%% of width)",
  composition_width,
  composition_height,
  panel_gutter,
  100 * panel_gutter / figure_width
))
message(sprintf(
  paste0(
    "Fig. 2 placed-box left edges (inches): A %.2f -> %.2f; ",
    "B %.2f -> %.2f; A-B offset %.2f -> %.2f"
  ),
  composition_x_centered * figure_width,
  composition_x * figure_width,
  tree_x * figure_width,
  tree_x * figure_width,
  (composition_x_centered - tree_x) * figure_width,
  (composition_x - tree_x) * figure_width
))
if (fig2_house_style_v3_run) {
  cat(sprintf(
    paste0(
      "Fig. 2B v3 optical tag alignment (canvas y from bottom): ",
      "B top %.6f npc / %.6f in before -> %.6f npc / %.6f in after; ",
      "Panel B title top %.6f npc / %.6f in; B x %.6f npc / %.6f in; ",
      "A x %.6f npc / %.6f in and A top %.6f npc / %.6f in unchanged\n"
    ),
    tag_b_y_v2,
    tag_b_y_v2 * figure_height,
    tag_b_y,
    tag_b_y * figure_height,
    tree_h * tree_metrics$title_top_ndc,
    tree_h * tree_metrics$title_top_ndc * figure_height,
    tag_b_x,
    tag_b_x * figure_width,
    tag_a_x,
    tag_a_x * figure_width,
    0.994,
    0.994 * figure_height
  ))
}
message(sprintf(
  "Fig. 2B retains %d abundant Crenarchaeota ASVs (%.1f%% of source-tree Crenarchaeota reads)",
  length(main_asvs),
  100 * main_share
))
message("Wrote vector PDF: ", out_candidate)
if (!house_style_run) {
  message("Copied vector PDF: ", out_release)
  message("Copied vector PDF: ", out_provisional)
}
