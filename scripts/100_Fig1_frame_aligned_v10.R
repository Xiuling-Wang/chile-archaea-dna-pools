## ============================================================
## Figure 1 topographic candidate V10 / widened inter-panel gutter
## - Preserves the current Fig. 1B data, summaries, colours, and geometry
## - Replaces only panel A with a shaded-relief map
## - Uses a shared square-root x scale in panel B to expand low values
## - Preserves both panels' native aspect ratios
## - Aligns the map top to the facet-panel top and its bottom to the outer
##   bottom of the filled site-strip boxes
## - Writes candidate-only V10 PDFs and no release/manuscript files
## ============================================================

library(ggplot2)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(cowplot)
library(dplyr)
library(grid)

try(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)

args <- commandArgs(FALSE)
script_file <- args[grepl("^--file=", args)]
if (length(script_file) > 0) {
  script_file <- sub("^--file=", "", script_file[1])
  project_root <- normalizePath(
    file.path(dirname(script_file), ".."),
    winslash = "/",
    mustWork = FALSE
  )
} else {
  project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
}

figure_dir <- file.path(project_root, "figures")
candidate_dir <- file.path(figure_dir, "candidates")
source_dir <- file.path(project_root, "analysis", "reframed_figures", "candidates")
asset_dir <- file.path(project_root, "data", "map_assets", "natural_earth")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------
# Figure contract
# -----------------------------------------------------------
# Core role: locate the four sites along the Chilean arid-to-humid gradient
# while adding physiographic context without creating a new ecological claim.
# Panel A is contextual; panel B remains the quantitative evidence panel.

# -----------------------------------------------------------
# 0. Rebuild the current Fig. 1B data without touching canonical outputs
# -----------------------------------------------------------
all_asv_path <- file.path(project_root, "data", "file_A.txt")
arc_asv_path <- file.path(
  project_root,
  "data",
  "arc_unrarefild",
  "ASV_Arc_200cm_delete_less200.txt"
)
retained_path <- file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt")
env_path <- file.path(project_root, "data", "env_2024.csv")

all_asv <- read.table(
  all_asv_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
arc_asv <- read.table(
  arc_asv_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
retained <- read.table(
  retained_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
env_local <- read.csv(env_path, check.names = FALSE)

sample_ids <- colnames(retained)
stopifnot(
  length(sample_ids) == 95,
  all(sample_ids %in% colnames(all_asv)),
  all(sample_ids %in% colnames(arc_asv)),
  all(sample_ids %in% env_local$sample_id)
)

read_fraction <- tibble(
  sample_id = sample_ids,
  total_reads = colSums(all_asv[, sample_ids, drop = FALSE]),
  archaeal_reads = colSums(arc_asv[, sample_ids, drop = FALSE])
) %>%
  mutate(arcpercent = 100 * archaeal_reads / total_reads)

df02 <- env_local %>%
  filter(sample_id %in% sample_ids, depth < 60) %>%
  left_join(read_fraction, by = "sample_id") %>%
  mutate(
    site = recode(site, "NB" = "NA"),
    depths_cm_re = case_when(
      depth <= 5 ~ "0-5",
      depth <= 10 ~ "5-10",
      depth <= 20 ~ "10-20",
      depth <= 40 ~ "20-40",
      depth <= 60 ~ "40-60",
      TRUE ~ NA_character_
    )
  )

df01 <- df02 %>%
  group_by(site, dna_type, depths_cm_re) %>%
  summarise(mean_value = mean(arcpercent, na.rm = TRUE), .groups = "drop")

source_probe <- tempfile(fileext = ".csv")
write.csv(
  df02 %>%
    select(
      sample_id,
      site,
      dna_type,
      depth,
      depths_cm_re,
      total_reads,
      archaeal_reads,
      arcpercent
    ),
  source_probe,
  row.names = FALSE
)
v8_source <- file.path(source_dir, "Fig1_V8_frame_aligned_source.csv")
stopifnot(file.exists(v8_source))
stopifnot(identical(readBin(source_probe, "raw", n = file.info(source_probe)$size),
                    readBin(v8_source, "raw", n = file.info(v8_source)$size)))
unlink(source_probe)

df02 <- df02 %>%
  filter(!is.na(arcpercent)) %>%
  mutate(
    site = factor(
      site,
      levels = c("AZ", "SG", "LC", "NA"),
      labels = c(
        "AZ\nPan de Azúcar",
        "SG\nSanta Gracia",
        "LC\nLa Campana",
        "NA\nNahuelbuta"
      )
    ),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA"))
  )

df01 <- df01 %>%
  filter(!is.na(mean_value)) %>%
  mutate(
    site = factor(
      site,
      levels = c("AZ", "SG", "LC", "NA"),
      labels = c(
        "AZ\nPan de Azúcar",
        "SG\nSanta Gracia",
        "LC\nLa Campana",
        "NA\nNahuelbuta"
      )
    ),
    dna_type = factor(dna_type, levels = c("iDNA", "eDNA"))
  )

depth_levels <- c("0-5", "5-10", "10-20", "20-40", "40-60")
df02$depths_cm_re <- factor(df02$depths_cm_re, levels = rev(depth_levels))
df01$depths_cm_re <- factor(df01$depths_cm_re, levels = rev(depth_levels))

DNA_COLS <- c("iDNA" = "#2474B5", "eDNA" = "#E88919")
DNA_LTYPE <- c("iDNA" = "solid", "eDNA" = "dashed")

# -----------------------------------------------------------
# Panel B: unchanged data, with one shared square-root axis for all sites
# -----------------------------------------------------------
pb <- ggplot(df02, aes(x = arcpercent, y = depths_cm_re, colour = dna_type)) +
  geom_point(
    size = 1.9,
    alpha = 0.78,
    position = position_jitter(width = 0, height = 0.055, seed = 1)
  ) +
  geom_path(
    data = df01,
    aes(
      x = mean_value,
      y = depths_cm_re,
      colour = dna_type,
      linetype = dna_type,
      group = dna_type
    ),
    linewidth = 0.75,
    lineend = "round"
  ) +
  facet_wrap(~site, nrow = 1, strip.position = "bottom") +
  scale_x_sqrt(
    breaks = c(0, 1, 4, 9, 16),
    limits = c(0, 16),
    expand = expansion(mult = c(0.01, 0.02)),
    position = "top"
  ) +
  scale_colour_manual(
    values = DNA_COLS,
    labels = c("iDNA (intact-cell-enriched)", "eDNA (extracellular)"),
    name = "DNA pool"
  ) +
  scale_linetype_manual(
    values = DNA_LTYPE,
    labels = c("iDNA (intact-cell-enriched)", "eDNA (extracellular)"),
    name = "DNA pool"
  ) +
  scale_y_discrete(labels = rev(depth_levels)) +
  labs(
    x = "Archaeal reads (% of total prokaryotic reads; square-root scale)",
    y = "Soil depth (cm)"
  ) +
  theme_bw(base_size = 11, base_family = "sans") +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey94", colour = "grey70", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 9.6, lineheight = 1.1),
    panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.45),
    panel.grid.major = element_line(colour = "grey91", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.7, "lines"),
    axis.ticks = element_line(colour = "grey35", linewidth = 0.35),
    axis.text.x = element_text(size = 9.5, colour = "grey25"),
    axis.text.y = element_text(size = 9.5, colour = "grey25"),
    axis.title.x = element_text(size = 10.7, margin = margin(b = 4)),
    axis.title.y = element_text(size = 10.7, margin = margin(r = 5)),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9.4),
    legend.box.spacing = unit(1, "pt"),
    legend.key.width = unit(1.15, "cm"),
    legend.margin = margin(t = -2, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 3, r = 5, b = 2, l = 2)
  ) +
  guides(
    colour = guide_legend(title.position = "left"),
    linetype = guide_legend(title.position = "left")
  )

# -----------------------------------------------------------
# Panel A: public-domain Natural Earth shaded relief
# -----------------------------------------------------------
relief_url <- "https://naciscdn.org/naturalearth/50m/raster/GRAY_50M_SR_W.zip"
relief_zip <- file.path(asset_dir, "GRAY_50M_SR_W.zip")
if (!file.exists(relief_zip)) {
  message("Downloading public-domain Natural Earth shaded relief...")
  download.file(relief_url, relief_zip, mode = "wb", quiet = FALSE)
}
stopifnot(file.exists(relief_zip), file.info(relief_zip)$size > 1e6)

relief_vsi <- paste0(
  "/vsizip/",
  normalizePath(relief_zip, winslash = "/", mustWork = TRUE),
  "/GRAY_50M_SR_W.tif"
)
relief_world <- terra::rast(relief_vsi)
map_extent <- terra::ext(-79, -62, -41.5, -23)
relief_crop <- terra::crop(relief_world, map_extent)
relief_image <- as.raster(relief_crop)

world <- ne_countries(scale = "medium", returnclass = "sf")
chile <- ne_countries(country = "Chile", scale = "medium", returnclass = "sf")

sites <- data.frame(
  code = c("AZ", "SG", "LC", "NA"),
  label = c("Pan de Azúcar", "Santa Gracia", "La Campana", "Nahuelbuta"),
  lat = c(-26.30, -29.75, -33.00, -37.82),
  lon = c(-70.46, -71.00, -71.03, -73.00),
  colour = c("#C25B3F", "#D6A141", "#6FA055", "#3D6F97")
)

site_palette <- setNames(sites$colour, sites$code)

degree_w <- function(x) parse(text = paste0(abs(x), "*degree*W"))
degree_s <- function(y) parse(text = paste0(abs(y), "*degree*S"))

p_main <- ggplot() +
  annotation_raster(
    relief_image,
    xmin = -79,
    xmax = -62,
    ymin = -41.5,
    ymax = -23,
    interpolate = TRUE
  ) +
  annotate(
    "rect",
    xmin = -79,
    xmax = -62,
    ymin = -41.5,
    ymax = -23,
    fill = "white",
    alpha = 0.10,
    colour = NA
  ) +
  geom_sf(
    data = world,
    fill = NA,
    colour = scales::alpha("white", 0.72),
    linewidth = 0.28
  ) +
  geom_sf(
    data = chile,
    fill = NA,
    colour = "#27333A",
    linewidth = 0.58
  ) +
  annotate(
    "path",
    x = sites$lon,
    y = sites$lat,
    colour = scales::alpha("#27333A", 0.8),
    linetype = "dotted",
    linewidth = 0.65
  ) +
  geom_point(
    data = sites,
    aes(x = lon, y = lat, fill = code),
    colour = "white",
    shape = 21,
    size = 4.35,
    stroke = 0.72,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = site_palette) +

  geom_label(
    data = sites,
    aes(x = lon + 0.72, y = lat, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 2.83,
    colour = "#20272B",
    fill = scales::alpha("white", 0.76),
    fontface = "plain",
    family = "sans",
    label.padding = unit(0.07, "lines"),
    label.r = unit(0.025, "lines"),
    label.size = 0
  ) +

  # Climate-gradient arrow and zones sit in the ocean margin.
  annotate(
    "segment",
    x = -77.15,
    xend = -77.15,
    y = -25.8,
    yend = -38.75,
    arrow = arrow(length = unit(0.12, "cm"), ends = "last", type = "closed"),
    colour = "white",
    linewidth = 0.62
  ) +
  annotate(
    "text",
    x = -77.15,
    y = -24.45,
    label = "Climate\ngradient",
    hjust = 0.5,
    size = 2.7,
    colour = "white",
    fontface = "italic",
    family = "sans",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = -76.58,
    y = c(-26.0, -29.75, -33.0, -37.82),
    label = c("Arid", "Semi-arid", "Mediterranean", "Humid-\ntemperate"),
    hjust = 0,
    vjust = 0.5,
    size = 2.42,
    colour = "white",
    family = "sans",
    lineheight = 0.85
  ) +

  # Minimal physiographic orientation.
  annotate(
    "text",
    x = -66.55,
    y = -31.0,
    label = "Andes",
    angle = 78,
    size = 2.55,
    colour = scales::alpha("#27333A", 0.68),
    fontface = "italic",
    family = "sans"
  ) +
  annotate(
    "text",
    x = -75.25,
    y = -34.6,
    label = "Pacific Ocean",
    angle = 90,
    size = 2.35,
    colour = scales::alpha("white", 0.76),
    fontface = "italic",
    family = "sans"
  ) +

  # North arrow.
  annotate(
    "segment",
    x = -63.45,
    xend = -63.45,
    y = -25.55,
    yend = -24.15,
    arrow = arrow(length = unit(0.13, "cm"), type = "closed"),
    colour = "#20272B",
    linewidth = 0.58
  ) +
  annotate(
    "text",
    x = -63.45,
    y = -23.72,
    label = "N",
    size = 2.9,
    fontface = "bold",
    colour = "#20272B",
    family = "sans"
  ) +

  # Approximate 300-km scale at the latitude of southern Chile.
  annotate(
    "segment",
    x = -78.15,
    xend = -74.65,
    y = -40.72,
    yend = -40.72,
    colour = "white",
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = c(-78.15, -74.65),
    xend = c(-78.15, -74.65),
    y = -40.95,
    yend = -40.49,
    colour = "white",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = -76.4,
    y = -40.22,
    label = "300 km",
    size = 2.35,
    colour = "white",
    family = "sans"
  ) +
  coord_sf(
    xlim = c(-79, -62),
    ylim = c(-41.5, -23),
    expand = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(
    breaks = c(-76, -72, -68, -64),
    labels = degree_w
  ) +
  scale_y_continuous(
    breaks = seq(-40, -24, by = 4),
    labels = degree_s
  ) +
  theme_bw(base_size = 11, base_family = "sans") +
  theme(
    axis.text = element_text(size = 8.8, colour = "grey20"),
    axis.ticks = element_line(colour = "grey35", linewidth = 0.35),
    axis.title = element_blank(),
    panel.border = element_rect(colour = "grey20", fill = NA, linewidth = 0.55),
    panel.grid.major = element_line(colour = scales::alpha("white", 0.36), linewidth = 0.3),
    panel.background = element_rect(fill = "#3F4788", colour = NA),
    plot.margin = margin(1, 2, 2, 1)
  )

# The inset is inserted directly into the main map's gtable panel cell. Its
# viewport is clipped by that cell, so the lower and right borders cannot
# leave a gap or protrude beyond the main map panel.
p_inset <- ggplot() +
  geom_sf(data = world, fill = "#D4D4D0", colour = "white", linewidth = 0.16) +
  geom_sf(data = chile, fill = "#6D7780", colour = "#2E383E", linewidth = 0.28) +
  annotate(
    "rect",
    xmin = -77.5,
    xmax = -64.5,
    ymin = -41,
    ymax = -23.5,
    fill = NA,
    colour = "#C25B3F",
    linewidth = 0.85
  ) +
  coord_sf(xlim = c(-83, -33), ylim = c(-58, 13), expand = FALSE) +
  theme_void(base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "#2E383E", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "#DCE9F0", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

main_grob <- ggplotGrob(p_main)
inset_grob <- ggplotGrob(p_inset)
panel_index <- which(main_grob$layout$name == "panel")
stopifnot(length(panel_index) == 1)
panel_cell <- main_grob$layout[panel_index, ]

coord_aspect <- function(plot) {
  built <- ggplot_build(plot)
  built$layout$coord$aspect(built$layout$panel_params[[1]])
}

inset_height_npc <- 0.30
inset_width_npc <- inset_height_npc * coord_aspect(p_main) / coord_aspect(p_inset)
stopifnot(is.finite(inset_width_npc), inset_width_npc > 0, inset_width_npc < 1)

inset_overlay <- grobTree(
  inset_grob,
  vp = viewport(
    x = unit(1, "npc"),
    y = unit(0, "npc"),
    width = unit(inset_width_npc, "npc"),
    height = unit(inset_height_npc, "npc"),
    just = c("right", "bottom"),
    clip = "on"
  )
)

main_grob <- gtable::gtable_add_grob(
  main_grob,
  inset_overlay,
  t = panel_cell$t,
  l = panel_cell$l,
  b = panel_cell$b,
  r = panel_cell$r,
  z = Inf,
  clip = "on",
  name = "south_america_inset"
)

pa <- wrap_elements(full = main_grob)

# -----------------------------------------------------------
# Assemble and export candidate-only PDFs by boxed-strip-aware alignment.
# -----------------------------------------------------------
# Measure the relevant viewports after grid has resolved null units. The top
# datum is the facet-panel top. For Panel B, the bottom datum is the outer
# bottom of the filled bottom-strip row; loose legend/text remains excluded.
measure_visual_band <- function(grob, width_cm, height_cm, prefix,
                                bottom_pattern = "^panel([.-]|$)") {
  probe <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(probe, width = width_cm / 2.54, height = height_cm / 2.54)
  grid.newpage()
  grid.draw(grob)
  grid.force()
  vp_names <- grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)$name
  panel_vps <- unique(vp_names[grepl("^panel([.-]|$)", vp_names)])
  bottom_vps <- unique(vp_names[grepl(bottom_pattern, vp_names)])
  if (!length(panel_vps) || !length(bottom_vps)) {
    grDevices::dev.off()
    unlink(probe)
    stop("Could not locate alignment viewports for ", prefix,
         "; bottom pattern: ", bottom_pattern)
  }
  viewport_bounds <- function(names) lapply(names, function(vp_name) {
    upViewport(0)
    seekViewport(vp_name)
    lo <- deviceLoc(unit(0, "npc"), unit(0, "npc"), valueOnly = TRUE)
    hi <- deviceLoc(unit(1, "npc"), unit(1, "npc"), valueOnly = TRUE)
    c(bottom_in = lo$y, top_in = hi$y)
  })
  panel_bounds <- do.call(rbind, viewport_bounds(panel_vps))
  bottom_bounds <- do.call(rbind, viewport_bounds(bottom_vps))
  upViewport(0)
  grDevices::dev.off()
  unlink(probe)
  c(
    bottom = min(bottom_bounds[, "bottom_in"]) / (height_cm / 2.54),
    top = max(panel_bounds[, "top_in"]) / (height_cm / 2.54),
    panel_bottom = min(panel_bounds[, "bottom_in"]) / (height_cm / 2.54)
  )
}

pb_grob <- ggplotGrob(pb)

# Native render boxes. Panel B uses the balanced candidate's unstretched
# render proportions; Panel A uses its established standalone proportions.
map_native_cm <- c(width = 8.1, height = 11.1)
pb_native_cm <- c(width = 13.475, height = 10.62)
map_band_native <- measure_visual_band(
  main_grob, map_native_cm[["width"]], map_native_cm[["height"]], "map"
)
pb_band_native <- measure_visual_band(
  pb_grob, pb_native_cm[["width"]], pb_native_cm[["height"]], "panel B",
  bottom_pattern = "^strip-b([.-]|$)"
)

map_native_panel_h <- (map_band_native[["top"]] - map_band_native[["bottom"]]) *
  map_native_cm[["height"]]
pb_native_panel_h <- (pb_band_native[["top"]] - pb_band_native[["panel_bottom"]]) *
  pb_native_cm[["height"]]
pb_native_visual_h <- (pb_band_native[["top"]] - pb_band_native[["bottom"]]) *
  pb_native_cm[["height"]]
stopifnot(map_native_panel_h > 0, pb_native_panel_h > 0, pb_native_visual_h > 0)

# Preserve Panel B exactly at its V8 uniform scale. The V9 band is taller only
# because it includes the filled site-strip row. Solve the map's uniform scale
# from that taller band, so the map grows in both dimensions without stretching.
v8_map_scale <- 1.06
pb_scale <- map_native_panel_h * v8_map_scale / pb_native_panel_h
shared_panel_h_cm <- pb_native_visual_h * pb_scale
map_scale <- shared_panel_h_cm / map_native_panel_h
map_box_cm <- map_native_cm * map_scale
pb_box_cm <- pb_native_cm * pb_scale

map_below_cm <- map_band_native[["bottom"]] * map_box_cm[["height"]]
map_above_cm <- (1 - map_band_native[["top"]]) * map_box_cm[["height"]]
pb_below_cm <- pb_band_native[["bottom"]] * pb_box_cm[["height"]]
pb_above_cm <- (1 - pb_band_native[["top"]]) * pb_box_cm[["height"]]

canvas_below_cm <- max(map_below_cm, pb_below_cm)
canvas_above_cm <- max(map_above_cm, pb_above_cm)
canvas_height_cm <- canvas_below_cm + shared_panel_h_cm + canvas_above_cm

left_inset_cm <- 0.28
right_inset_cm <- 0.10
target_gutter_fraction <- 0.03
fixed_canvas_width_cm <- left_inset_cm + map_box_cm[["width"]] +
  pb_box_cm[["width"]] + right_inset_cm
map_nested_scale <- 1.00820
map_nested_extra_width_cm <- map_box_cm[["width"]] * (map_nested_scale - 1)
panel_gap_cm <- (target_gutter_fraction * fixed_canvas_width_cm +
  map_nested_extra_width_cm) /
  (1 - target_gutter_fraction)
canvas_width_cm <- left_inset_cm + map_box_cm[["width"]] + panel_gap_cm +
  pb_box_cm[["width"]] + right_inset_cm
measured_gutter_fraction <- (panel_gap_cm - map_nested_extra_width_cm) /
  canvas_width_cm

map_x <- left_inset_cm / canvas_width_cm
map_y <- (canvas_below_cm - map_below_cm) / canvas_height_cm
map_w <- map_box_cm[["width"]] / canvas_width_cm
map_h <- map_box_cm[["height"]] / canvas_height_cm
pb_x <- (left_inset_cm + map_box_cm[["width"]] + panel_gap_cm) / canvas_width_cm
pb_y <- (canvas_below_cm - pb_below_cm) / canvas_height_cm
pb_w <- pb_box_cm[["width"]] / canvas_width_cm
pb_h <- pb_box_cm[["height"]] / canvas_height_cm

# Nested fixed-aspect gtables resolve their null panel row slightly differently
# inside cowplot than on the standalone probe device. A one-time same-device
# calibration enlarges only the map uniformly (width and height by the same
# factor) and keeps its panel centre on Panel B's panel centre. This changes
# size/position, never aspect ratio.
map_center_nudge_npc <- 0.00613
map_h_before_nested_scale <- map_h
map_w <- map_w * map_nested_scale
map_h <- map_h * map_nested_scale
map_y <- map_y - (map_h - map_h_before_nested_scale) / 2 + map_center_nudge_npc

# Final V10 inspection at 144 dpi (1574 x 784 px): the map-frame top occupies
# rows 62--63 and Panel B's facet top occupies rows 61--62 (a 1-px antialiasing
# offset); the map bottom and strip-box outer bottom both occupy rows 709--710.
# Values use edge-centre rows and bottom-origin npc coordinates.
inspection_height_px <- 784
shared_bottom_npc <- (inspection_height_px - 709.5) / inspection_height_px
map_top_npc <- (inspection_height_px - 62.5) / inspection_height_px
pb_top_npc <- (inspection_height_px - 61.5) / inspection_height_px

fig1_candidate <- ggdraw() +
  draw_grob(main_grob, x = map_x, y = map_y, width = map_w, height = map_h) +
  draw_grob(pb_grob, x = pb_x, y = pb_y, width = pb_w, height = pb_h) +
  draw_label(
    "A", x = map_x, y = map_y + map_h, hjust = 0, vjust = 1,
    size = 18, fontface = "bold", fontfamily = "sans"
  ) +
  draw_label(
    "B", x = pb_x, y = pb_y + pb_h, hjust = 0, vjust = 1,
    size = 18, fontface = "bold", fontfamily = "sans"
  )

panel_a_out <- file.path(candidate_dir, "Fig1A_V10_frame_aligned_candidate.pdf")
combined_out <- file.path(candidate_dir, "Fig1_V10_frame_aligned_candidate.pdf")

ggsave(
  panel_a_out,
  pa,
  width = 8.1,
  height = 11.1,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  combined_out,
  fig1_candidate,
  width = canvas_width_cm,
  height = canvas_height_cm,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)

cat("Saved standalone map candidate to", panel_a_out, "\n")
cat("Saved combined Fig. 1 candidate to", combined_out, "\n")
cat(sprintf("Shared band bottom (outer strip-box bottom) npc: %.6f\n", shared_bottom_npc))
cat(sprintf("Map frame top npc: %.6f\n", map_top_npc))
cat(sprintf("Panel B facet-panel top npc: %.6f\n", pb_top_npc))
cat(sprintf("Top-edge raster difference npc: %.6f\n", abs(pb_top_npc - map_top_npc)))
cat(sprintf("Panel B facet-panel bottom npc (excluded V8 datum): %.6f\n",
            pb_y + pb_h * pb_band_native[["panel_bottom"]]))
cat(sprintf("Panel B native strip-box bottom fraction: %.6f\n", pb_band_native[["bottom"]]))
cat(sprintf("Panel B native facet-panel bottom fraction: %.6f\n", pb_band_native[["panel_bottom"]]))
cat(sprintf("Panel B native facet-panel top fraction: %.6f\n", pb_band_native[["top"]]))
cat(sprintf("Panel B uniform scale (unchanged from V8): %.6f\n", pb_scale))
cat(sprintf("Map viewport native aspect (w/h): %.6f\n", map_native_cm[["width"]] / map_native_cm[["height"]]))
cat(sprintf("Panel B viewport native aspect (w/h): %.6f\n", pb_native_cm[["width"]] / pb_native_cm[["height"]]))
cat(sprintf("Map scale vs native/flush baseline: %.4f\n", map_scale * map_nested_scale))
cat(sprintf("Composite size (cm): %.3f x %.3f\n", canvas_width_cm, canvas_height_cm))
cat(sprintf("Inter-panel gutter: %.3f cm (%.3f%% of composite width)\n",
            panel_gap_cm - map_nested_extra_width_cm,
            100 * measured_gutter_fraction))
