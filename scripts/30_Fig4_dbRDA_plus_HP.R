## ============================================================
## Fig. 4: dbRDA ordinations plus environmental HP panels
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggrepel)
  library(patchwork)
  library(rdacca.hp)
})

script_path <- commandArgs(FALSE) |>
  grep(pattern = "^--file=", value = TRUE) |>
  sub(pattern = "^--file=", replacement = "")
if (length(script_path) == 0) {
  script_path <- file.path("scripts", "30_Fig4_dbRDA_plus_HP.R")
}

PROJECT_DIR <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
house_style_run <- identical(Sys.getenv("ROUND48_HOUSESTYLE_OUTPUT"), "1")
house_style_v2_run <- house_style_run &&
  identical(Sys.getenv("ROUND48_FIG4_HOUSESTYLE_V2"), "1")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
CANDIDATE_DIR <- file.path(FIG_DIR, "candidates")
PROVISIONAL_DIR <- file.path(FIG_DIR, "provisional")
dir.create(PROVISIONAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CANDIDATE_DIR, recursive = TRUE, showWarnings = FALSE)
ASV_PATH <- file.path(PROJECT_DIR, "data", "arc_60cm_202rare", "asv_202.txt")
ENV_PATH <- file.path(PROJECT_DIR, "data", "env_2024.csv")
OUT_PDF <- if (house_style_v2_run) {
  file.path(CANDIDATE_DIR, "Fig4_round48_housestyle_v2.pdf")
} else if (house_style_run) {
  file.path(CANDIDATE_DIR, "Fig4_round48_housestyle.pdf")
} else {
  file.path(FIG_DIR, "04_Fig4_dbRDA_plus_HP_round42.pdf")
}
RELEASE_PDF <- file.path(FIG_DIR, "Fig4_dbRDA_HP_round42.pdf")
PROVISIONAL_PDF <- file.path(PROVISIONAL_DIR, "Fig4_dbRDA_HP_round42.pdf")

SITE_LEVELS <- c("AZ", "SG", "LC", "NB")
SITE_LABELS <- c("AZ", "SG", "LC", "NA")
DNA_LEVELS <- c("iDNA", "eDNA")
DEPTH_60_LABELS <- c("0-5", "5-10", "10-20", "20-40", "40-60")

env_vars <- c(
  "pH", "Conductivity", "moisture", "CN",
  "Feo", "Alo", "Mno", "Sio",
  "NH4", "NO3", "Po", "Pi"
)

theme_clean <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.45),
      axis.text = element_text(colour = "grey20"),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
}

label_for_parse <- function(variable) {
  case_when(
    variable == "Feo" ~ "Fe[ox]",
    variable == "Alo" ~ "Al[ox]",
    variable == "Mno" ~ "Mn[ox]",
    variable == "Sio" ~ "Si[ox]",
    variable == "NH4" ~ "NH[4]^'+'",
    variable == "NO3" ~ "NO[3]^'-'",
    variable == "Pi" ~ "P[i]",
    variable == "Po" ~ "P[o]",
    variable == "pH" ~ "'pH'",
    variable == "CN" ~ "'C/N'",
    variable == "moisture" ~ "'moisture'",
    variable == "Conductivity" ~ "'Conductivity'",
    TRUE ~ paste0("'", variable, "'")
  )
}

asv_raw <- read.table(
  ASV_PATH,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  sep = "\t"
)
env_raw <- read.csv(ENV_PATH, header = TRUE, check.names = FALSE)

env_60 <- env_raw %>%
  filter(depth < 60, sample_id %in% colnames(asv_raw)) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    site_label = factor(site, levels = SITE_LEVELS, labels = SITE_LABELS),
    dna_type = factor(dna_type, levels = DNA_LEVELS),
    depth_order = as.numeric(depth_order),
    depth_mid = as.numeric(depth),
    depth_label = factor(
      depth_order,
      levels = 1:5,
      labels = paste(DEPTH_60_LABELS, "cm")
    )
  ) %>%
  arrange(match(sample_id, colnames(asv_raw)))

asv <- asv_raw[, env_60$sample_id, drop = FALSE]
asv <- asv[rowSums(asv) > 0, colSums(asv) > 0, drop = FALSE]
env_60 <- env_60 %>% filter(sample_id %in% colnames(asv))
asv <- asv[, env_60$sample_id, drop = FALSE]
stopifnot(identical(colnames(asv), env_60$sample_id))

run_pool_dbrda <- function(pool) {
  env_pool <- env_60 %>%
    filter(dna_type == pool) %>%
    select(all_of(env_vars), site, site_label, dna_type, depth_order,
           depth_mid, depth_label, sample_id) %>%
    na.omit()

  asv_pool <- asv[, env_pool$sample_id, drop = FALSE]
  asv_pool <- asv_pool[rowSums(asv_pool) > 0, , drop = FALSE]
  asv_pool_hell <- decostand(t(asv_pool), method = "hellinger")

  mod <- dbrda(
    as.formula(paste("asv_pool_hell ~", paste(env_vars, collapse = " + "))),
    data = env_pool,
    distance = "bray",
    scale = TRUE
  )

  site_scores <- scores(mod, display = "sites", choices = 1:2, scaling = 1) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    left_join(
      env_pool %>% select(sample_id, site_label, depth_label),
      by = "sample_id"
    )

  bp_scores <- scores(mod, display = "bp", choices = 1:2, scaling = 1) %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    mutate(
      label = label_for_parse(variable)
    )

  list(
    pool = pool,
    model = mod,
    sites = site_scores,
    biplot = bp_scores,
    axis_pct = summary(mod)$cont$importance["Proportion Explained", 1:2] * 100,
    r2 = RsquareAdj(mod),
    n = nrow(env_pool)
  )
}

make_dbrda_panel <- function(db, panel_label, title) {
  arrow_mul <- 1.65
  stat_label <- sprintf(
    "paste(R^2, ' = %.3f; adj. ', R^2, ' = %.3f; ', italic(n), ' = %d')",
    db$r2$r.squared,
    db$r2$adj.r.squared,
    db$n
  )
  bp <- db$biplot %>%
    mutate(
      dbRDA1 = dbRDA1 * arrow_mul,
      dbRDA2 = dbRDA2 * arrow_mul
    )

  ggplot(db$sites, aes(dbRDA1, dbRDA2)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey78") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey78") +
    geom_point(aes(colour = site_label, alpha = depth_label), size = 2.7) +
    geom_segment(
      data = bp,
      aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
      inherit.aes = FALSE,
      arrow = arrow(angle = 18, length = unit(0.22, "cm"), type = "closed"),
      linewidth = 0.55,
      colour = "grey10"
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = stat_label,
      parse = TRUE,
      hjust = -0.06,
      vjust = 1.25,
      size = 3.35,
      family = "Helvetica",
      colour = "grey15"
    ) +
    geom_text_repel(
      data = bp,
      aes(dbRDA1, dbRDA2, label = label),
      inherit.aes = FALSE,
      size = 3.45,
      fontface = "bold",
      parse = TRUE,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.1,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE,
      colour = "#1F2A70"
    ) +
    scale_colour_manual(
      values = c("AZ" = "#F8766D", "SG" = "#7CAE00", "LC" = "#00BFC4", "NA" = "#C77CFF"),
      breaks = SITE_LABELS,
      name = "Site"
    ) +
    scale_alpha_manual(
      values = c(
        "0-5 cm" = 0.30, "5-10 cm" = 0.46, "10-20 cm" = 0.62,
        "20-40 cm" = 0.80, "40-60 cm" = 1.00
      ),
      name = "Depth"
    ) +
    guides(
      alpha = guide_legend(
        order = 1,
        title.position = "left",
        nrow = 1,
        override.aes = list(colour = "grey25", size = 3.2)
      ),
      colour = guide_legend(
        order = 2,
        title.position = "left",
        nrow = 1,
        override.aes = list(alpha = 1, size = 3.2)
      )
    ) +
    labs(
      x = sprintf("dbRDA1 (%.1f%%)", db$axis_pct[1]),
      y = sprintf("dbRDA2 (%.1f%%)", db$axis_pct[2]),
      title = title,
      tag = panel_label
    ) +
    theme_clean(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12.5),
      plot.tag = element_text(
        face = "bold", size = 18, family = "Helvetica", colour = "black"
      ),
      plot.tag.position = c(0, 1.015),
      plot.margin = margin(t = 10, r = 6, b = 4, l = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.spacing.y = unit(1.5, "pt"),
      legend.spacing.x = unit(8, "pt"),
      legend.margin = margin(t = -4),
      legend.title = element_text(face = "bold", size = 9.8),
      legend.text = element_text(size = 9.8)
    )
}

run_pool_hp <- function(pool) {
  env_pool <- env_60 %>%
    filter(dna_type == pool) %>%
    select(sample_id, all_of(env_vars)) %>%
    na.omit()

  asv_pool <- asv[, env_pool$sample_id, drop = FALSE]
  asv_pool <- asv_pool[rowSums(asv_pool) > 0, , drop = FALSE]
  asv_pool_hell <- decostand(t(asv_pool), method = "hellinger")
  bray_dist <- vegdist(asv_pool_hell, method = "bray")

  hp <- rdacca.hp(
    bray_dist,
    env_pool %>% select(all_of(env_vars)),
    method = "dbRDA",
    type = "adjR2"
  )

  as.data.frame(hp$Hier.part) %>%
    rownames_to_column("variable") %>%
    as_tibble() %>%
    transmute(
      variable,
      individual = Individual,
      percent = `I.perc(%)`,
      highlighted = FALSE
    ) %>%
    arrange(desc(percent)) %>%
    mutate(variable = factor(variable, levels = variable))
}

make_hp_panel <- function(hp_tbl, panel_label, title) {
  ggplot(hp_tbl, aes(variable, percent, fill = highlighted)) +
    geom_col(width = 0.82, colour = "grey20", linewidth = 0.18) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey35") +
    scale_fill_manual(values = c(`TRUE` = "grey65", `FALSE` = "grey65"), guide = "none") +
    scale_x_discrete(
      labels = function(x) parse(text = label_for_parse(x))
    ) +
    scale_y_continuous(
      limits = c(min(-1, hp_tbl$percent, na.rm = TRUE), 27),
      breaks = seq(0, 25, 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Variables",
      y = expression("% Individual effect to " * R^2 * " (%)"),
      title = title,
      tag = panel_label
    ) +
    theme_clean(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12.5),
      plot.tag = element_text(
        face = "bold", size = 18, family = "Helvetica", colour = "black"
      ),
      plot.tag.position = c(0, 1.015),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.title = element_text(size = 12),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 5)),
      plot.margin = margin(t = 10, r = 8, b = 4, l = 10)
    )
}

db_i <- run_pool_dbrda("iDNA")
db_e <- run_pool_dbrda("eDNA")
hp_i <- run_pool_hp("iDNA")
hp_e <- run_pool_hp("eDNA")

p_a <- make_dbrda_panel(db_i, "A", "iDNA dbRDA")
p_b <- make_dbrda_panel(db_e, "B", "eDNA dbRDA")
p_c <- make_hp_panel(hp_i, "C", "iDNA environmental HP")
p_d <- make_hp_panel(hp_e, "D", "eDNA environmental HP")

shared_figure_theme <-
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y = unit(2, "pt")
  )

if (house_style_v2_run) {
  # Four independent analyses receive explicit spacer cells. The canvas grows
  # enough to preserve the established panel-frame widths after patchwork's
  # guide collection and spacer-column allocation.
  horizontal_gutter_in <- 0.36
  vertical_gutter_in <- 0.25
  figure_width <- 11.25
  figure_height <- 9.9 + vertical_gutter_in
  column_unit_in <- (figure_width - horizontal_gutter_in) / 2
  row_unit_in <- (figure_height - vertical_gutter_in) / (1 + 0.98)
  fig4_design <- c(
    area(t = 1, l = 1, b = 1, r = 1),
    area(t = 1, l = 3, b = 1, r = 3),
    area(t = 3, l = 1, b = 3, r = 1),
    area(t = 3, l = 3, b = 3, r = 3)
  )
  fig4 <- (
    p_a + p_b + p_c + p_d +
      plot_layout(
        design = fig4_design,
        widths = c(1, horizontal_gutter_in / column_unit_in, 1),
        heights = c(1, vertical_gutter_in / row_unit_in, 0.98),
        guides = "collect"
      )
  ) & shared_figure_theme
} else {
  horizontal_gutter_in <- 0
  vertical_gutter_in <- 0
  figure_width <- 10.8
  figure_height <- 9.9
  fig4 <- (
    (p_a + p_b) / (p_c + p_d) +
      plot_layout(heights = c(1, 0.98), guides = "collect")
  ) & shared_figure_theme
}

ggsave(
  OUT_PDF,
  fig4,
  width = figure_width,
  height = figure_height,
  device = cairo_pdf,
  bg = "white"
)
if (!house_style_run) {
  file.copy(OUT_PDF, RELEASE_PDF, overwrite = TRUE)
  file.copy(OUT_PDF, PROVISIONAL_PDF, overwrite = TRUE)
}

message(sprintf(
  paste0(
    "Fig. 4 explicit gutters: horizontal %.2f in (%.3f%% of canvas); ",
    "vertical %.2f in (%.3f%% of canvas)"
  ),
  horizontal_gutter_in,
  100 * horizontal_gutter_in / figure_width,
  vertical_gutter_in,
  100 * vertical_gutter_in / figure_height
))
message("Wrote combined Figure 4 PDF to: ", OUT_PDF)
if (!house_style_run) {
  message("Copied combined Figure 4 PDF to: ", RELEASE_PDF)
  message("Copied combined Figure 4 PDF to: ", PROVISIONAL_PDF)
}
