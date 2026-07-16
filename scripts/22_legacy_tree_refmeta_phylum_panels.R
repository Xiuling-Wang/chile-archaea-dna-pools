#!/usr/bin/env Rscript

# Build taxonomic zoom-in panels and iTOL-compatible annotation files for the
# legacy archaeal RAxML/ARB tree. This script emphasizes ASV occurrence across
# Site x DNA-pool groups rather than only specialist status.

suppressPackageStartupMessages({
  library(ape)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
house_style_run <- identical(Sys.getenv("ROUND48_HOUSESTYLE_OUTPUT"), "1")
house_style_v2_run <- house_style_run &&
  identical(Sys.getenv("ROUND48_FIGS45_HOUSESTYLE_V2"), "1")
house_style_v3_run <- house_style_run &&
  identical(Sys.getenv("ROUND48_FIGS45_HOUSESTYLE_V3"), "1")
house_style_v4_run <- house_style_run &&
  identical(Sys.getenv("ROUND48_FIGS45_HOUSESTYLE_V4"), "1")

# Render-time layout measurements are retained in memory for compositors and
# QA messages; no scientific object or persistent analysis output is changed.
round48_tree_layout_metrics <- list()

tree_path <- file.path(project_root, "data", "tree arch", "Tree archaea 1.tree")
asv_path <- file.path(project_root, "data", "arc_unrarefild", "ASV_Arc_200cm_delete_less200.txt")
retained_samples_path <- file.path(project_root, "data", "arc_60cm_202rare", "asv_202.txt")
env_path <- file.path(project_root, "data", "env_2024.csv")
ann_path <- file.path(project_root, "analysis", "tree_audit", "legacy_tree_asv_annotation.csv")
ref_metadata_path <- file.path(project_root, "analysis", "tree_audit", "taxonomic_panels", "legacy_tree_reference_metadata_ncbi.csv")

if (house_style_run) {
  out_dir <- file.path(tempdir(), "chile_archaea_round48_housestyle", "tree_panels")
  itol_dir <- file.path(tempdir(), "chile_archaea_round48_housestyle", "itol")
  figure_dir <- file.path(project_root, "figures", "candidates")
  provisional_dir <- file.path(tempdir(), "chile_archaea_round48_housestyle", "provisional")
} else {
  out_dir <- file.path(project_root, "analysis", "tree_audit", "taxonomic_panels")
  itol_dir <- file.path(project_root, "analysis", "tree_audit", "itol_annotations")
  figure_dir <- file.path(project_root, "figures")
  provisional_dir <- file.path(figure_dir, "provisional")
}

# Reference-tip thinning knob: max informative reference tips kept per nearest
# project ASV in each panel (1 = keep only the single closest reference). Lower =
# sparser, more merging of adjacent same-signature ASVs. See plot_panel().
max_refs_per_asv <- 1L
collapse_adjacent_refs <- TRUE
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(itol_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provisional_dir, recursive = TRUE, showWarnings = FALSE)

site_levels <- c("AZ", "SG", "LC", "NB")
pool_levels <- c("iDNA", "eDNA")
site_pool_levels <- as.vector(t(outer(site_levels, pool_levels, paste, sep = "_")))
site_display <- c(AZ = "AZ", SG = "SG", LC = "LC", NB = "NA")
site_pool_labels <- paste0(site_display[sub("_.*$", "", site_pool_levels)], "_", sub("^.*_", "", site_pool_levels))

site_cols <- c(AZ = "#D95F02", SG = "#1B9E77", LC = "#7570B3", NB = "#66A61E")
class_cols <- c(
  Nitrososphaeria = "#B2182B",
  Thermoplasmata = "#542788",
  Thermoplasmatota_unclassified = "#2166AC",
  Thermoplasmatota_cl = "#E08214",
  Thermoprotei = "#1B7837",
  Bathyarchaeia = "#C51B7D",
  Nanoarchaeia = "#4D4D4D",
  unclassified = "#D6604D",
  Minor_archaea = "#636363",
  Reference = "#BDBDBD"
)
ref_source_cols <- c(
  soil = "#8C510A",
  "root/rhizo" = "#1B9E77",
  thermal = "#D73027",
  "marine sed." = "#2B8CBE",
  marine = "#2B8CBE",
  sediment = "#6A51A3",
  freshwater = "#1F78B4",
  "peat/wetland" = "#7B3294",
  acidic = "#E66101",
  alkaline = "#DFC27D",
  saline = "#80CDC1",
  environmental = "#8C8C8C"
)
phylum_panel_cols <- c(
  Crenarchaeota = "#D73027",
  Thermoplasmatota_plus_Nanoarchaeota = "#7B3294"
)
panel_cols <- c(
  Thermoplasmatota_cl_and_minor = "#E08214",
  Thermoplasmata_plus_unclassified = "#2166AC",
  Thermoprotei_and_minor_Crenarchaeota = "#1B7837",
  Nitrososphaeria = "#B2182B"
)

extract_newick <- function(path) {
  x <- readLines(path, warn = FALSE)
  start <- grep("^\\s*\\(", x)[1]
  if (is.na(start)) stop("No Newick tree found in ", path)
  nwk_lines <- x[start:length(x)]
  end <- grep(";\\s*$", nwk_lines)[1]
  if (is.na(end)) stop("No terminating semicolon found in ", path)
  paste(nwk_lines[seq_len(end)], collapse = "")
}

clean_tax <- function(x) {
  x <- sub("^[a-z]__", "", x)
  x[is.na(x) | x == "" | x == "NA"] <- "unclassified"
  x
}

clean_tax_display <- function(x) {
  x <- clean_tax(x)
  x <- sub("_(cl|or|fa|ge|sp)$", "", x)
  x
}

get_tax_col <- function(df, col) {
  if (col %in% names(df)) return(df[[col]])
  rep(NA_character_, nrow(df))
}

rescale01 <- function(x) {
  if (all(is.na(x)) || diff(range(x, na.rm = TRUE)) == 0) {
    return(rep(0, length(x)))
  }
  (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
}

infer_note <- function(class, order, genus, tip_label, tip_type) {
  text <- paste(class, order, genus, tip_label, collapse = " ")
  if (grepl("Nitroso|Nitro", text, ignore.case = TRUE)) {
    return(c("putative ammonia-oxidizing archaeal lineage (AOA; taxonomy/name-based)", "medium"))
  }
  if (grepl("Methano|Methanomassiliicoccales", text, ignore.case = TRUE)) {
    return(c("putative methanogen lineage (taxonomy/name-based)", "medium"))
  }
  if (grepl("Bathy", text, ignore.case = TRUE)) {
    return(c("Bathyarchaeia; broad anaerobic/sediment-associated lineage, function unresolved here", "low"))
  }
  if (grepl("Woes|Nanoarchae", text, ignore.case = TRUE)) {
    return(c("Nanoarchaeota/DPANN-related lineage; function unresolved from this 16S tree", "low"))
  }
  if (grepl("Thermoplas", text, ignore.case = TRUE)) {
    return(c("Thermoplasmatota-related lineage; function unresolved from this 16S tree", "low"))
  }
  if (grepl("Thermoprotei", text, ignore.case = TRUE)) {
    return(c("Thermoprotei-related lineage; function unresolved from this 16S tree", "low"))
  }
  if (tip_type == "Reference") {
    return(c("ARB/SILVA reference tip; no functional inference from short label alone", "none"))
  }
  c("function unresolved from available taxonomy", "none")
}

infer_guild_code <- function(note) {
  if (grepl("ammonia-oxidizing", note, ignore.case = TRUE)) return("putative AOA")
  if (grepl("methanogen", note, ignore.case = TRUE)) return("putative methanogen")
  if (grepl("Bathyarchaeia", note, ignore.case = TRUE)) return("Bathyarchaeia")
  if (grepl("Nanoarchaeota|DPANN", note, ignore.case = TRUE)) return("DPANN/Nano")
  if (grepl("Thermoplasmatota", note, ignore.case = TRUE)) return("Thermoplasmatota unresolved")
  if (grepl("Thermoprotei", note, ignore.case = TRUE)) return("Thermoprotei unresolved")
  "unresolved"
}

make_ref_source_tag <- function(keyword) {
  if (is.na(keyword) || keyword == "") return("")
  tags <- trimws(strsplit(keyword, ";", fixed = TRUE)[[1]])
  has <- function(z) z %in% tags
  if (has("hot spring/thermal")) return("thermal")
  if (has("root/rhizosphere")) return("root/rhizo")
  if (has("marine") && has("sediment")) return("marine sed.")
  if (has("soil")) return("soil")
  if (has("freshwater")) return("freshwater")
  if (has("peat/wetland")) return("peat/wetland")
  if (has("marine")) return("marine")
  if (has("sediment")) return("sediment")
  if (has("acidic")) return("acidic")
  if (has("alkaline")) return("alkaline")
  if (has("saline")) return("saline")
  if (has("environmental")) return("environmental")
  ""
}

make_ref_context_tag <- function(keyword) {
  if (is.na(keyword) || keyword == "") return("")
  tags <- trimws(strsplit(keyword, ";", fixed = TRUE)[[1]])
  if ("putative AOA context" %in% tags) return("AOA ctx")
  if ("methanogen context" %in% tags) return("methano ctx")
  if ("enrichment" %in% tags) return("enrich.")
  ""
}

tree <- read.tree(text = extract_newick(tree_path))
write.tree(tree, file = file.path(itol_dir, "legacy_tree_clean_newick_for_itol.nwk"))

asv_counts <- read.delim(asv_path, row.names = 1, check.names = FALSE)
retained_samples <- colnames(read.delim(retained_samples_path, row.names = 1, check.names = FALSE))
env <- read.csv(env_path, check.names = FALSE)
ann <- read.csv(ann_path, check.names = FALSE)

sample_cols <- retained_samples[retained_samples %in% colnames(asv_counts) & retained_samples %in% env$sample_id]
if (length(sample_cols) != length(retained_samples)) {
  stop("Expected all 95 retained 0-60 cm samples in the unrarefied table and metadata.")
}
asv_counts <- asv_counts[, sample_cols, drop = FALSE]
env <- env[match(sample_cols, env$sample_id), ]
if (!identical(env$sample_id, sample_cols)) {
  stop("Sample metadata and ASV table columns do not align.")
}
env$site <- factor(env$site, levels = site_levels)
env$dna_type <- factor(env$dna_type, levels = pool_levels)
env$site_pool <- paste(env$site, env$dna_type, sep = "_")

rel_abund <- sweep(asv_counts, 2, colSums(asv_counts), "/") * 100
tree_asvs <- grep("^ASV_[0-9]+$", tree$tip.label, value = TRUE)
tree_refs <- setdiff(tree$tip.label, tree_asvs)

occ_rows <- list()
for (asv in tree_asvs) {
  for (sp in site_pool_levels) {
    idx <- env$site_pool == sp
    raw_values <- as.numeric(asv_counts[asv, idx])
    rel_values <- as.numeric(rel_abund[asv, idx])
    occ_rows[[length(occ_rows) + 1]] <- data.frame(
      ASV = asv,
      site_pool = sp,
      site_pool_display = paste0(site_display[sub("_.*$", "", sp)], "_", sub("^.*_", "", sp)),
      site = sub("_.*$", "", sp),
      site_display = site_display[sub("_.*$", "", sp)],
      dna_type = sub("^.*_", "", sp),
      n_samples = sum(idx),
      n_detected = sum(raw_values > 0, na.rm = TRUE),
      occupancy_fraction = ifelse(sum(idx) > 0, sum(raw_values > 0, na.rm = TRUE) / sum(idx), NA_real_),
      mean_relative_abundance_percent = mean(rel_values, na.rm = TRUE),
      max_relative_abundance_percent = max(rel_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
occ_long <- do.call(rbind, occ_rows)
write.csv(occ_long, file.path(out_dir, "legacy_tree_asv_site_pool_occurrence_long.csv"), row.names = FALSE)
ra_scale <- as.numeric(quantile(
  occ_long$mean_relative_abundance_percent[occ_long$mean_relative_abundance_percent > 0],
  0.95,
  na.rm = TRUE
))
if (!is.finite(ra_scale) || ra_scale <= 0) {
  ra_scale <- max(occ_long$mean_relative_abundance_percent, na.rm = TRUE)
}

occ_wide <- reshape(
  occ_long[, c("ASV", "site_pool", "occupancy_fraction")],
  idvar = "ASV",
  timevar = "site_pool",
  direction = "wide"
)
names(occ_wide) <- sub("^occupancy_fraction\\.", "occ_", names(occ_wide))
ra_wide <- reshape(
  occ_long[, c("ASV", "site_pool", "mean_relative_abundance_percent")],
  idvar = "ASV",
  timevar = "site_pool",
  direction = "wide"
)
names(ra_wide) <- sub("^mean_relative_abundance_percent\\.", "mean_ra_", names(ra_wide))

ann$Class_clean <- clean_tax(get_tax_col(ann, "Class"))
ann$Order_clean <- clean_tax(get_tax_col(ann, "Order"))
ann$Family_clean <- clean_tax(get_tax_col(ann, "Family"))
ann$Genus_clean <- clean_tax(get_tax_col(ann, "Genus"))
ann$Species_clean <- clean_tax(get_tax_col(ann, "Species"))
ann$Class_display <- clean_tax_display(get_tax_col(ann, "Class"))
ann$Order_display <- clean_tax_display(get_tax_col(ann, "Order"))
ann$Family_display <- clean_tax_display(get_tax_col(ann, "Family"))
ann$Genus_display <- clean_tax_display(get_tax_col(ann, "Genus"))
ann$Species_display <- clean_tax_display(get_tax_col(ann, "Species"))
ann$Phylum_clean <- clean_tax(get_tax_col(ann, "Phylum"))
ann$phylum_panel <- ifelse(
  ann$Phylum_clean == "Crenarchaeota",
  "Crenarchaeota",
  "Thermoplasmatota_plus_Nanoarchaeota"
)
ann$taxonomic_panel <- ifelse(
  ann$Class_clean == "Nitrososphaeria", "Nitrososphaeria",
  ifelse(
    ann$Class_clean == "Thermoplasmata", "Thermoplasmata_plus_unclassified",
    ifelse(
      ann$Phylum_clean == "Thermoplasmatota" & ann$Class_clean == "unclassified", "Thermoplasmata_plus_unclassified",
      ifelse(
        ann$Class_clean == "Thermoplasmatota_cl", "Thermoplasmatota_cl_and_minor",
        ifelse(
          ann$Class_clean == "Thermoprotei", "Thermoprotei_and_minor_Crenarchaeota",
          "Thermoprotei_and_minor_Crenarchaeota"
        )
      )
    )
  )
)

ann <- merge(ann, occ_wide, by = "ASV", all.x = TRUE, sort = FALSE)
ann <- merge(ann, ra_wide, by = "ASV", all.x = TRUE, sort = FALSE)
ann$overall_detected_samples <- rowSums(asv_counts[ann$ASV, , drop = FALSE] > 0)
ann$overall_samples <- ncol(asv_counts)
ann$overall_occupancy_fraction <- ann$overall_detected_samples / ann$overall_samples

notes <- t(mapply(
  infer_note,
  ann$Class_clean,
  ann$Order_clean,
  ann$Genus_clean,
  ann$ASV,
  "ASV"
))
ann$inferred_ecological_note <- notes[, 1]
ann$inference_confidence <- notes[, 2]
ann$potential_guild_code <- vapply(ann$inferred_ecological_note, infer_guild_code, character(1))

write.csv(ann, file.path(out_dir, "legacy_tree_asv_occurrence_taxonomy_function_notes.csv"), row.names = FALSE)

tip_meta <- data.frame(label = tree$tip.label, tip_type = ifelse(tree$tip.label %in% tree_asvs, "ASV", "Reference"), stringsAsFactors = FALSE)
tip_meta <- merge(tip_meta, ann, by.x = "label", by.y = "ASV", all.x = TRUE, sort = FALSE)
tip_meta$Phylum_clean[is.na(tip_meta$Phylum_clean)] <- "Reference"
tip_meta$Class_clean[is.na(tip_meta$Class_clean)] <- "Reference"
tip_meta$Family_clean[is.na(tip_meta$Family_clean)] <- "Reference"
tip_meta$Genus_clean[is.na(tip_meta$Genus_clean)] <- "Reference"
tip_meta$Species_clean[is.na(tip_meta$Species_clean)] <- "Reference"
tip_meta$Class_display[is.na(tip_meta$Class_display)] <- "Reference"
tip_meta$Order_display[is.na(tip_meta$Order_display)] <- "Reference"
tip_meta$taxonomic_panel[is.na(tip_meta$taxonomic_panel)] <- "Reference"
tip_meta$phylum_panel[is.na(tip_meta$phylum_panel)] <- "Reference"
tip_meta$display_label <- ifelse(
  tip_meta$tip_type == "ASV",
  paste0(tip_meta$label, "  ", tip_meta$Class_display, "; ", tip_meta$Order_display),
  tip_meta$label
)
tip_meta$ref_source_tag <- ""
tip_meta$ref_context_tag <- ""
tip_meta$ref_accession <- ""
tip_meta$ref_definition <- ""
tip_meta$ref_isolation_source <- ""
tip_meta$ref_geo_loc_name <- ""
if (file.exists(ref_metadata_path)) {
  ref_metadata <- read.csv(ref_metadata_path, check.names = FALSE)
  ref_match <- match(tip_meta$label, ref_metadata$reference_tip)
  ref_keyword <- ref_metadata$reference_keyword[ref_match]
  ref_accession <- ref_metadata$primary_accession[ref_match]
  ref_definition <- ref_metadata$definition[ref_match]
  ref_isolation_source <- ref_metadata$isolation_source[ref_match]
  ref_geo_loc_name <- ref_metadata$geo_loc_name[ref_match]
  has_ref_acc <- tip_meta$tip_type == "Reference" &
    !is.na(ref_accession) &
    ref_accession != "" &
    ref_accession != "NA"
  tip_meta$display_label[has_ref_acc] <- paste0(tip_meta$label[has_ref_acc], " (", ref_accession[has_ref_acc], ")")
  tip_meta$ref_source_tag <- vapply(ref_keyword, make_ref_source_tag, character(1))
  tip_meta$ref_context_tag <- vapply(ref_keyword, make_ref_context_tag, character(1))
  tip_meta$ref_accession <- ifelse(is.na(ref_accession), "", ref_accession)
  tip_meta$ref_definition <- ifelse(is.na(ref_definition), "", ref_definition)
  tip_meta$ref_isolation_source <- ifelse(is.na(ref_isolation_source), "", ref_isolation_source)
  tip_meta$ref_geo_loc_name <- ifelse(is.na(ref_geo_loc_name), "", ref_geo_loc_name)
}

ref_plot_metadata <- tip_meta[tip_meta$tip_type == "Reference", c(
  "label", "display_label", "ref_accession", "ref_source_tag", "ref_context_tag",
  "ref_definition", "ref_isolation_source", "ref_geo_loc_name"
)]
write.csv(ref_plot_metadata, file.path(out_dir, "legacy_tree_reference_metadata_for_plot.csv"), row.names = FALSE)

ref_notes <- tip_meta[tip_meta$tip_type == "Reference", c("label", "tip_type")]
ref_notes$inferred_ecological_note <- vapply(
  ref_notes$label,
  function(z) infer_note(NA, NA, NA, z, "Reference")[1],
  character(1)
)
ref_notes$inference_confidence <- "none"
ref_notes$potential_guild_code <- "reference-only"
write.csv(ref_notes, file.path(out_dir, "legacy_tree_reference_tip_notes.csv"), row.names = FALSE)

tip_dist <- cophenetic.phylo(tree)
nearest_ref_rows <- lapply(tree_refs, function(ref) {
  nearest_asv <- tree_asvs[which.min(tip_dist[ref, tree_asvs])]
  nearest_meta <- ann[match(nearest_asv, ann$ASV), ]
  data.frame(
    reference_tip = ref,
    nearest_asv = nearest_asv,
    distance_to_nearest_asv = tip_dist[ref, nearest_asv],
    phylum_panel = nearest_meta$phylum_panel,
    taxonomic_panel = nearest_meta$taxonomic_panel,
    Phylum_clean = nearest_meta$Phylum_clean,
    Class_clean = nearest_meta$Class_clean,
    stringsAsFactors = FALSE
  )
})
ref_assign <- do.call(rbind, nearest_ref_rows)
write.csv(ref_assign, file.path(out_dir, "legacy_tree_reference_tip_nearest_asv_assignment.csv"), row.names = FALSE)

normalize_tax_value <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x[x %in% c("", "unclassified", "None", "NA")] <- ""
  x
}

short_display_taxonomy <- function(df) {
  class_vals <- unique(as.character(df$Class_display))
  class_vals <- class_vals[!is.na(class_vals)]
  order_vals <- unique(as.character(df$Order_display))
  order_vals <- order_vals[!is.na(order_vals)]
  class_label <- if (length(class_vals) == 1) class_vals else "mixed class"
  order_label <- if (length(order_vals) == 1) order_vals else "mixed order"
  paste(class_label, order_label, sep = "; ")
}

merge_signature <- function(df) {
  sig <- short_display_taxonomy(df)
  ifelse(is.na(sig), "", sig)
}

collapse_group_guild <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return("")
  if (length(x) == 1) return(x)
  "mixed"
}

build_panel_display_rows <- function(tip_xy, panel_name) {
  rows <- list()
  audit_rows <- list()
  idx <- 1
  group_id <- 1

  while (idx <= nrow(tip_xy)) {
    current <- tip_xy[idx, , drop = FALSE]

    if (current$tip_type != "ASV") {
      current$display_n_asvs <- 1L
      current$display_members <- current$label
      current$display_taxonomy_label <- ""
      current$display_group_type <- "reference"
      rows[[length(rows) + 1]] <- current
      idx <- idx + 1L
      next
    }

    end_idx <- idx
    current_signature <- merge_signature(current)
    while (end_idx < nrow(tip_xy) && tip_xy$tip_type[end_idx + 1L] == "ASV") {
      candidate <- tip_xy[idx:(end_idx + 1L), , drop = FALSE]
      candidate_signature <- merge_signature(candidate)
      if (candidate_signature == current_signature && candidate_signature != "") {
        end_idx <- end_idx + 1L
      } else {
        break
      }
    }

    block <- tip_xy[idx:end_idx, , drop = FALSE]
    if (nrow(block) >= 2L && current_signature != "") {
      members <- block$label
      pooled_counts <- asv_counts[members, , drop = FALSE]
      pooled_rel <- rel_abund[members, , drop = FALSE]
      aggregated <- block[1, , drop = FALSE]
      aggregated$y <- mean(block$y)
      aggregated$label <- paste0(panel_name, "_ASV_group_", group_id)
      aggregated$display_n_asvs <- nrow(block)
      aggregated$display_members <- paste(members, collapse = ";")
      aggregated$display_taxonomy_label <- current_signature
      aggregated$display_group_type <- "merged_asv"
      aggregated$display_label <- paste0(nrow(block), " ASVs  ", aggregated$display_taxonomy_label)
      aggregated$overall_detected_samples <- sum(colSums(pooled_counts, na.rm = TRUE) > 0)
      aggregated$overall_samples <- ncol(asv_counts)
      aggregated$overall_occupancy_fraction <- aggregated$overall_detected_samples / aggregated$overall_samples
      aggregated$potential_guild_code <- collapse_group_guild(block$potential_guild_code)
      aggregated$inferred_ecological_note <- collapse_group_guild(block$inferred_ecological_note)

      for (sp in site_pool_levels) {
        sp_idx <- env$site_pool == sp
        pooled_counts_sp <- colSums(pooled_counts[, sp_idx, drop = FALSE], na.rm = TRUE)
        pooled_rel_sp <- colSums(pooled_rel[, sp_idx, drop = FALSE], na.rm = TRUE)
        aggregated[[paste0("occ_", sp)]] <- ifelse(sum(sp_idx) > 0, sum(pooled_counts_sp > 0) / sum(sp_idx), NA_real_)
        aggregated[[paste0("mean_ra_", sp)]] <- mean(pooled_rel_sp, na.rm = TRUE)
      }

      rows[[length(rows) + 1]] <- aggregated
      audit_rows[[length(audit_rows) + 1]] <- data.frame(
        panel_name = panel_name,
        display_label = aggregated$display_label,
        n_asvs = nrow(block),
        members = paste(members, collapse = ";"),
        taxonomy_label = aggregated$display_taxonomy_label,
        stringsAsFactors = FALSE
      )
      group_id <- group_id + 1L
    } else {
      current$display_n_asvs <- 1L
      current$display_members <- current$label
      current$display_taxonomy_label <- short_display_taxonomy(current)
      current$display_group_type <- "single_asv"
      rows[[length(rows) + 1]] <- current
    }

    idx <- end_idx + 1L
  }

  plot_rows <- do.call(rbind, rows)
  audit_df <- if (length(audit_rows) > 0) do.call(rbind, audit_rows) else data.frame()
  list(plot_rows = plot_rows, audit = audit_df)
}

get_rendered_tip_xy <- function(sub_tree) {
  # plot.phylo(plot = FALSE) still consults a graphics device. Keep standalone
  # v2/v3/v4 renders from refreshing the project-root Rplots.pdf side file.
  opened_null_device <- grDevices::dev.cur() == 1L
  if (opened_null_device) {
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  plot.phylo(
    sub_tree,
    type = "phylogram",
    show.tip.label = FALSE,
    no.margin = TRUE,
    direction = "rightwards",
    plot = FALSE
  )
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  data.frame(
    label = sub_tree$tip.label,
    x = pp$xx[seq_along(sub_tree$tip.label)],
    y = pp$yy[seq_along(sub_tree$tip.label)],
    stringsAsFactors = FALSE
  )
}

find_adjacent_reference_drops <- function(sub_tree) {
  tip_xy <- get_rendered_tip_xy(sub_tree)
  tip_xy$tip_type <- tip_meta$tip_type[match(tip_xy$label, tip_meta$label)]
  tip_xy <- tip_xy[order(tip_xy$y, decreasing = TRUE), ]

  drops <- character()
  idx <- 1L
  while (idx <= nrow(tip_xy)) {
    if (tip_xy$tip_type[idx] != "Reference") {
      idx <- idx + 1L
      next
    }

    end_idx <- idx
    while (end_idx < nrow(tip_xy) && tip_xy$tip_type[end_idx + 1L] == "Reference") {
      end_idx <- end_idx + 1L
    }

    run_labels <- tip_xy$label[idx:end_idx]
    if (length(run_labels) >= 2L) {
      ref_idx <- match(run_labels, ref_assign$reference_tip)
      ref_dist <- ref_assign$distance_to_nearest_asv[ref_idx]
      keep_idx <- if (all(is.na(ref_dist))) 1L else which.min(ifelse(is.na(ref_dist), Inf, ref_dist))
      drops <- c(drops, run_labels[-keep_idx])
    }

    idx <- end_idx + 1L
  }

  unique(drops)
}

has_adjacent_references <- function(sub_tree) {
  tip_xy <- get_rendered_tip_xy(sub_tree)
  tip_type <- tip_meta$tip_type[match(tip_xy$label, tip_meta$label)]
  tip_type <- tip_type[order(tip_xy$y, decreasing = TRUE)]
  any(tip_type[-length(tip_type)] == "Reference" & tip_type[-1L] == "Reference")
}

panel_defs <- list(
  Thermoplasmatota_cl_and_minor = list(
    title = "Thermoplasmatota_cl and minor archaeal lineages",
    asv_filter = quote(taxonomic_panel == "Thermoplasmatota_cl_and_minor")
  ),
  Thermoplasmata_plus_unclassified = list(
    title = "Thermoplasmata and unclassified Thermoplasmatota",
    asv_filter = quote(taxonomic_panel == "Thermoplasmata_plus_unclassified")
  ),
  Thermoprotei_and_minor_Crenarchaeota = list(
    title = "Thermoprotei and minor Crenarchaeota/Nanoarchaeota",
    asv_filter = quote(taxonomic_panel == "Thermoprotei_and_minor_Crenarchaeota")
  ),
  Nitrososphaeria = list(
    title = "Nitrososphaeria",
    asv_filter = quote(taxonomic_panel == "Nitrososphaeria")
  )
)

phylum_panel_defs <- list(
  Crenarchaeota = list(
    title = "Crenarchaeota",
    asv_filter = quote(phylum_panel == "Crenarchaeota")
  ),
  Thermoplasmatota_plus_Nanoarchaeota = list(
    title = "Thermoplasmatota and Nanoarchaeota",
    asv_filter = quote(phylum_panel == "Thermoplasmatota_plus_Nanoarchaeota")
  )
)

plot_panel <- function(panel_name, file, defs = panel_defs, title_cols = panel_cols) {
  def <- defs[[panel_name]]
  is_main_tree <- panel_name == "Crenarchaeota" && identical(Sys.getenv("FIG3_MAIN_TREE"), "1")
  use_main_readability_v2 <- is_main_tree &&
    identical(Sys.getenv("FIG2_MAIN_TREE_READABILITY_V2"), "1")
  tighten_track_spacing <- (house_style_v2_run || house_style_v3_run || house_style_v4_run) && !is_main_tree
  use_optical_track_spacing_v3 <- house_style_v3_run && !is_main_tree
  use_optical_source_spacing_v4 <- house_style_v4_run && !is_main_tree
  panel_asvs <- ann$ASV[eval(def$asv_filter, ann)]
  panel_refs <- unique(ref_assign$reference_tip[ref_assign$nearest_asv %in% panel_asvs])
  panel_refs <- panel_refs[!is.na(panel_refs) & panel_refs %in% tree$tip.label]
  panel_ref_meta <- tip_meta[match(panel_refs, tip_meta$label), ]
  informative_ref <- !is.na(panel_ref_meta$ref_source_tag) &
    panel_ref_meta$ref_source_tag != "" &
    panel_ref_meta$ref_source_tag != "environmental"
  panel_refs <- panel_refs[informative_ref]

  # Thin reference tips: keep only the closest `max_refs_per_asv` informative
  # reference(s) per nearest project ASV. This declutters the tree (many ASVs
  # carried 2-6 redundant references) while keeping a representative reference
  # anchor for each lineage -- references are reduced, NOT removed. Pruning the
  # interrupting references lets adjacent same-signature project ASVs become
  # neighbours, so build_panel_display_rows() (run below on the pruned tree)
  # re-merges them and shows the ASV count per merged row.
  if (length(panel_refs) > 0) {
    ref_keep <- ref_assign[ref_assign$reference_tip %in% panel_refs,
                           c("reference_tip", "nearest_asv", "distance_to_nearest_asv")]
    ref_keep <- ref_keep[order(ref_keep$nearest_asv, ref_keep$distance_to_nearest_asv), ]
    ref_keep$rank <- ave(ref_keep$distance_to_nearest_asv, ref_keep$nearest_asv, FUN = seq_along)
    panel_refs <- ref_keep$reference_tip[ref_keep$rank <= max_refs_per_asv]
  }
  message("Panel ", panel_name, ": ", length(panel_asvs), " project ASVs, ",
          length(panel_refs), " reference tips kept (max ", max_refs_per_asv, " per ASV)")

  keep <- unique(c(panel_asvs, panel_refs))
  sub_tree <- keep.tip(tree, keep)

  n_refs_before_collapse <- length(panel_refs)
  refs_dropped_for_adjacency <- character()
  if (collapse_adjacent_refs) {
    refs_dropped_for_adjacency <- find_adjacent_reference_drops(sub_tree)
    if (length(refs_dropped_for_adjacency) > 0) {
      sub_tree <- drop.tip(sub_tree, refs_dropped_for_adjacency)
      panel_refs <- setdiff(panel_refs, refs_dropped_for_adjacency)
    }
    if (has_adjacent_references(sub_tree)) {
      stop("Adjacent reference tips remain after collapse in panel ", panel_name)
    }
  }
  message("Panel ", panel_name, ": adjacent-reference collapse ",
          n_refs_before_collapse, " -> ", length(panel_refs),
          " reference tips; dropped ", length(refs_dropped_for_adjacency),
          "; adjacent refs after collapse: ", has_adjacent_references(sub_tree))

  panel_tip_meta <- tip_meta[match(sub_tree$tip.label, tip_meta$label), ]

  # The Crenarchaeota panel carries more displayed rows than the complementary
  # panel. Give it additional vertical room so labels and bubbles do not read
  # as compressed when the two full trees are viewed side by side.
  panel_height_mult <- if (panel_name == "Crenarchaeota") 0.273 else 0.420
  panel_width <- 19.2
  panel_height <- max(if (panel_name == "Crenarchaeota") 17.2 else 24.0, length(sub_tree$tip.label) * panel_height_mult)
  if (!is.null(file)) {
    grDevices::cairo_pdf(filename = file, width = panel_width, height = panel_height, family = "Helvetica", bg = "white")
    on.exit(dev.off(), add = TRUE)
  }

  par(mar = c(2.2, 1.0, 2.6, 0.7), xpd = NA, family = "Helvetica")
  xlim <- c(0, max(node.depth.edgelength(sub_tree), na.rm = TRUE) + if (panel_name == "Crenarchaeota") 5.90 else 6.00)
  # The revised main-tree layout tightens unused coordinate margins so the
  # fixed journal-height panel gains real row spacing without shrinking text.
  ylim <- if (use_main_readability_v2) {
    c(-18.0, length(sub_tree$tip.label) + 15.0)
  } else {
    c(-27, length(sub_tree$tip.label) + 17.5)
  }
  plot.phylo(
    sub_tree,
    type = "phylogram",
    show.tip.label = FALSE,
    no.margin = TRUE,
    x.lim = xlim,
    y.lim = ylim,
    edge.color = "#9E9E9E",
    edge.width = 0.7,
    direction = "rightwards"
  )

  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_xy <- data.frame(
    label = sub_tree$tip.label,
    x = pp$xx[seq_along(sub_tree$tip.label)],
    y = pp$yy[seq_along(sub_tree$tip.label)],
    stringsAsFactors = FALSE
  )
  tip_xy <- merge(tip_xy, panel_tip_meta, by = "label", all.x = TRUE, sort = FALSE)
  display_build <- build_panel_display_rows(tip_xy, panel_name)
  plot_rows <- display_build$plot_rows
  # The compact main-figure subset reuses this plotting function but must not
  # overwrite the full supplementary-tree merge audit.
  if (nrow(display_build$audit) > 0 && !is_main_tree) {
    write.csv(
      display_build$audit,
      file.path(out_dir, paste0("legacy_tree_panel_", panel_name, "_merged_asv_groups.csv")),
      row.names = FALSE
    )
  }

  node_ids <- (length(sub_tree$tip.label) + 1):(length(sub_tree$tip.label) + sub_tree$Nnode)
  node_labels <- suppressWarnings(as.numeric(sub_tree$node.label))
  node_xy <- data.frame(node = node_ids, x = pp$xx[node_ids], y = pp$yy[node_ids], support = node_labels)
  node_xy <- node_xy[!is.na(node_xy$support), ]
  node_xy$cex <- ifelse(node_xy$support >= 95, 0.65, ifelse(node_xy$support >= 85, 0.52, ifelse(node_xy$support >= 75, 0.40, ifelse(node_xy$support >= 65, 0.32, NA))))
  node_xy$fill <- ifelse(node_xy$support >= 95, "#111111", ifelse(node_xy$support >= 85, "#666666", ifelse(node_xy$support >= 75, "#BDBDBD", "white")))
  node_xy <- node_xy[!is.na(node_xy$cex), ]
  points(node_xy$x, node_xy$y, pch = 21, cex = node_xy$cex, bg = node_xy$fill, col = "#555555", lwd = 0.45)

  x_max <- max(pp$xx, na.rm = TRUE)
  label_x <- x_max + 0.025
  main_column_shift <- if (use_main_readability_v2) 0.28 else 0
  source_gap_default <- if (panel_name == "Crenarchaeota") 1.650 else 2.600
  source_gap <- if (tighten_track_spacing) source_gap_default * 0.50 else source_gap_default
  right_track_shift <- source_gap - source_gap_default
  source_x <- x_max + source_gap + main_column_shift
  occ_x0 <- x_max +
    if (panel_name == "Crenarchaeota") 2.300 else 3.250
  occ_x0 <- occ_x0 + main_column_shift + right_track_shift
  occ_dx <- if (panel_name == "Crenarchaeota") 0.135 else 0.100
  occ_x <- occ_x0 + seq(0, by = occ_dx, length.out = length(site_pool_levels))
  overall_x <- max(occ_x) + if (panel_name == "Crenarchaeota") 0.520 else 0.440
  note_gap_default <- if (panel_name == "Crenarchaeota") 0.850 else 0.650
  note_gap <- if (tighten_track_spacing) note_gap_default / 3 else note_gap_default
  note_x <- overall_x + note_gap
  bubble_base_cex <- if (panel_name == "Crenarchaeota") 1.05 else 0.92
  bubble_gain_cex <- if (panel_name == "Crenarchaeota") 1.45 else 1.32
  bubble_zero_cex <- if (panel_name == "Crenarchaeota") 0.72 else 0.64
  body_text_scale <- if (is_main_tree) 1.45 else if (panel_name == "Crenarchaeota") 1.18 else 1.42
  context_text_scale <- if (is_main_tree) 1.35 else if (panel_name == "Crenarchaeota") 1.12 else 1.36
  header_text_scale <- if (is_main_tree) 1.40 else if (panel_name == "Crenarchaeota") 1.28 else 1.40
  legend_text_scale <- if (is_main_tree) 1.42 else if (panel_name == "Crenarchaeota") 1.36 else 1.42

  asv_idx <- plot_rows$tip_type == "ASV"
  ref_idx <- plot_rows$tip_type == "Reference"
  ref_source_idx <- ref_idx & plot_rows$ref_source_tag != ""
  ref_label_cex <- (if (panel_name == "Crenarchaeota") 0.84 else 0.88) * body_text_scale
  asv_label_cex <- (if (panel_name == "Crenarchaeota") 0.90 else 0.96) * body_text_scale
  ref_tip_right <- label_x + strwidth(
    plot_rows$display_label[ref_idx],
    units = "user",
    cex = ref_label_cex,
    font = 3
  )
  asv_tip_right <- label_x + strwidth(
    plot_rows$display_label[asv_idx],
    units = "user",
    cex = asv_label_cex,
    font = 2
  )
  max_tip_right <- max(c(ref_tip_right, asv_tip_right), na.rm = TRUE)
  max_bold_asv_right <- max(asv_tip_right, na.rm = TRUE)
  source_marker_x <- source_x - 0.055
  source_text_x <- source_x + 0.020

  if (use_optical_track_spacing_v3) {
    # Hold the perceived label-to-column whitespace constant across panels and
    # back-solve each panel's coordinate shift on its own active Cairo device.
    target_gap_in <- 0.35
    inches_per_user_x <- par("pin")[1] / diff(par("usr")[1:2])
    target_gap_user <- target_gap_in / inches_per_user_x
    first_bubble_x_v2 <- occ_x[1]
    source_marker_x_v2 <- source_marker_x
    marker_half_width <- strwidth("M", units = "user", cex = 0.78) / 2

    bubble_shift_required <- max_tip_right + target_gap_user - first_bubble_x_v2
    source_shift_required <- if (any(ref_source_idx)) {
      max(
        label_x + strwidth(
          plot_rows$display_label[ref_source_idx],
          units = "user",
          cex = ref_label_cex,
          font = 3
        ) + target_gap_user + marker_half_width - source_marker_x_v2,
        na.rm = TRUE
      )
    } else {
      -Inf
    }
    optical_block_shift <- max(bubble_shift_required, source_shift_required)

    # Shift the source, all abundance columns, detected count, and final
    # context column as one block. This preserves the v2 within-block geometry,
    # including the accepted one-third detected-to-context gap.
    source_x <- source_x + optical_block_shift
    occ_x <- occ_x + optical_block_shift
    overall_x <- overall_x + optical_block_shift
    note_x <- note_x + optical_block_shift
    source_marker_x <- source_x - 0.055
    source_text_x <- source_x + 0.020

    bold_bubble_clearance_v2 <- first_bubble_x_v2 - max_bold_asv_right
    bold_bubble_clearance_v3 <- occ_x[1] - max_bold_asv_right
    all_label_bubble_clearance_v3 <- occ_x[1] - max_tip_right
    source_label_clearance_v3 <- if (any(ref_source_idx)) {
      min(
        source_marker_x - marker_half_width -
          (label_x + strwidth(
            plot_rows$display_label[ref_source_idx],
            units = "user",
            cex = ref_label_cex,
            font = 3
          )),
        na.rm = TRUE
      )
    } else {
      Inf
    }
    spacing_tolerance_user <- 1e-7 / inches_per_user_x
    if (all_label_bubble_clearance_v3 + spacing_tolerance_user < target_gap_user) {
      stop("Optical label-to-first-bubble gap is below target in panel ", panel_name)
    }
    if (source_label_clearance_v3 + spacing_tolerance_user < target_gap_user) {
      stop("Optical reference-label-to-source-marker gap is below target in panel ", panel_name)
    }

    round48_tree_layout_metrics[[panel_name]] <<- list(
      max_tip_right_user = max_tip_right,
      max_tip_right_device_in = grconvertX(max_tip_right, from = "user", to = "inches"),
      max_bold_asv_right_user = max_bold_asv_right,
      max_bold_asv_right_device_in = grconvertX(max_bold_asv_right, from = "user", to = "inches"),
      target_gap_in = target_gap_in,
      target_gap_user = target_gap_user,
      first_bubble_x_v2_user = first_bubble_x_v2,
      first_bubble_x_v3_user = occ_x[1],
      first_bubble_x_v3_device_in = grconvertX(occ_x[1], from = "user", to = "inches"),
      optical_block_shift_user = optical_block_shift,
      optical_block_shift_in = optical_block_shift * inches_per_user_x,
      bold_bubble_clearance_v2_user = bold_bubble_clearance_v2,
      bold_bubble_clearance_v2_in = bold_bubble_clearance_v2 * inches_per_user_x,
      bold_bubble_clearance_v3_user = bold_bubble_clearance_v3,
      bold_bubble_clearance_v3_in = bold_bubble_clearance_v3 * inches_per_user_x,
      all_label_bubble_clearance_v3_in = all_label_bubble_clearance_v3 * inches_per_user_x,
      source_label_clearance_v3_in = source_label_clearance_v3 * inches_per_user_x,
      source_gap_v2_user = source_gap,
      note_gap_v2_user = note_gap,
      note_gap_v2_in = note_gap * inches_per_user_x
    )
    figure_label <- if (panel_name == "Crenarchaeota") "Fig. S4" else "Fig. S5"
    cat(sprintf(
      paste0(
        "%s v3 optical spacing: max tip right %.6f user / %.6f device in; ",
        "target %.3f in; first bubble %.6f user / %.6f device in; ",
        "v2->v3 block shift %+.6f user / %+.6f in; ",
        "longest bold ASV->bubble %.6f -> %.6f in; ",
        "all-label->bubble %.6f in; reference-label->marker %.6f in; ",
        "detected/context gap retained at %.6f user / %.6f in\n"
      ),
      figure_label,
      max_tip_right,
      grconvertX(max_tip_right, from = "user", to = "inches"),
      target_gap_in,
      occ_x[1],
      grconvertX(occ_x[1], from = "user", to = "inches"),
      optical_block_shift,
      optical_block_shift * inches_per_user_x,
      bold_bubble_clearance_v2 * inches_per_user_x,
      bold_bubble_clearance_v3 * inches_per_user_x,
      all_label_bubble_clearance_v3 * inches_per_user_x,
      source_label_clearance_v3 * inches_per_user_x,
      note_gap,
      note_gap * inches_per_user_x
    ))
  }
  if (use_optical_source_spacing_v4) {
    # The source-marker square, not the first abundance bubble, is the true
    # column immediately adjacent to the tip labels. Measure all rendered tip
    # labels (bold project ASVs and italic references) on this panel's active
    # Cairo device, then translate the complete downstream block as one unit.
    target_gap_in <- 0.20
    pdf_bbox_safety_in <- 0.005
    inches_per_user_x <- par("pin")[1] / diff(par("usr")[1:2])
    target_gap_user <- target_gap_in / inches_per_user_x
    placement_gap_user <- (target_gap_in + pdf_bbox_safety_in) / inches_per_user_x

    # Cairo renders base-graphics pch 15 as a square with side length
    # 0.45 * pointsize * cex points. Convert its half-width to user units so
    # the collision floor is bound to the marker's leftmost ink, not its anchor.
    source_marker_cex <- 0.78
    source_marker_half_width_in <-
      0.5 * 0.45 * par("ps") * source_marker_cex / 72
    source_marker_half_width_user <-
      source_marker_half_width_in / inches_per_user_x
    source_marker_left_v2 <- source_marker_x - source_marker_half_width_user

    source_shift_required <-
      max_tip_right + placement_gap_user - source_marker_left_v2
    optical_block_shift <- max(0, source_shift_required)

    source_x <- source_x + optical_block_shift
    occ_x <- occ_x + optical_block_shift
    overall_x <- overall_x + optical_block_shift
    note_x <- note_x + optical_block_shift
    source_marker_x <- source_x - 0.055
    source_text_x <- source_x + 0.020

    source_marker_left_v4 <-
      source_marker_x - source_marker_half_width_user
    all_tip_right <- c(ref_tip_right, asv_tip_right)
    all_tip_labels <- c(
      plot_rows$display_label[ref_idx],
      plot_rows$display_label[asv_idx]
    )
    max_tip_index <- which.max(all_tip_right)
    max_tip_label <- all_tip_labels[max_tip_index]
    per_tip_source_clearance_v4 <- source_marker_left_v4 - all_tip_right
    min_tip_source_clearance_v4 <- min(per_tip_source_clearance_v4, na.rm = TRUE)
    spacing_tolerance_user <- 1e-7 / inches_per_user_x

    if (min_tip_source_clearance_v4 + spacing_tolerance_user < target_gap_user) {
      stop("Optical all-label-to-source-marker gap is below target in panel ", panel_name)
    }

    round48_tree_layout_metrics[[panel_name]] <<- list(
      max_tip_label = max_tip_label,
      max_tip_right_user = max_tip_right,
      max_tip_right_device_in = grconvertX(max_tip_right, from = "user", to = "inches"),
      target_gap_in = target_gap_in,
      pdf_bbox_safety_in = pdf_bbox_safety_in,
      source_marker_left_v2_user = source_marker_left_v2,
      source_marker_left_v4_user = source_marker_left_v4,
      source_marker_left_v4_device_in = grconvertX(source_marker_left_v4, from = "user", to = "inches"),
      optical_block_shift_user = optical_block_shift,
      optical_block_shift_in = optical_block_shift * inches_per_user_x,
      min_tip_source_clearance_v4_user = min_tip_source_clearance_v4,
      min_tip_source_clearance_v4_in = min_tip_source_clearance_v4 * inches_per_user_x,
      note_gap_v2_user = note_gap,
      note_gap_v2_in = note_gap * inches_per_user_x,
      panel_width_in = dev.size("in")[1]
    )
    figure_label <- if (panel_name == "Crenarchaeota") "Fig. S4" else "Fig. S5"
    cat(sprintf(
      paste0(
        "%s v4 source spacing: longest tip '%s'; max tip right %.6f user / %.6f device in; ",
        "source-marker left %.6f user / %.6f device in; target %.3f in; ",
        "achieved %.6f in; v2->v4 block shift %+.6f user / %+.6f in; ",
        "detected/context gap retained at %.6f user / %.6f in; canvas %.3f in\n"
      ),
      figure_label,
      max_tip_label,
      max_tip_right,
      grconvertX(max_tip_right, from = "user", to = "inches"),
      source_marker_left_v4,
      grconvertX(source_marker_left_v4, from = "user", to = "inches"),
      target_gap_in,
      min_tip_source_clearance_v4 * inches_per_user_x,
      optical_block_shift,
      optical_block_shift * inches_per_user_x,
      note_gap,
      note_gap * inches_per_user_x,
      dev.size("in")[1]
    ))
  }
  text(label_x, plot_rows$y[ref_idx], plot_rows$display_label[ref_idx], adj = c(0, 0.5), cex = (if (panel_name == "Crenarchaeota") 0.84 else 0.88) * body_text_scale, col = "#4D4D4D", font = 3)
  asv_cols <- class_cols[plot_rows$Class_clean[asv_idx]]
  asv_cols[is.na(asv_cols)] <- "#B2182B"
  text(label_x, plot_rows$y[asv_idx], plot_rows$display_label[asv_idx], adj = c(0, 0.5), cex = (if (panel_name == "Crenarchaeota") 0.90 else 0.96) * body_text_scale, col = asv_cols, font = 2)

  ref_source_col <- ref_source_cols[plot_rows$ref_source_tag[ref_source_idx]]
  ref_source_col[is.na(ref_source_col)] <- "#8C8C8C"
  if (tighten_track_spacing && !use_optical_track_spacing_v3 && !use_optical_source_spacing_v4) {
    marker_half_width <- strwidth("M", units = "user", cex = 0.78) / 2
    ref_source_tip_right <- label_x + strwidth(
      plot_rows$display_label[ref_source_idx],
      units = "user",
      cex = ref_label_cex,
      font = 3
    )
    tree_source_clearance <- min(
      source_marker_x - marker_half_width - ref_source_tip_right,
      na.rm = TRUE
    )
    detected_labels <- paste0(
      plot_rows$overall_detected_samples[asv_idx],
      "/",
      plot_rows$overall_samples[asv_idx]
    )
    detected_right <- overall_x + strwidth(
      detected_labels,
      units = "user",
      cex = 0.90 * context_text_scale
    ) / 2
    detected_context_clearance <- min(note_x - detected_right, na.rm = TRUE)
    if (!is.finite(tree_source_clearance) || tree_source_clearance <= 0) {
      stop("Tightened tree-to-reference track overlaps tip text in panel ", panel_name)
    }
    if (!is.finite(detected_context_clearance) || detected_context_clearance <= 0) {
      stop("Tightened detected/context tracks overlap in panel ", panel_name)
    }
    inches_per_user_x <- par("pin")[1] / diff(par("usr")[1:2])
    message(sprintf(
      paste0(
        "Panel %s v2 gaps: tree/ref %.3f -> %.3f user units ",
        "(%.3f -> %.3f in; min text clearance %.3f in); ",
        "detected/context %.3f -> %.3f user units ",
        "(%.3f -> %.3f in; min text clearance %.3f in)"
      ),
      panel_name,
      source_gap_default,
      source_gap,
      source_gap_default * inches_per_user_x,
      source_gap * inches_per_user_x,
      tree_source_clearance * inches_per_user_x,
      note_gap_default,
      note_gap,
      note_gap_default * inches_per_user_x,
      note_gap * inches_per_user_x,
      detected_context_clearance * inches_per_user_x
    ))
  }
  points(
    rep(source_marker_x, sum(ref_source_idx)),
    plot_rows$y[ref_source_idx],
    pch = 15,
    cex = 0.78,
    col = ref_source_col
  )
  text(
    source_text_x,
    plot_rows$y[ref_source_idx],
    plot_rows$ref_source_tag[ref_source_idx],
    adj = c(0, 0.5),
    cex = (if (panel_name == "Crenarchaeota") 0.82 else 0.90) * body_text_scale,
    col = ref_source_col
  )

  ra_cols <- paste0("mean_ra_", site_pool_levels)
  for (j in seq_along(site_pool_levels)) {
    vals <- plot_rows[[ra_cols[j]]]
    vals[is.na(vals)] <- 0
    scaled_vals <- sqrt(pmin(vals / ra_scale, 1))
    site <- sub("_.*$", "", site_pool_levels[j])
    bg <- ifelse(
      vals > 0,
      vapply(scaled_vals, function(v) adjustcolor(site_cols[site], alpha.f = 0.18 + 0.82 * v), character(1)),
      "white"
    )
    cex <- ifelse(vals > 0, bubble_base_cex + bubble_gain_cex * scaled_vals, bubble_zero_cex)
    points(rep(occ_x[j], sum(asv_idx)), plot_rows$y[asv_idx], pch = 21, cex = cex[asv_idx], bg = bg[asv_idx], col = "#8C8C8C", lwd = 0.25)
  }

  text(
    overall_x,
    plot_rows$y[asv_idx],
    paste0(plot_rows$overall_detected_samples[asv_idx], "/", plot_rows$overall_samples[asv_idx]),
    adj = c(0.5, 0.5),
    cex = 0.90 * context_text_scale,
    col = "#404040"
  )
  text(
    note_x,
    plot_rows$y[asv_idx],
    plot_rows$potential_guild_code[asv_idx],
    adj = c(0, 0.5),
    cex = 0.88 * context_text_scale,
    col = "#404040"
  )
  ref_context_idx <- ref_idx & plot_rows$ref_context_tag != ""
  text(
    note_x,
    plot_rows$y[ref_context_idx],
    plot_rows$ref_context_tag[ref_context_idx],
    adj = c(0, 0.5),
    cex = 0.86 * context_text_scale,
    col = "#404040",
    font = 3
  )

  title_col <- title_cols[[panel_name]]
  if (is.null(title_col)) title_col <- "#333333"
  n_display_groups <- sum(asv_idx)
  title_y <- max(plot_rows$y) + 12.4
  title_cex <- 1.40
  title_height_user <- strheight(def$title, units = "user", cex = title_cex, font = 2)
  title_top_user <- title_y + title_height_user / 2
  panel_metrics <- round48_tree_layout_metrics[[panel_name]]
  if (is.null(panel_metrics)) panel_metrics <- list()
  panel_metrics$title_center_user <- title_y
  panel_metrics$title_height_user <- title_height_user
  panel_metrics$title_top_user <- title_top_user
  panel_metrics$title_top_ndc <- grconvertY(title_top_user, from = "user", to = "ndc")
  panel_metrics$title_top_device_in <- grconvertY(title_top_user, from = "user", to = "inches")
  round48_tree_layout_metrics[[panel_name]] <<- panel_metrics
  text(0, title_y, def$title, adj = c(0, 0.5), cex = title_cex, font = 2, col = title_col)
  text(0, max(plot_rows$y) + 10.2, paste0(length(panel_asvs), " ASVs shown as ", n_display_groups, " display rows; ", length(panel_refs), " nearest reference tips retained"), adj = c(0, 0.5), cex = 0.84 * header_text_scale, col = "#404040")
  text(label_x, max(plot_rows$y) + 8.1, "Tip label", adj = c(0, 0.5), cex = 0.90 * header_text_scale, font = 2)
  text(source_x, max(plot_rows$y) + 8.1, "Ref.\nsource", adj = c(0, 0.5), cex = 0.95 * header_text_scale, font = 2)
  text(mean(range(occ_x)), max(plot_rows$y) + 8.1, "ASV mean relative abundance\nby site x DNA pool", cex = 0.90 * header_text_scale, font = 2)
  text(overall_x, max(plot_rows$y) + 8.1, "detected\nsamples", cex = 0.90 * header_text_scale, font = 2)
  text(note_x, max(plot_rows$y) + 8.1, "ASV guild /\nref. context", adj = c(0, 0.5), cex = 0.88 * header_text_scale, font = 2)
  text(occ_x, rep(max(plot_rows$y) + 1.4, length(occ_x)), site_pool_labels, srt = 55, cex = 0.86 * header_text_scale, adj = c(0, 0.5))

  legend_y <- min(plot_rows$y) - 5.8
  legend_x <- 0
  text(legend_x, legend_y, "Abundance circles: mean relative abundance by site x DNA pool; merged ASV rows sum member ASVs before averaging.", adj = c(0, 0.5), cex = 0.78 * legend_text_scale)
  legend_block_y <- legend_y - 2.1
  site_legend_x <- legend_x + 0.12
  site_legend_y <- legend_block_y
  site_legend_gap <- 0.92
  points(rep(site_legend_x, 4), site_legend_y - c(0, site_legend_gap, 2 * site_legend_gap, 3 * site_legend_gap), pch = 21, cex = 1.25, bg = site_cols[site_levels], col = "#8C8C8C", lwd = 0.25)
  text(rep(site_legend_x + 0.04, 4), site_legend_y - c(0, site_legend_gap, 2 * site_legend_gap, 3 * site_legend_gap), site_display[site_levels], adj = c(0, 0.5), cex = 0.84 * legend_text_scale)
  size_legend_x <- if (panel_name == "Crenarchaeota") note_x + 0.10 else legend_x + 0.74
  size_legend_y <- if (panel_name == "Crenarchaeota") legend_block_y - 0.65 else legend_block_y - 0.20
  text(size_legend_x, size_legend_y + 1.05, "Rel. abund. (%)", adj = c(0, 0.5), cex = 0.82 * legend_text_scale, font = 2, col = "#404040")
  legend_demo_x <- rep(size_legend_x, 3)
  legend_demo_y <- size_legend_y - c(0.05, 1.10, 2.25)
  legend_demo_vals <- c(0.1, 1, 5)
  legend_demo_scaled <- sqrt(pmin(legend_demo_vals / ra_scale, 1))
  points(
    legend_demo_x,
    legend_demo_y,
    pch = 21,
    cex = bubble_base_cex + bubble_gain_cex * legend_demo_scaled,
    bg = vapply(legend_demo_scaled, function(v) adjustcolor("#4D4D4D", alpha.f = 0.18 + 0.82 * v), character(1)),
    col = "#8C8C8C",
    lwd = 0.25
  )
  text(legend_demo_x + 0.135, legend_demo_y, paste0(legend_demo_vals, "%"), adj = c(0, 0.5), cex = 0.80 * legend_text_scale, col = "#404040")

  ref_legend_tags <- names(ref_source_cols)[names(ref_source_cols) %in% unique(plot_rows$ref_source_tag[plot_rows$tip_type == "Reference"])]
  if (length(ref_legend_tags) > 0) {
    ref_legend_tags <- ref_legend_tags[seq_len(min(length(ref_legend_tags), 8))]
    ref_legend_x <- legend_x + if (panel_name == "Crenarchaeota") 1.82 else 2.15
    ref_legend_y <- legend_block_y - 0.40
    text(ref_legend_x, ref_legend_y + 0.92, "Reference source tags", adj = c(0, 0.5), cex = 0.90 * legend_text_scale, font = 2)
    points(
      rep(ref_legend_x, length(ref_legend_tags)),
      ref_legend_y - seq(0, by = 0.82, length.out = length(ref_legend_tags)),
      pch = 15,
      cex = 0.85,
      col = ref_source_cols[ref_legend_tags]
    )
    text(
      rep(ref_legend_x + 0.065, length(ref_legend_tags)),
      ref_legend_y - seq(0, by = 0.82, length.out = length(ref_legend_tags)),
      ref_legend_tags,
      adj = c(0, 0.5),
      cex = 0.80 * legend_text_scale,
      col = ref_source_cols[ref_legend_tags]
    )
  }

  scale_len <- 0.10
  scale_y <- min(plot_rows$y) - 15.0
  segments(0, scale_y, scale_len, scale_y, lwd = 1.2)
  segments(0, scale_y - 0.35, 0, scale_y + 0.35, lwd = 1.0)
  segments(scale_len, scale_y - 0.35, scale_len, scale_y + 0.35, lwd = 1.0)
  text(0, scale_y - 1.1, "0.10 substitutions/site", adj = c(0, 0.5), cex = 0.85 * legend_text_scale)
}

plot_class_panels <- FALSE
write_itol_files <- FALSE

if (plot_class_panels) {
  for (panel_name in names(panel_defs)) {
    plot_panel(panel_name, file.path(out_dir, paste0("legacy_tree_panel_", panel_name, "_refmeta.pdf")))
  }
}

if (!identical(Sys.getenv("FIG3_TREE_SOURCE_ONLY"), "1")) {
  for (panel_name in names(phylum_panel_defs)) {
    panel_file <- file.path(out_dir, paste0("legacy_tree_phylum_panel_", panel_name, "_refmeta.pdf"))
    plot_panel(
      panel_name,
      panel_file,
      phylum_panel_defs,
      phylum_panel_cols
    )
    if (panel_name == "Crenarchaeota") {
      if (house_style_run) {
        target_file <- if (house_style_v4_run) {
          file.path(figure_dir, "FigS4_round48_housestyle_v4.pdf")
        } else if (house_style_v3_run) {
          file.path(figure_dir, "FigS4_round48_housestyle_v3.pdf")
        } else if (house_style_v2_run) {
          file.path(figure_dir, "FigS4_round48_housestyle_v2.pdf")
        } else {
          file.path(figure_dir, "FigS4_round48_housestyle.pdf")
        }
        file.copy(panel_file, target_file, overwrite = TRUE)
      } else {
        file.copy(panel_file, file.path(figure_dir, "14_FigS14_phylogenetic_context_Crenarchaeota_AOA.pdf"), overwrite = TRUE)
        file.copy(panel_file, file.path(figure_dir, "FigS4_Crenarchaeota_AOA_phylogenetic_context.pdf"), overwrite = TRUE)
      }
    } else if (panel_name == "Thermoplasmatota_plus_Nanoarchaeota") {
      if (house_style_run) {
        target_file <- if (house_style_v4_run) {
          file.path(figure_dir, "FigS5_round48_housestyle_v4.pdf")
        } else if (house_style_v3_run) {
          file.path(figure_dir, "FigS5_round48_housestyle_v3.pdf")
        } else if (house_style_v2_run) {
          file.path(figure_dir, "FigS5_round48_housestyle_v2.pdf")
        } else {
          file.path(figure_dir, "FigS5_round48_housestyle.pdf")
        }
        file.copy(panel_file, target_file, overwrite = TRUE)
      } else {
        file.copy(panel_file, file.path(figure_dir, "FigS5_Thermoplasmatota_Nanoarchaeota_phylogenetic_context.pdf"), overwrite = TRUE)
      }
    }
  }
}

write_itol_heatmap <- function(df, value_prefix, out_file, label) {
  cols <- paste0(value_prefix, site_pool_levels)
  lines <- c(
    "DATASET_HEATMAP",
    "SEPARATOR TAB",
    paste0("DATASET_LABEL\t", label),
    "COLOR\t#4575B4",
    paste0("FIELD_LABELS\t", paste(site_pool_labels, collapse = "\t")),
    paste0("FIELD_COLORS\t", paste(rep("#999999", length(site_pool_levels)), collapse = "\t")),
    "SHOW_INTERNAL\t0",
    "DATA"
  )
  dat <- df[, c("ASV", cols)]
  dat[is.na(dat)] <- 0
  dat_lines <- apply(dat, 1, function(z) paste(z, collapse = "\t"))
  writeLines(c(lines, dat_lines), out_file)
}

write_itol_colorstrip <- function(df, out_file) {
  panel_colors <- c(
    Nitrososphaeria = "#B2182B",
    Thermoplasmata_plus_unclassified = "#2166AC",
    Thermoplasmatota_cl_and_minor = "#E08214",
    Thermoprotei_and_minor_Crenarchaeota = "#1B7837"
  )
  lines <- c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tTaxonomic panel",
    "COLOR\t#000000",
    "DATA"
  )
  dat_lines <- paste(df$ASV, panel_colors[df$taxonomic_panel], df$taxonomic_panel, sep = "\t")
  writeLines(c(lines, dat_lines), out_file)
}

write_itol_text_note <- function(df, out_file) {
  lines <- c(
    "DATASET_TEXT",
    "SEPARATOR TAB",
    "DATASET_LABEL\tEcological note",
    "COLOR\t#000000",
    "DATA"
  )
  dat_lines <- paste(df$ASV, df$inferred_ecological_note, sep = "\t")
  writeLines(c(lines, dat_lines), out_file)
}

write_itol_guild_note <- function(df, out_file) {
  lines <- c(
    "DATASET_TEXT",
    "SEPARATOR TAB",
    "DATASET_LABEL\tPotential guild",
    "COLOR\t#000000",
    "DATA"
  )
  dat_lines <- paste(df$ASV, df$potential_guild_code, sep = "\t")
  writeLines(c(lines, dat_lines), out_file)
}

if (write_itol_files) {
  write_itol_heatmap(ann, "occ_", file.path(itol_dir, "itol_heatmap_site_pool_occupancy_fraction.txt"), "Site-pool occupancy")
  write_itol_heatmap(ann, "mean_ra_", file.path(itol_dir, "itol_heatmap_site_pool_mean_relative_abundance.txt"), "Site-pool mean relative abundance")
  write_itol_colorstrip(ann, file.path(itol_dir, "itol_colorstrip_taxonomic_panel.txt"))
  write_itol_text_note(ann, file.path(itol_dir, "itol_text_ecological_note.txt"))
  write_itol_guild_note(ann, file.path(itol_dir, "itol_text_potential_guild.txt"))
}

message("Reference-metadata phylum panels written to: ", out_dir)
if (write_itol_files) message("iTOL annotation files written to: ", itol_dir)
