# Chilean soil archaea iDNA/eDNA manuscript: Round03 reproducibility release

This repository contains the processed data and current R code supporting the
manuscript **“Archaeal communities in intracellular and extracellular DNA
fractions share broad regional patterns but diverge locally in Chilean
soils.”**

The analysis compares operational intracellular-DNA (iDNA) and
extracellular-DNA (eDNA) fractions across four Chilean sites and five soil
depths. DNA fraction is not treated as a direct live/dead classification, and
archaeal sequence proportion is not interpreted as absolute abundance,
activity, viability, or DNA mass.

Raw demultiplexed 16S rRNA gene sequences are available from the European
Nucleotide Archive under [PRJEB73502](https://www.ebi.ac.uk/ena/browser/view/PRJEB73502).
The repository provides the processed tables needed by the scripts below.

## Start here

Run commands from the repository root. The main analysis backbone is:

```bash
Rscript scripts/13_submission_statistics.R
Rscript scripts/97_round48_paired_concordance_and_control_sensitivity.R
```

The first script creates the core diversity and community summaries. The
second reconstructs the 39 complete iDNA/eDNA horizon pairs, paired-community
concordance, restricted fraction test, and conservative negative-control
sensitivity analyses. Other analysis scripts are independent audits or use
the same processed inputs.

## Current figure map

| Manuscript label | Scientific role | Script(s) |
| --- | --- | --- |
| Fig. 1 | Site context and archaeal sequence proportion | `100_Fig1_frame_aligned_v10.R` |
| Fig. 2 | Class composition and Crenarchaeota phylogeny | `87_Fig2_tree_readability_round45.R`, sourcing `24_...R` and `22_...R` |
| Fig. 3 | Paired concordance-divergence evidence chain | `101_Fig3_concordance_divergence_round50.R`, after `97_...R` |
| Fig. 4 | dbRDA and hierarchical partitioning | `30_Fig4_dbRDA_plus_HP.R` |
| Fig. S1 | Phylum/genus composition | `33_SI_taxonomic_composition_phylum_genus.R` |
| Fig. S2 | Rarefaction and 500-read sensitivity | `32_SI_QC_rarefaction_plus_500sensitivity.R` |
| Figs. S3-S4 | Complete phylogenetic tracks | `22_legacy_tree_refmeta_phylum_panels.R` |
| Fig. S5 | Archaeal diversity and descriptive cross-domain trajectories | `95_Fig3_crossdomain_diversity_round47.R`, after `13_...R` |
| Fig. S6 | iDNA/eDNA ASV overlap | `68_FigS3_iDNA_eDNA_ASV_overlap.R` |
| Fig. S7 | Profile-aware paired sequence-proportion offsets | `104_round55_profile_read_share_offsets.R`, after `97_...R` |

Several script filenames retain historical figure numbers. The manuscript
labels in this table are authoritative.

To reproduce the frozen house-style variants used in the manuscript:

```bash
ROUND48_HOUSESTYLE_OUTPUT=1 ROUND48_FIG2_HOUSESTYLE_V3=1 \
  Rscript scripts/87_Fig2_tree_readability_round45.R
ROUND48_HOUSESTYLE_OUTPUT=1 ROUND48_FIG4_HOUSESTYLE_V2=1 \
  Rscript scripts/30_Fig4_dbRDA_plus_HP.R
ROUND48_HOUSESTYLE_OUTPUT=1 Rscript scripts/33_SI_taxonomic_composition_phylum_genus.R
ROUND48_HOUSESTYLE_OUTPUT=1 Rscript scripts/32_SI_QC_rarefaction_plus_500sensitivity.R
ROUND48_HOUSESTYLE_OUTPUT=1 ROUND48_FIGS45_HOUSESTYLE_V4=1 \
  Rscript scripts/22_legacy_tree_refmeta_phylum_panels.R
ROUND48_HOUSESTYLE_OUTPUT=1 Rscript scripts/95_Fig3_crossdomain_diversity_round47.R
ROUND48_HOUSESTYLE_OUTPUT=1 ROUND48_FIGS3_HOUSESTYLE_V2=1 \
  Rscript scripts/68_FigS3_iDNA_eDNA_ASV_overlap.R
```

Fig. 1 retrieves Natural Earth terrain assets when they are not already
cached, so its first run requires network access. The statistical analyses and
the other figures use repository data only.

## Analysis and sensitivity scripts

| Script | Role |
| --- | --- |
| `13_submission_statistics.R` | Alpha diversity, composition, PERMANOVA, and dbRDA summaries |
| `36_reviewer2_transparency_tables.R` | Coverage and complete-case retention audit |
| `44_restricted_permutation_audit.R` | Pair-restricted permutation sensitivity |
| `61_profile_level_indval_audit.R` | Profile-level indicator-ASV analysis |
| `62_design_aware_site_depth_audit.R` | Site and within-profile depth sensitivities |
| `65_unrarefied_composition_and_aitchison_audit.R` | Unrarefied composition and robust-Aitchison sensitivity |
| `79_round40_reporting_audit.R` | Dispersion, rare-ASV, and descriptive reporting checks |
| `81_round41_negative_control_transparency.R` | Negative-control and raw-to-pre-rarefaction provenance audit |
| `97_round48_paired_concordance_and_control_sensitivity.R` | Current paired analysis and control-filter sensitivity backbone |
| `161_round92_LC_control_filter_class_sensitivity.R` | La Campana class-composition sensitivity |

## Repository layout

```text
scripts/                         current analysis and figure scripts
data/                            processed count, taxonomy, metadata, and tree inputs
analysis/reframed_figures/       fixed Fig. 1 source data
analysis/submission_statistics/  fixed alpha-diversity source for direct Fig. S5 use
analysis/tree_audit/             ASV and reference-tip annotations
analysis/reproducibility/        recorded R session information
figures/                         generated outputs (created by scripts)
```

The canonical rarefied table is
`data/arc_60cm_202rare/asv_202.txt`. The same file at repository root is a
historical compatibility duplicate. `scripts/00_config.R` and
`scripts/15_iDNA_eDNA_paired_composition_similarity.R`, if present in the
repository history, are not part of the Round03 analysis map above.

## Design and interpretation contract

- Independent soil profiles/pits are the biological replicates.
- Paired fraction comparisons use the explicit `site × pit × depth` key.
- The retained archaeal dataset contains 95 fraction libraries; 39 horizons
  have one iDNA and one eDNA library.
- The profile-aware sequence-proportion analysis uses 12 profile means for
  inference; horizon pairs remain descriptive display units.
- Environmental predictors in the dbRDA/hierarchical-partitioning analysis
  are correlated. Contributions are model-dependent, not isolated causal
  effects.
- The control filter is a deliberately conservative sensitivity analysis, not
  a contaminant classification and not a replacement primary pipeline.

## Software

The frozen release was parsed with R 4.5.1. Required packages across the
scripts include `ape`, `cowplot`, `dplyr`, `ggplot2`, `ggrepel`, `ggtree`,
`grid`, `gridGraphics`, `gtable`, `indicspecies`, `patchwork`, `permute`,
`purrr`, `rdacca.hp`, `readr`, `rnaturalearth`, `rnaturalearthdata`, `scales`,
`sf`, `stringr`, `terra`, `tibble`, `tidyr`, `tidyverse`, and `vegan`.
Recorded session information is in
`analysis/reproducibility/round48_R_session_info.txt`.

## Public-release notes

`data/file_A.txt` is a scope-reduced table constructed for this archaeal
manuscript. It preserves all 1,166 archaeal ASV rows and their counts exactly.
The 59,488 non-archaeal ASV rows are represented by one
`NON_ARCHAEAL_TOTAL` row, preserving the exact total prokaryotic read count in
each library. This retains every quantity used by Fig. 1, rarefaction, and the
negative-control sensitivity scripts, but it cannot be used to analyse
non-archaeal taxonomic composition. See `DATA_NOTES.md` for validation details.

Personal filesystem paths, workstation logs, and the original person-named
tree filename were removed from the public release. Scientific sample IDs,
site labels, and site coordinates remain because they are analysis variables,
not personal identifiers. File-level SHA-256 hashes are listed in
`REPRODUCIBILITY_MANIFEST.tsv`.

No software or data reuse license has yet been declared in this repository.

