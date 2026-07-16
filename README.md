# Reproducibility code — Chile soil archaea iDNA/eDNA manuscript

These R scripts reproduce the figures and figure-underlying statistics of the manuscript. Raw sequences are not included here. The demultiplexed 16S rRNA gene sequences are available in the European Nucleotide Archive (ENA; PRJEB73502; https://www.ebi.ac.uk/ena/browser/view/PRJEB73502). This accession contains the same universal-primer libraries analysed in the companion bacterial study [32]. The rarefied (202 reads per sample) processed archaeal ASV count table used by these scripts is included in this repository as `asv_202.txt`. Fraction-library metadata and soil physicochemical variables are provided in Supplementary Table S1; alpha-diversity, indicator-ASV, coverage, and dbRDA-retention summaries are provided in Supplementary Tables S2 and S5-S8. The ASV taxonomy and sample-metadata tables will be provided through the same repository at submission.

## Figure -> script

| Figure | Script(s) |
| --- | --- |
| Fig 1 | `scripts/100_Fig1_frame_aligned_v10.R` |
| Fig 2 | `scripts/24_Fig3_reframed_phylum_phylo_context.R` (panel A) + `scripts/22_legacy_tree_refmeta_phylum_panels.R` (phylogeny tracks) + `scripts/87_Fig2_tree_readability_round45.R` (compose/house-style) |
| Fig 3 | `scripts/101_Fig3_concordance_divergence_round50.R` (analysis backbone: `scripts/97_round48_paired_concordance_and_control_sensitivity.R`, `scripts/15_iDNA_eDNA_paired_composition_similarity.R`) |
| Fig 4 | `scripts/30_Fig4_dbRDA_plus_HP.R` |
| Fig S1 | `scripts/32_SI_QC_rarefaction_plus_500sensitivity.R` |
| Fig S2 | `scripts/33_SI_taxonomic_composition_phylum_genus.R` |
| Fig S3 | `scripts/68_FigS3_iDNA_eDNA_ASV_overlap.R` |
| Fig S4 | `scripts/22_legacy_tree_refmeta_phylum_panels.R` |
| Fig S5 | `scripts/22_legacy_tree_refmeta_phylum_panels.R` |
| Fig S6 | `scripts/95_Fig3_crossdomain_diversity_round47.R` |
| Fig S7 | `scripts/104_round55_profile_read_share_offsets.R` |
| shared | `scripts/00_config.R` |

Some historical filenames (e.g. `24_Fig3_...`, `95_Fig3_...`) predate the current figure numbering; the table above is authoritative.

Run the scripts from the repository root. Expected inputs belong under `data/`; generated analysis and figure outputs are written under `analysis/` and `figures/`. The processed ASV count table (`asv_202.txt`) is provided at the repository root; other tables are as described above.

Requirements: R packages `ape`, `cowplot`, `dplyr`, `ggplot2`, `ggrepel`, `ggtree`, `grDevices`, `grid`, `gridGraphics`, `gtable`, `patchwork`, `permute`, `purrr`, `rdacca.hp`, `readr`, `rnaturalearth`, `rnaturalearthdata`, `scales`, `sf`, `stringr`, `terra`, `tibble`, `tidyr`, `tidyverse`, `tools`, and `vegan`.
