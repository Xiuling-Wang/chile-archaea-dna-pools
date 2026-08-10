# Processed data notes

## Main archaeal inputs

- `data/arc_60cm_202rare/asv_202.txt`: 202-read rarefied archaeal ASV table
  used for the primary community analyses.
- `data/arc_unrarefild/ASV_Arc_200cm_delete_less200.txt`: pre-rarefaction
  archaeal count table.
- `data/arc_100_202rarefild/tax_202.txt`: SILVA v138 taxonomy used by the
  retained archaeal ASVs.
- `data/env_2024.csv`: fraction-library metadata and environmental variables.
- `data/tree arch/Tree archaea 1.tree`: Newick tree used by the phylogenetic
  panels. RAxML/ARB console logs and personal filesystem paths were removed.

## Scope-reduced total-library table

The archived source `file_A.txt` contained 60,654 ASV rows from all prokaryotic
libraries. The public `data/file_A.txt` contains:

- 1,166 unchanged archaeal ASV rows;
- one `NON_ARCHAEAL_TOTAL` row containing the per-library sum of 59,488
  non-archaeal ASV rows;
- the original 320 biological-library columns, 27 negative/no-template
  control columns, and six positive-control columns.

The release builder verified that every per-library column sum is identical
before and after aggregation. The archaeal rows used for taxonomy,
rarefaction, and control sensitivity are unchanged. This table therefore
reproduces the manuscript's archaeal sequence proportions and control audits,
while deliberately withholding taxon-level bacterial/non-archaeal counts that
are outside this manuscript's scope.

## Companion bacterial inputs

`data/companion_bacteria/processed/bac_200_9435rarefild/` contains only the two
tables used for the descriptive, sample-matched cross-domain trajectories in
Fig. S5. They are not used for the archaeal hypothesis tests.

## Provenance

Raw sequence data are archived at ENA under PRJEB73502. SHA-256 hashes for
every public-release file are listed in `REPRODUCIBILITY_MANIFEST.tsv`.

