# Script map

Run scripts from the repository root. Historical numbers in filenames record
development chronology and do not define the current figure order.

Recommended dependency order:

1. `13_submission_statistics.R`
2. `97_round48_paired_concordance_and_control_sensitivity.R`
3. `101_Fig3_concordance_divergence_round50.R`
4. `104_round55_profile_read_share_offsets.R`

Fig. S5 can be run directly because its frozen archaeal alpha-diversity input
is included; rerunning `13_submission_statistics.R` recreates that table.
Other figure and audit scripts can be run independently from the processed
inputs documented in the top-level README.

The Round03 public set contains 12 figure scripts and nine additional analysis
or sensitivity scripts. Local manuscript-packaging utilities are not included
because they operate on author files rather than scientific data.
