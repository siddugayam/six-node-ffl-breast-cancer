#!/bin/bash
## Full new-cohort replication pipeline (v2). Run from the revision directory.
set -e
cd /path/to/revision
R=/usr/bin/Rscript
$R scripts/07_expression_validation/external_cohorts/N1_build_cohorts.R   > logs/v2/N1_build.log    2>&1; echo "N1 ok"
$R scripts/07_expression_validation/external_cohorts/N1b_build_geo.R      > logs/v2/N1b_geo.log     2>&1; echo "N1b ok"
$R scripts/07_expression_validation/external_cohorts/N1c_build_geo2.R     > logs/v2/N1c_geo2.log    2>&1; echo "N1c ok"
if [ -f cache/newcohorts/GSE96058_subset.csv ]; then
  $R scripts/07_expression_validation/external_cohorts/N1e_build_scanb.R  > logs/v2/N1e_scanb.log   2>&1; echo "N1e ok"
else echo "N1e SKIPPED (no SCAN-B subset)"; fi
$R scripts/07_expression_validation/external_cohorts/N2_cohort_analysis.R > logs/v2/N2_analysis.log 2>&1; echo "N2 ok"
$R scripts/07_expression_validation/external_cohorts/N3_meta.R            > logs/v2/N3_meta.log     2>&1; echo "N3 ok"
$R scripts/07_expression_validation/external_cohorts/N4_normal_vs_tumour.R> logs/v2/N4_contrast.log 2>&1; echo "N4 ok"
