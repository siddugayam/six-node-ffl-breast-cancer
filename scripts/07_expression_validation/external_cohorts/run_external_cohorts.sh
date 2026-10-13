#!/bin/bash
## Full external-cohort replication pipeline. Runs from the repository root.
set -e
cd /path/to/six-node-ffl-breast-cancer
R=/usr/bin/Rscript
$R scripts/07_expression_validation/external_cohorts/N1_build_cohorts.R   > logs/expression_validation/N1_build.log    2>&1; echo "N1 ok"
$R scripts/07_expression_validation/external_cohorts/N1b_build_geo.R      > logs/expression_validation/N1b_geo.log     2>&1; echo "N1b ok"
$R scripts/07_expression_validation/external_cohorts/N1c_build_geo_gse20685_gse22219.R     > logs/expression_validation/N1c_geo2.log    2>&1; echo "N1c ok"
if [ -f cache/external_cohorts/GSE96058_subset.csv ]; then
  $R scripts/07_expression_validation/external_cohorts/N1e_build_scanb.R  > logs/expression_validation/N1e_scanb.log   2>&1; echo "N1e ok"
else echo "N1e SKIPPED (no SCAN-B subset)"; fi
$R scripts/07_expression_validation/external_cohorts/N2_cohort_analysis.R > logs/expression_validation/N2_analysis.log 2>&1; echo "N2 ok"
$R scripts/07_expression_validation/external_cohorts/N3_meta.R            > logs/expression_validation/N3_meta.log     2>&1; echo "N3 ok"
$R scripts/07_expression_validation/external_cohorts/N4_normal_vs_tumour.R> logs/expression_validation/N4_contrast.log 2>&1; echo "N4 ok"
