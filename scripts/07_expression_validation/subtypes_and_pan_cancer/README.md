# Subtypes and pan-cancer generalisation

Two versions of these analyses ran. Scripts `30`–`36` are the second version and wrote the deposited tables
`results/v2/subtype_stratified.csv`, `pancancer_axis.csv` and `cptac_pancancer_mir29.csv`; scripts `20`–`23` are the
first version, and their other outputs are deposited beside them. The pan-cancer values of the paper come from
`scripts/10_survival_and_clinical/mirna_survival_meta_analysis/c16_pancancer_axes.R`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `20_pam50_subtype.R` | PAM50 subtype stratification of the miR-29 -> collagen axis and the TF -> collagen axis, in TCGA-BRCA; replication in METABRIC. |
| `21_pancancer_axis.py` | Pan-cancer generalisation of the miR-29 -> collagen axis and of the TF -> collagen (stromal-confounded) axis, across all TCGA cohorts. |
| `22_cptac_pancancer.R` | CPTAC PAN-CANCER PROTEOME Replicate the miR-29 -> collagen PROTEIN-level result in every CPTAC tumour type that has both a mature-miRNA matrix and a harmonised proteome. |
| `23_consolidate_generalisation.py` | Consolidation of the subtype, pan-cancer and CPTAC pan-cancer analyses. |
| `30_pam50_subtype.R` | PAM50 subtype stratification of the miR-29 and TF → collagen axes (second implementation, which wrote the deposited tables). |
| `31_pancancer_axis.py` | Pan-cancer generalisation of the miR-29 -> collagen and TF -> collagen axes. |
| `32_cptac_pancancer.R` | CPTAC pan-cancer: miR-29 -> collagen at PROTEIN vs mRNA level. |
| `33_subtype_contrasts.R` | Pairwise subtype contrasts + within-subtype pooled estimates (TCGA-BRCA, METABRIC) |
| `34_consolidate_generalisation.py` | Consolidation / reporting for Parts A-C (fresh run 2026-09-09). |
| `35_pancancer_ecological_product.py` | Ecological decomposition: per TCGA cohort compute rho(TF, CAF), rho(CAF, COL1A1), their product, and the observed rho(TF, COL1A1). |
| `36_meta_analysis.py` | Random-effects (DerSimonian-Laird) meta-analysis across TCGA cohorts and CPTAC cohorts, plus BRCA's rank within the pan-cancer distribution. |
