# Cox models

## Scripts in this folder

| Script | What it does |
|---|---|
| `03_cox_crosscheck.R` | Cross-check HPA pathology-atlas prognostic calls against our own TCGA-BRCA and METABRIC Cox / best-cutoff results for the same genes. |
| `10_survival_cox_hubs.R` | Univariable Cox regression for every node of the canonical FFL network in TCGA-BRCA primary tumours. |
| `12_module_score_survival.R` | Module-level scores (GSVA + mean-z) -> Cox, KM median split, multivariable Cox adjusted for age and AJCC stage. |
| `13_circuit_level_3node_vs_6node.R` | The aggregate-score test (script 12) asks whether the *union* of higher-order nodes carries information. |
| `15_external_metabric.R` | EXTERNAL COHORT 2: METABRIC (cBioPortal datahub study brca_metabric), Illumina HT-12 v3 microarray, n=1980 tumours with expression, OS and RFS. |
| `18_metabric_treatment.R` | METABRIC records which patients actually received chemotherapy, hormone therapy and radiotherapy, so the closest legitimate computational test is: (i) is the hub / module score prognostic WITHIN the treated subgroup, and (ii) is … |
