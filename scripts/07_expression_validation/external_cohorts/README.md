# External cohorts

## Scripts in this folder

| Script | What it does |
|---|---|
| `14_external_met500.R` | EXTERNAL COHORT 1: MET500 metastatic breast cancer (local files). normal DE in primary TCGA tumours cannot support a statement about metastasis. |
| `16_external_geo_de.R` | EXTERNAL COHORT 3 (+4): GEO breast tumour-vs-normal series on Affymetrix GPL570. |
| `17_external_cohort_validation.R` | Consolidates every external-cohort comparison into one table and computes replication rates: of the features significant in TCGA-BRCA, how many replicate in the same direction in an independent cohort. |
| `N0_geo_utils.R` | Self-contained GEO series-matrix parser and probe collapse for the external-cohort replication. |
| `N1_build_cohorts.R` | Standardise every transcriptomic cohort into list(name, accession, platform, X = symbols x samples (log2), pheno, notes) and save cache/newcohorts/cohorts.rds Every number printed here is computed in this run. |
| `N1b_build_geo.R` | Parse the GEO array series into standardised cohorts, collapse probes to gene symbols, extract survival + clinical covariates. |
| `N1c_build_geo2.R` | N1c: add GSE20685 (GPL570) and GSE22219 (GPL6098) to the cohort list |
| `N1d_build_scanb.py` | Stream the 592 MB SCAN-B (GSE96058) expression CSV and keep only the rows needed for the replication analysis (hubs + CAF signatures + collagens/TFs). |
| `N1e_build_scanb.R` | Add SCAN-B / GSE96058 (3,273 tumours, RNA-seq) to the cohort list, from the streamed gene subset + the two series matrices. |
| `N1f_fix_covariates.R` | Repair three clinical-covariate defects found by auditing the Cox convergence warnings emitted by N2. |
| `N2_cohort_analysis.R` | Replication of the TCGA-BRCA findings in independent transcriptomic cohorts. |
| `N3_meta.R` | Random-effects (DerSimonian-Laird) meta-analysis across the independent cohorts of the ETS1->COL1A1 / NFKB1->COL1A1 correlations (unadjusted and CAF-adjusted) and of COL1A1~COL3A1. |
| `N4_normal_vs_tumour.R` | Contrast the GTEx normal-breast correlations with the random-effects pooled tumour estimate for the same edge, and append the rows to newcohorts_summary.csv. |
| `N5_headline.R` | Prints the headline results of the external-cohort replication from results/v2/newcohorts_*.csv (cohort inventory, direction concordance of hubs and edges, meta-analysis). |
| `N6_forest.R` | Forest plots from newcohorts_meta.csv (base graphics, no deps) |
| `N8_inventory.R` | Append an explicit per-cohort inventory to results/v2/newcohorts_summary.csv: accession, n, platform, genes measured, survival endpoints and event counts, covariates usable, and an explicit statement of which of analyses (a)-(e) … |
| `N9_gtex_sex.R` | GTEx "breast - mammary tissue" contains MALE donors. |
| `run_newcohorts.sh` | Full new-cohort replication pipeline (v2). |
