# Mirna survival meta analysis

## Scripts in this folder

| Script | What it does |
|---|---|
| `c01_geo_search.py` | Systematic GEO (E-utilities) + ArrayExpress/BioStudies search for breast cancer miRNA datasets. |
| `c02_geo_triage.py` | Triage the GEO search hits: keep human breast miRNA series with n>=40, then probe each for outcome fields. |
| `c03b_fetch_soft_bulk.py` | Downloads the GEO SOFT descriptions of the series listed in the file given as argument into cache/v6/cohorts/soft/, skipping files already present. |
| `c04_soft_summary.py` | Screens the GEO SOFT descriptions in cache/v6/cohorts/soft/ for outcome annotation (survival, relapse, response) and prints a summary for each series. |
| `c05_fetch_matrices.py` | Download GEO series matrices + platform annotations for the usable miRNA cohorts. |
| `c06_extract_gpl.py` | Extract probe -> miRNA-name maps from GEO platform SOFT files. |
| `c07_build_mirna_cohorts.py` | Parse GEO series matrices for breast-cancer miRNA cohorts, map probes to the 11 target miRNAs, extract clinical/outcome, and write tidy per-cohort CSVs. |
| `c08_make_cohorts.py` | Builds the miRNA survival cohorts from GEO series matrices (expression, probe mapping, endpoint and event) and records the probe mapping in cache/v6/cohorts/built/probe_report.json. |
| `c09_make_cohorts2.py` | Builds the remaining miRNA survival cohorts in the same way as c08_make_cohorts.py and adds them to the probe report. |
| `c10_dump_full_matrices.py` | Dump collapsed mature-miRNA matrices + pheno for the cohorts new to v6. |
| `c11_mirna_cox_meta.R` | Cox models for the 10 prioritised miRNAs (+ miR-29a control) in every breast-cancer miRNA cohort with an outcome, then DerSimonian-Laird random-effects meta-analysis per miRNA. |
| `c12_mirna_binary.R` | Secondary: breast-cancer miRNA cohorts whose outcome is a binary flag with no follow-up time -> logistic OR per SD, meta-analysed separately from the Cox HRs. |
| `c13_fetch_mrna_cohorts.py` | Downloads the GEO series matrices of the mRNA cohorts into cache/v6/cohorts/matrix/, skipping files already present. |
| `c15_pancancer_prep.sh` | Extract the rows we need from the 331MB pan-cancer Xena gene-expression matrix. |
| `c16_pancancer_axes.R` | miR-29a -\| COL1A1/COL3A1 and miR-101 -\| EZH2 in every TCGA cohort with enough samples; per-cohort CAF score; does axis strength track CAF content? |
| `c17_mrna_cohorts.R` | Seven additional breast-cancer mRNA cohorts with survival, none of which is already in the study: MAINZ, TRANSBIG, VDX, UPP, UNT, NKI (Bioconductor breastCancer* experiment packages) and GSE1456 (Stockholm) from GEO. |
| `c18_liquid_biopsy.R` | Are the prioritised miRNAs detectable and differential in serum / plasma of breast cancer patients? |
| `c19_liquid_biopsy2.R` | GSE73002 serum: detection and case–control discrimination for the prioritised miRNAs, benchmarked against all 2,540 assayed miRNAs and repeated after within-array rank normalisation. |
| `c20_mirna_inventory.py` | Inventory of the breast cancer miRNA cohorts found in GEO, with their outcome fields and whether each could be used. |
| `c21_liquid_other.R` | Three further circulating-miRNA cohorts. |
| `c22_autoclassify.py` | Auto-classify every probed breast miRNA series by the outcome fields actually deposited in GEO: time+event (Cox-able), event-only (logistic), or none. |
| `c23_addq.R` | Adds Benjamini-Hochberg q values across the miRNAs tested to the pooled survival meta-analysis (results/v6/mirna_meta_pooled.csv) and to the forest-plot table. |
| `c24_sensitivity.R` | Sensitivity analyses of the miRNA survival meta-analysis: pooled estimates for alternative cohort sets, fixed- against random-effects estimates and Egger tests; writes results/v6/mirna_meta_sensitivity.csv and … |
