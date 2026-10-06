# Multi cohort

## Scripts in this folder

| Script | What it does |
|---|---|
| `M1_build_genesets.R` | Define every gene set used by the multi-cohort expression / survival meta-analysis and write them to results/v4/multicohort_genesets.csv (+ an RDS for the analysis scripts). |
| `M2_stream_scanb.py` | Stream the 592 MB SCAN-B (GSE96058) expression CSV and keep the rows for the 1,901 symbols the multi-cohort module analysis needs. |
| `M3_build_mirna_cohorts.R` | Assemble every breast cohort on disk that carries a miRNA assay, with its phenotype/outcome table. |
| `M4_build_mrna_cohorts.R` | Take the 8 cohorts assembled by the parallel v2 workflow (cache/newcohorts/cohorts_base.rds), replace the 173-gene SCAN-B slice with the 1,897-gene slice streamed in M2, and attach ER / subtype / covariates. |
| `M5_analysis_A.R` | A) DE-direction concordance with TCGA for the network hubs (mRNA + miRNA) |
| `M6_analysis_C_D_E.R` | C/D) Cox + random-effects meta-analysis across every mRNA cohort with an outcome for the miR-29 target / ECM module and for the 3-node vs higher-order FFL module scores, univariate and adjusted for age/grade/size/node E) … |
| `M7_analysis_E_extras.R` | E) subtype-stratified module scores (basal/TNBC vs LumA) + head-to-head 3-node vs higher-order-only FFL modules in one Cox model + module scores in tumour vs non-tumour breast (GTEx / TCGA normals) |
| `M8_analysisA_forest.R` | A3) hub DE direction: TCGA tumour vs GTEx normal breast (cross-study), re-derived here rather than taken from the v2 workflow A4) consolidate every direction-concordance test into one table Z ) one forest-plot-ready table … |
| `M9_link22216_figures.R` | (i) prove the GSE22219 (mRNA) <-> GSE22216 (miRNA) patient linkage from the deposited clinical fields (ii) forest-plot figures for the headline meta-analyses |
| `M12_cohort_inventory.R` | M12 -- one consolidated cohort inventory for the multi-cohort arm: accession, n, platform, survival availability, miRNA availability. |
