# Cell type sources

## Scripts in this folder

| Script | What it does |
|---|---|
| `ct_01_hpa.py` | Human Protein Atlas single-cell RNA: cell-type source of network genes. |
| `ct_03_stromal_scores.R` | Stromal / fibroblast / epithelial content scores for TCGA-BRCA primary tumours. |
| `ct_04_celltype_source.py` | Per-cell-type expression of the network's protein-coding nodes in 26 primary breast tumours (GSE176078, Wu et al. 2021, 100,064 cells), combined with the Human Protein Atlas normal-breast single-cell reference. |
| `ct_05_caf_signature.py` | Derive a data-driven CAF signature from GSE176078 that is deliberately COLLAGEN-FREE and NETWORK-NODE-FREE, so scoring TCGA bulk with it and then partialling it out of the collagen correlations cannot be circular. |
| `ct_06_partial_corr.R` | Does the collagen signal survive adjustment for stromal / CAF content? |
| `ct_07_mirna_hostgenes.py` | Cellular source of the miRNA loci, approached through their pri-miRNA HOST GENES. |
| `ct_08_network_stromal.R` | Network-wide: does the sign concordance of the edges survive adjustment for CAF content? |
| `ct_09_external_corr.R` | Independent-cohort replication of (i) the collagen/TF correlations and (ii) the miR-29 -> collagen axis, each with and without adjustment for CAF content. |
| `ct_10_metabric_survival.R` | C) Survival replication in METABRIC: hub set + collagen/TF module scores, with and without adjustment for CAF (stromal) content. |
| `ct_11_power_downsample.R` | Is the TCGA "no hub is prognostic" null a power artefact? |
| `ct_12_replication_summary.R` | Consolidated independent-cohort replication table. |
| `ct_13_stratified.R` | Stratified check: partial correlation on a covariate as collinear as the CAF score (rho 0.88 with COL1A1) can over-adjust. |
