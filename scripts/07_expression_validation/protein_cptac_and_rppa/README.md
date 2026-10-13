# Protein cptac and rppa

## Scripts in this folder

| Script | What it does |
|---|---|
| `06_cptac.R` | CPTAC-BRCA replication of the miR-29 -> collagen axis at PROTEIN level, and global miRNA-target concordance mRNA vs protein. |
| `20_cptac_prep.R` | Harmonise CPTAC-BRCA proteome / phosphoproteome / RNA / miRNA Output: data-side RDS bundle in results/multiomics/cptac_bundle.rds (kept inside multiomics/) |
| `21_cptac_protein_validation.R` | A) miRNA->target: mRNA vs PROTEIN correlation, decile-matched null, by evidence tier B) TF->target at protein level, split by TRRUST mode C) named axes at protein level D) per-gene mRNA-protein concordance E) NF-kB phospho-site … |
| `22_cptac_axes_mrnaprot_phospho.R` | C) named axes at protein level D) per-gene mRNA-protein concordance E) NF-kB phospho-site activity vs collagen protein plus sensitivity analyses for A/B (TF mRNA -> target protein ; complete-protein subset) |
| `23_rppa_second_platform.R` | TCGA RPPA as an independent protein platform |
| `24_nfkb_activity_score.R` | NF-kB TRANSCRIPTIONAL ACTIVITY (target-gene signature) vs collagen, because RELA has no phospho-site in the CPTAC-BRCA phosphoproteome. |
