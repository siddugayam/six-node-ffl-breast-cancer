# Metabolic

## Scripts in this folder

| Script | What it does |
|---|---|
| `20_metabolic_genesets.R` | Build the metabolic gene-set space for the metabolic-arms analysis. |
| `21_metabolic_scores.R` | Per-sample ssGSEA + GSVA scores for the metabolic gene-set space, over all 1,211 TCGA-BRCA samples (1,097 primary tumours + 114 normals). |
| `22_metabolic_tumour_vs_normal.R` | Tumour vs normal for every metabolic score, with FDR. |
| `23_metabolic_associations.R` | Correlate every metabolic score with the miR-29 family, the let-7 family and the FFL 3-node / higher-order module scores. |
| `24_metabolic_three_claims.R` | Tests of three metabolic associations: glycolysis and let-7, folate metabolism, and AGE–RAGE signalling partners. |
| `25_metabolic_module_overlap_control.R` | The FFL module scores correlate with metabolic scores at rho up to 0.87. |
| `26_metabolic_survival.R` | Cox models for every metabolic score (OS, PFI, DSS), per SD, adjusted for age at diagnosis and AJCC pathologic stage. |
| `27_metabolite_consistency.R` | Metabolite-level evidence. Re-derived independently from the cached Metabolomics Workbench ST000054 / AN000092 mwTab already on disk (cache/metabolomics/an92.mwtab.txt) -- no re-download. |
| `29_metabolic_summary.R` | Summary tables of the three metabolic associations (one row each) and the headline tables. |
| `E1_metabolomics.R` | MEASURED metabolomics (not expression-inferred) |
| `E2_pathway_gsva.R` | Metabolic pathway activity inferred from mRNA expression. |
| `E3_nrf2_keap1.R` | KEAP1–NFE2L2 in the network and in TCGA-BRCA: network membership, NRF2 target-signature scores against random signatures, expression and mutation, and mutation frequency compared with lung cancers. |
| `E4_let7b_hk2_control.R` | Specificity control for the association of let-7b with the glycolysis score (let-7b and HK2). |
