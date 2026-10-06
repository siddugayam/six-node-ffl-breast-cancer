# Edge correlations

## Scripts in this folder

| Script | What it does |
|---|---|
| `02_edge_corr_caf.R` | Independent verification of N6 / N8 / N9 (tasks B and C) Edge-level Spearman correlations in TCGA-BRCA primary tumours, expression-decile-matched null, and non-circular CAF (stromal) adjustment. |
| `03_null_constructions.R` | Sensitivity of the N8 null to how the decile-matched random pairs are drawn. |
| `03b_null_within_network.R` | 4th null construction: random pairs drawn from the NETWORK'S OWN node set, decile-matched. |
| `04_caf_adjusted_vs_matched_null.R` | N9 under the SAME (within-network-node, decile-matched) null used by the paper. |
| `07_consolidate.R` | Consolidated table of the expression tests of the named regulatory axes (results/v2/verify_expression.csv). |
| `07_genomewide_direction_concordance.R` | Is the direction disagreement between HPA's breast prognostic calls and our own TCGA/METABRIC Cox specific to our 26 genes, or genome-wide? |
| `08_expression_validation.R` | Expression-based validation of PREDICTED REGULATORY RELATIONSHIPS in TCGA-BRCA. |
| `09_anova_ffl_named_axes.R` | Comparison of TF and miRNA log fold changes, sign coherence of the three-node FFL cores, and correlations of the named regulatory axes with partial-correlation sensitivity analyses. |
| `10_directional_evidence_and_anova_control.R` | Supplementary to 08/09: (i) per-stratum counts of edges whose correlation is SIGNIFICANT and in the predicted direction vs significant and in the OPPOSITE direction, real edges vs matched nulls -> … |
| `10a_export_control_data.R` | 10a: export, as flat TSVs, the three edge sources used to rebuild BOTH the breast-cancer network and the disease-unrelated matched control networks under an IDENTICAL pipeline: TRRUST TF -> target (local) TransmiR TF -> miRNA … |
| `11_edge_rho_all_edges.R` | Independent re-computation (this run) of edge-level Spearman rho in TCGA-BRCA for every miRNA_target and TF_target edge of the canonical network. |
