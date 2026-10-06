# Gsea

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_build_ranked_lists.R` | Build the three ranked lists for GSEA (results/v4/rank_*.csv) and their provenance table. |
| `02_build_genesets.R` | Assemble MSigDB collections + the FFL-class node sets as a custom collection |
| `03_run_gsea.R` | GSEA of every ranked list against every collection (fgsea, BH-adjusted) |
| `04_leading_edge_overlap.R` | (A) Leading-edge overlap analysis across the ranked lists |
| `05_ssgsea_scores.R` | (B) ssGSEA / GSVA per-sample scores for the top pathways (survival-arm input) |
| `06_figures.R` | Figures: dot plots, enrichment curves, ridge plots, FFL-class NES |
| `07_ffl_class_gsea.R` | (C) GSEA comparison of the FFL motif classes against the disease signature |
| `08_master_table_and_controls.R` | (D) one tidy master table + the C3:MIR internal positive control |
| `10_final_tables_and_figures.R` | (D) one consolidated tidy results table across every analysis in v4, plus the figure that sets the FFL classes against both backgrounds. |
| `11_gsea_signed_influence.R` | Signed random walk with restart from the top 60 hubs (composite centrality), giving every node a signed influence score, and GSEA on the ranked network nodes under two sign schemes. |
