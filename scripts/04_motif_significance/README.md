# Motif significance (Results 3.2)

`08_motif_significance.py` / `30_motif_significance.py`: three-node FFL counts in the observed network against randomised
networks (first analysis round); the three null models of the paper (Table 2, Fig. 3a) are run in
`analyses/census_and_motif_nulls/` (`nulls/v2_null.c`, `nulls/summarise_nulls.py`). `08c_reciprocity_audit.py`: reciprocal TF↔miRNA pairs in the network and in the source databases;
`08_motif_node_sets.py`: the node sets of each FFL class; `10_*`, `11_*`: class comparison and specificity. The
analysed-network re-run is `analyses/census_and_motif_nulls/nulls/`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `08_motif_node_sets.py` | Build the node sets used for functional enrichment of each FFL motif class. |
| `08_motif_significance.py` | Motif significance testing of the canonical breast-cancer regulatory network against degree- and type-preserving randomised networks. |
| `08a_control_pool_fetch.R` | 08a: build the node pools for the DISEASE-UNRELATED matched control networks and retrieve their multiMiR *validated* miRNA-target edges (genome-wide per gene). |
| `08c_reciprocity_audit.py` | Why the Composite-FFL result depends entirely on the null model. |
| `10_motif_class_comparison.R` | Produces results/motif_class_comparison.csv - a long tidy table, one row per comparison statistic, with columns: block, metric, class_a, class_b, value, statistic, p_value, p_adjust, n_a, n_b, note |
| `10_motif_specificity.py` | SPECIFICITY control: is the FFL structure specific to breast cancer, or does it merely reflect the general density of known molecular interactions? |
| `11_specificity_enrichment.py` | The scale-free half of the specificity control. |
| `30_motif_significance.py` | Motif over-representation against randomised networks. |
