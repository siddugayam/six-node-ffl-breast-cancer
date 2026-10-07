# E5

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_coherence_corrected_E5.py` | Coherent / incoherent typing (Mangan & Alon, PNAS 2003) on the unique cores (01a convention). |
| `e5a_hypergeometric_filter.py` | E5(a): which TF list entered the hypergeometric miRNA-TF pair filter, and are any of the 22 TF-typed nodes outside the Lambert census among the retained TFs? |
| `e5b_null_summary.py` | E5(b) nulls: NULL-A, NULL-B, NULL-C (analyses/census_and_motif_nulls/nulls/v2_null.c, 1,000 randomisations, seed 20250908, 100 swaps per edge) on the analysed network (6,829 edges) with the 22 non-Lambert TF-typed nodes re-typed … |
| `e5bc_retype_graphs.py` | E5(b)/(c) inputs: the 22 TF-typed nodes absent from the Lambert census re-typed as genes (type code 1 -> 2) in null_dep_nolegacy_retyped22.txt = analyses/census_and_motif_nulls/graphs/null_dep_nolegacy.txt (analysed network, … |
