# S7

## Scripts in this folder

| Script | What it does |
|---|---|
| `ffl_census_composition.c` | Counting-only n-node FFL census (C), used for the census re-runs. |
| `s7b_restricted_gene_gene.py` | S7(b): BHAT6 (all classes) and MODEL6 with the G1-G2 link restricted to (i) STRING v12.0 physical-subnetwork links with score >= 0.900 (ii) STRING v12.0 co-expression-channel scores >= 0.900 Downloaded public files (not in the … |
| `s7c_build_validated_graphs.py` | S7(c) validated-only census graphs (analyses/six_node_pattern). |
| `s7c_summarise.py` | S7(c) summary: exhaustive census (ffl_census_composition, max-flow D4) at n = 3-6 on the validated-only graphs with and without STRING; NULL-B / NULL-C (v2_null.c, 1,000 randomisations, seed 20250908, 100 swaps per edge) for … |
| `v2_null.c` | Motif significance for 3-node FFLs under three null models. |
