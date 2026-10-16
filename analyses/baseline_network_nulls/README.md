# A baseline network without FFL pre-selection

The curated gene and miRNA lists from which the network was assembled (1,057 genes, 421 miRNAs) were not saved. This analysis keeps the 587 nodes of the analysed network and takes every edge among them directly from the
primary interaction tables, not from the exported FFL networks:
- TF -> gene/TF: TRRUST version 2;
- TF -> miRNA: TransmiR version 2.0;
- miRNA -> target: the multiMiR validated tables, and, in a second network, also the TargetScan 8.0 default
  predictions.

The three-node FFL motif test is then run on each network under the degree-preserving, reciprocity-preserving and
bipartite-layer null models, with the published TF typing and with the 22 co-regulators typed as genes.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `build_baseline_network.R` | Builds the networks and checks the TRRUST and TransmiR arcs against the project's layers |
| `motif_null_models_large.c` | `motif_null_models.c` with a larger arc array (the networks have about 33,000 arcs) |
| `run_baseline_nulls.sh` | Builds the networks and runs the 12 null runs (1,000 replicates each) |
| `summarise_baseline_nulls.py` | Summary and decision, beside the analysed network |
| `write_values.py` | Writes `VALUES.md` |
| `network_*.txt` | The four network files in the null-model format (node types and arcs by index). The build also writes each file's node names to `<file>.meta.json`; those lists are not included here |
| `network_arc_counts.csv`, `analysed_network_overlap.csv` | Network composition and overlap with the analysed network |
| `null_model_summary.csv`, `decision.csv` | The motif test |
| `runs/` | Replicate files and logs |

**Run:** `FFL_ROOT=<repository checkout> bash run_baseline_nulls.sh`, then `python3 summarise_baseline_nulls.py` and
`python3 write_values.py`. The build reads the TRRUST and TransmiR downloads (`data/db/`), the multiMiR query cache
(`cache/multimir/`) and the TargetScan files (`cache/direct/ts/`).

## Scripts in this folder

| Script | What it does |
|---|---|
| `build_baseline_network.R` | a baseline regulatory network on the 587 nodes of the analysed network, with every edge among them taken directly from the primary interaction tables rather than from the exported FFL networks (see SETTINGS.md): TF -> gene/TF … |
| `motif_null_models_large.c` | Analyses/census_and_motif_nulls/nulls/motif_null_models.c copied unchanged except for the arc array size (MAXE 20000 -> 40000), so that the baseline networks of about 33,000 arcs fit. |
| `run_baseline_nulls.sh` | Three-node FFL motif test on the baseline networks (see SETTINGS.md): the degree-preserving, reciprocity-preserving and bipartite-layer null models of motif_null_models.c (program arguments A, B and C), 1,000 randomised networks … |
| `summarise_baseline_nulls.py` | Summary of the three-node FFL motif test on the baseline networks, beside the analysed network (see SETTINGS.md). |
| `write_values.py` | Writes VALUES.md for the baseline network analysis from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
