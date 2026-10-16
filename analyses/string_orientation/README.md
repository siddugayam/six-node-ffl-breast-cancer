# Orientation of the STRING associations in the census graph

The census graph enters each of the 941 STRING associations as a single arc in its listed orientation. This analysis:
- **(a)** orients them at random (500 replicates) and repeats the census of the four conditions at n = 3–6;
- **(b)** counts the census with TRRUST-only gene–gene links and with no gene–gene layer;
- **(c)** counts the modules that traverse no STRING arc.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `string_orientation_graphs.py` | Rebuilds the census graph from the network files, checks it against the published graph, and writes the variant and bootstrap graphs |
| `census_string_variants.c` | The published census enumerator, with per-module STRING counts and a look-up of each module among the published ones |
| `run_string_orientation.sh` | Runs everything |
| `summarise_string_orientation.py` | Checks the published counts and writes the tables |
| `write_values.py` | Writes `VALUES.md` |
| `string_bootstrap_summary.csv`, `string_bootstrap_counts.csv`, `string_bootstrap_graphs.csv` | (a) |
| `gene_gene_layer_variants.csv` | (b) |
| `orientation_robust_core.csv` | (c) |
| `counts/` | Raw counts of every enumeration |

**Run:** `FFL_ROOT=<repository checkout> bash run_string_orientation.sh`, then
`python3 summarise_string_orientation.py` and `python3 write_values.py`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `census_string_variants.c` | Census of the four conditions of the topological definition of a census graph, with STRING-related counts per module. |
| `run_string_orientation.sh` | Census counts for the STRING orientation analysis (see SETTINGS.md). |
| `string_orientation_graphs.py` | Census graphs with the STRING associations re-oriented, and the gene-gene layer variants. |
| `summarise_string_orientation.py` | Summary of the STRING orientation analysis (see SETTINGS.md), from counts/ written by run_string_orientation.sh. |
| `write_values.py` | Writes VALUES.md for the STRING orientation analysis from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
