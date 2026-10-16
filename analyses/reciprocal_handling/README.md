# Reciprocal TF<->miRNA handling in the topological-definition census

Condition codes in this folder follow the topological definition of the paper (Methods 2.2): D1, induced and weakly connected; D2, acyclic; D3, one source and one sink with every vertex on a source–sink path; D4, at least two internally vertex-disjoint source–sink paths.

The published census replaces each of the 1,223 reciprocal TF<->miRNA pairs by its TF -> miRNA arc. This analysis
compares that baseline with three alternatives, at n = 3–6 (exhaustive) and n = 7 (RAND-ESU, the published 8 seeds):
- **U:** no contraction, both arcs kept, so D2 rejects every module that holds both members of a pair;
- **S:** graph-level condensation of the strongly connected components of the census graph;
- **M:** module-level condensation, where the cycles inside each candidate node set are condensed and D3–D4 tested.

For each n it reports the counts, the instance-level overlap (Jaccard) with the baseline, and the fractions of baseline
instances lost and of variant instances gained. The same comparison is made for the six-node composite pattern.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `make_scc_file.py` | Strongly connected components and the condensation DAG of the census graph (`scc_structure.csv`) |
| `census_reciprocal_variants.c` | The published census enumerator with the four treatments; mode `base` reproduces the published counts |
| `run_reciprocal_variants.sh` | Runs every enumeration |
| `summarise_reciprocal_variants.py` | Compares the module sets exactly (packed vertex keys) and writes the tables |
| `six_node_pattern_variants.py` | The six-node composite pattern under each treatment (recount with `s1_null6.c`) |
| `check_small_n.py` | An independent Python implementation (undirected ESU, its own D2–D4 test) for n = 3 and 4 |
| `write_values.py` | Writes `VALUES.md` |
| `census_reciprocal_variants.csv` | Counts, shared, Jaccard, lost and gained per n and treatment |
| `lost_under_module_condensation_first_condition.csv` | Baseline modules lost, by first failing condition |
| `composite_three_node_under_module_condensation.csv` | The composite three-node FFLs under M |
| `six_node_pattern_variants.csv` | The six-node pattern under each treatment |
| `scc_structure.csv` | SCC sizes and node types |
| `check_small_n_output.txt` | Output of the independent check |
| `counts/` | Raw output of every enumeration |

**Run:** `FFL_ROOT=<repository checkout> bash run_reciprocal_variants.sh`, then
`python3 summarise_reciprocal_variants.py`, `python3 six_node_pattern_variants.py`,
`python3 check_small_n.py 3 4 > check_small_n_output.txt` and `python3 write_values.py`. The module key files, about
3 GB, are written to `work/` and are not part of the folder.

## Scripts in this folder

| Script | What it does |
|---|---|
| `census_reciprocal_variants.c` | the topological definition census of the census graph under four treatments of the reciprocal TF<->miRNA pairs (and, for the condensation variants, of every directed cycle). |
| `check_small_n.py` | Independent check of census_reciprocal_variants.c at small n. |
| `make_scc_file.py` | Strongly connected components of the census graph and its condensation DAG. |
| `run_reciprocal_variants.sh` | the topological definition census of the census graph under four treatments of reciprocal TF<->miRNA pairs (see SETTINGS.md): base (contracted, the published census), unc (uncontracted), scc (graph-level SCC condensation), mod (module-level … |
| `six_node_pattern_variants.py` | The six-node composite pattern under the four treatments of reciprocal pairs. |
| `summarise_reciprocal_variants.py` | Summary of the census under four treatments of reciprocal TF<->miRNA pairs (see SETTINGS.md). |
| `write_values.py` | Writes VALUES.md for the reciprocal-handling census from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
