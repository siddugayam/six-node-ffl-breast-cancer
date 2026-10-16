# The six-node pattern and acyclicity

Condition codes in this folder follow the topological definition of the paper (Methods 2.2): D1, induced and weakly connected; D2, acyclic; D3, one source and one sink with every vertex on a source–sink path; D4, at least two internally vertex-disjoint source–sink paths.

- **(a)** The 2,498 node sets of the six-node composite pattern are tested against the topological definition, giving the first failing
  condition and the kind of cycle.
- **(b)** Six-node modules containing the composite core plus a second TF are counted under D2′, which relaxes D2 to
  allow exactly one cyclic component among {TF1, TF2, miRNA}. The count is tested against the bipartite-layer null and the core-conditional null.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `six_node_d2prime_nulls.c` | `analyses/six_node_pattern/S1/s1_null6.c` with one count added (D2′); its other output lines are unchanged |
| `run_d2prime_nulls.sh` | Compiles it, counts the observed graph, and runs the bipartite-layer null and the core-conditional null (1,000 replicates each) |
| `pattern_first_failing.py` | Part (a) |
| `summarise_d2prime.py` | Part (b): summary and decision. It checks that the replicates are the graphs of the S1 runs |
| `write_values.py` | Writes `VALUES.md` |
| `pattern_first_failing_condition.csv` | (a) node sets by first failing condition and cycle class |
| `d2prime_counts_observed.csv`, `d2prime_null_summary.csv`, `d2prime_decision.csv` | (b) |
| `check_replicates_against_s1.txt` | The replicate check |
| `runs/` | Observed and per-replicate counts |

**Run:** `FFL_ROOT=<repository checkout> bash run_d2prime_nulls.sh`, then `python3 pattern_first_failing.py`,
`python3 summarise_d2prime.py` and `python3 write_values.py`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `pattern_first_failing.py` | Part (a): the six-node composite pattern against the topological definition (paper numbering: D1 induced and weakly connected; D2 acyclic; D3 one source and one sink, every vertex on a source-sink path; D4 two internally vertex-disjoint source-sink … |
| `run_d2prime_nulls.sh` | D2' six-node counts on the census graph (observed) and under the bipartite-layer null and the core-conditional null (see SETTINGS.md): s1_null6.c with the added D2' count, seed 20250908, 100 swaps or trades per item, R replicates per model split over P processes … |
| `six_node_d2prime_nulls.c` | s1_null6.c of analyses/six_node_pattern/S1 (copied below this note; four comments reworded) with one added count per graph, D2PRIME: six-node modules that contain the composite core plus a second TF (the CENSUS6c seed: … |
| `summarise_d2prime.py` | Part (b) summary: D2' six-node counts on the census graph against the bipartite-layer null and the core-conditional null (see SETTINGS.md). |
| `write_values.py` | Writes VALUES.md for the six-node pattern and acyclicity analysis from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
