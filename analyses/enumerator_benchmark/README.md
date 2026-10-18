# Benchmark of the topological-definition enumerator on planted modules

Condition codes follow the topological definition of the paper (Methods 2.2): D1, induced and weakly connected; D2, acyclic; D3, one source and one sink with every vertex on a source–sink path; D4, at least two internally vertex-disjoint source–sink paths.

The census enumerator (source-rooted ESU with acyclicity pruning and the topological-definition test of
`analyses/census_and_motif_nulls/census/ffl_census_composition.c`) is run on 200 synthetic graphs:
- **background:** the analysed network, rewired with its degree sequence kept and made acyclic;
- **planted modules:** 10 topological-definition modules per n = 3–6;
- **noise:** arcs equal to 10 % of the background edges.

The enumerator's output is compared with:
- the planted modules;
- an independent topological-definition test written in Python (networkx);
- an independent exhaustive enumeration, at n = 3 for all replicates and at n = 4 for the first five.

Two relaxed tests are run beside the full one: D3 without "every vertex on a source–sink path", and D4 with one
source–sink path instead of two vertex-disjoint ones.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `census_benchmark.c` | The census enumerator, with module-key output and the two relaxed tests selected by an argument |
| `benchmark.py` | Builds the graphs, plants the modules, runs the enumerator and the independent checks, and writes the tables |
| `write_values.py` | Writes `VALUES.md` |
| `planted_topologies.csv` | The number of the topological definition topologies per n |
| `benchmark_per_replicate.csv` | Every measure per replicate |
| `benchmark_summary.csv` | Median, min, max and mean of each measure over the 200 replicates |

**Run:** `FFL_ROOT=<repository checkout> python3 benchmark.py 200 <parallel processes>`, then `python3 write_values.py`.
The run needs about 27 processor-minutes and up to 8 GB of memory per replicate. Finished replicates are kept in
`work/`, so a stopped run resumes.

## Scripts in this folder

| Script | What it does |
|---|---|
| `benchmark.py` | Benchmark of the topological-definition census enumerator on planted modules (see SETTINGS.md). |
| `census_benchmark.c` | The census enumerator of analyses/census_and_motif_nulls/census/ffl_census_composition.c (source-rooted ESU over out-neighbourhoods, acyclicity pruning, D2-D4 test of check_module(); exhaustive), writing the vertex set of every … |
| `write_values.py` | Writes VALUES.md for the benchmark of the topological-definition enumerator from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
