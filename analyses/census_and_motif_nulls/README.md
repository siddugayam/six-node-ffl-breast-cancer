# Census of n-node FFL modules and three-node motif null models (Table 2, Figs 3a and S2, Supplementary Note S2)

`run_all.sh` compiles and runs everything here; `jobs.txt` lists every command it ran and `run_all.log` its log.

| Folder | What it holds |
|---|---|
| `graphs/` | The census graphs (with and without the 30 exemplar miRNA–miRNA edges, with and without STRING arcs) and the network graphs of the null-model runs, written by `make_graphs.py` |
| `census/` | `ffl_census_composition.c` (exhaustive census at n = 3–6, RAND-ESU sampling at n = 7), one subfolder of counts per graph, and `summarise_census.py`, which writes `census_all_graphs.csv` (the counts of Supplementary Note S2 and Fig. S2a) |
| `nulls/` | `v2_null.c`, the three-node motif test under three null models (1,000 randomised networks each), its outputs and `summarise_nulls.py`, which writes `motif_nulls_all_graphs.csv` |

`census_all_graphs.csv` and `motif_nulls_all_graphs.csv` are also in `results/v2/census_rerun/`, where the figure
scripts read them. Graph `nolegacy` is the census graph of the paper.
