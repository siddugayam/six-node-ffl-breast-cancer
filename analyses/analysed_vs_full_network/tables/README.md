# Tables

## Scripts in this folder

| Script | What it does |
|---|---|
| `flag_all_interactions.py` | Table S2 with the analysed-network flag: every row of results/prioritisation_and_stroma/tables/all_interactions.csv (6,859 deposited edges) plus 'in_analysed_network' (FALSE for the 30 GeneMANIA miRNA_miRNA edges of the exemplar circuits, TRUE for the 6,829 analysed edges). |
| `build_directed_hub_table.py` | Rebuild Table S6 (direction-aware hub table), for all 587 nodes ordered by undirected current-flow betweenness, from the Note S4 information-flow and controllability outputs (../systems_analyses/out_<mode>_hash0/). |
