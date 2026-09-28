# Three- to six-node composite FFL networks (Supplementary Note S2)

The FFL networks of Bhat et al. (2024) rebuilt on the census graph of this study (`build_bhat_networks.py`), with
independent recounts, topology, hub scores and MCODE clusters.

| File or folder | What it holds |
|---|---|
| `build_bhat_networks.py` | Builds the networks and counts their instances (`bhat_table3_counts.csv`) |
| `instances/` | Every six-node instance, by class |
| `sif/` | The networks in SIF format for Cytoscape, with node attributes |
| `check_counts_independent.py`, `check_topology_igraph.R` | Independent recounts and an igraph cross-check |
| `topology_mcc.py` | Topology and MCC hub scores of the merged network |
| `mcode_six_node_composite.R` | MCODE clusters of the six-node composite network |
| `b5_reconcile.py` | Reconciles the six-node composite count with the census graphs |
