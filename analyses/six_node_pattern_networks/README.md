# Three- to six-node composite FFL networks (Supplementary Note S2)

The FFL networks of Bhat et al. (2024) rebuilt on the census graph of this study (`build_bhat_networks.py`), with
topology, hub scores and MCODE clusters.

| File or folder | What it holds |
|---|---|
| `build_bhat_networks.py` | Builds the networks and counts their instances (`bhat_table3_counts.csv`) |
| `instances/` | Every six-node instance, by class |
| `sif/` | The networks in SIF format for Cytoscape, with node attributes |
| `topology_mcc.py` | Topology and MCC hub scores of the merged network |
| `mcode_six_node_composite.R` | MCODE clusters of the six-node composite network |

## Scripts in this folder

| Script | What it does |
|---|---|
| `build_bhat_networks.py` | 3-to-6-node FFL networks of Bhat et al. (2024, Curr Bioinform 19:73-90, doi:10.2174/1574893618666230731164002), rebuilt on the breast-cancer census graph. |
| `mcode_six_node_composite.R` | MCODE on the six-node composite FFL network, as in Bhat et al. |
| `topology_mcc.py` | Topology and MCC hubs of the merged 3-to-6-node composite FFL network, as in Bhat et al. |
