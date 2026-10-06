# V7

## Scripts in this folder

| Script | What it does |
|---|---|
| `cd_audit.R` | Community-detection results of Supplementary Note S6 and Table S9, with and without the 30 legacy GeneMANIA miRNA–miRNA edges. |
| `seed_frequency.R` | How often do COL1A1, COL3A1 and miR-29a share a module across the stochastic seeds that 04_igraph_battery.R stores (50 Louvain, 20 Infomap, 100 label-propagation), with and without the 30 legacy edges. |
| `setup_sandbox.sh` | Build sandbox_<mode>/ mirroring the paths scripts/14_module_detection/04_igraph_battery.R expects, so that the script runs unchanged except for its setwd() line. |
