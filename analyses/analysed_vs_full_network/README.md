# Analysed network compared with the full network

The analysed network excludes the 30 miRNA–miRNA edges of the exemplar circuits. These folders run the network
analyses without them (`noexemplar`) and, as a check, with them (`full`).

| Folder | What it holds |
|---|---|
| `systems_analyses/` | Degree distribution, architecture, controllability and information flow (Supplementary Note S4), with exact control roles |
| `hubboot/` | Hub stability by bootstrap (published Table S3) |
| `tables/` | The interaction table with its `in_analysed_network` flag (published Table S2) and the hub tables without the exemplar edges (published Tables S3 and S6); the file names follow the analysis numbering (see `results/prioritisation_and_stroma/tables/README.md`) |
| `literature/` | Literature volume against network topology |
| `clustermaker2_modules/` | Community detection (Supplementary Note S6, Table S9) |
| `jam_and_bionet_modules/` | BioNet, jActiveModules and GO/KEGG over-representation of the modules |
| `information_flow_figure/` | The information-flow figure of the paper (Fig. S5) |

The `sandbox_*` folders mirror the paths that the module-detection scripts (`scripts/14_module_detection/`) expect.
