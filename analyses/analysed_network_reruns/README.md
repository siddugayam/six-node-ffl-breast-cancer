# Re-runs on the analysed network

The analysed network excludes the 30 miRNA–miRNA edges of the exemplar circuits. These folders re-run the earlier
analyses without them (`nolegacy`) and, as a check, with them (`full`).

| Folder | What it holds |
|---|---|
| `v3/` | Degree distribution, architecture, controllability and information flow (Supplementary Note S4), with exact control roles |
| `hubboot/` | Hub stability by bootstrap (published Table S3) |
| `tables/` | The interaction table with its `in_analysed_network` flag (published Table S2) and the hub tables without the exemplar edges (published Tables S3 and S6); the file names follow the analysis numbering (see `results/v5/tables/README.md`) |
| `literature/` | Literature volume against network topology |
| `v7/` | Community detection (Supplementary Note S7, Table S9) |
| `v7b/` | BioNet, jActiveModules and GO/KEGG over-representation of the modules |
| `figure_fixes/` | The information-flow figure of the paper (Fig. S5) |

The `sandbox_*` folders mirror the paths that the module-detection scripts (`scripts/14_module_detection/`) expect.
