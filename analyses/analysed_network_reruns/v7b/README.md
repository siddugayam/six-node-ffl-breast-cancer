# V7b

## Scripts in this folder

| Script | What it does |
|---|---|
| `jam_converged_consensus.R` | jActiveModules consensus restricted to the primary seeds whose rank-1 module converged (>= 100 nodes), because without the legacy edges one of the ten seeds stalls at a 2-node module and the all-seed consensus rule … |
| `ora_modules.R` | GO-BP + KEGG over-representation of the module-detection outputs with and without the 30 legacy miRNA-miRNA edges, using the code, universes, thresholds and cached KEGG release of scripts/14_module_detection/go_02_ora.R, and the … |
| `setup_sandbox.sh` | Sandboxes for the BioNet (01, 02, 04) and jActiveModules (01, 02, 05) chains of scripts/14_module_detection. |
