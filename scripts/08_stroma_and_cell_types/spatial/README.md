# Spatial transcriptomics (Visium)

`05_spatial_full.py` wrote the six `results/v6/spatial_*` tables that Fig. 4 and `06_spatial_mediation_summary.py` read.

## Scripts in this folder

| Script | What it does |
|---|---|
| `05_spatial_full.py` | Spatial transcriptomics test of the compartment argument (full version). |
| `05b_spatial_sp1_rela_controls.py` | Visium compartment tests with SP1 and RELA added to the TF panel, and positive-control gene pairs. |
| `06_spatial_mediation_summary.py` | (1) Spot-level mediation: TF -> stromal content -> COL1A1/COL3A1/FN1, run inside each Visium section and pooled, mirroring the bulk mediation analysis but with the mediator measured in situ rather than deconvolved. |
