# Tier- and sign-aware null models for the three-node FFL test

This analysis runs the three-node FFL motif test (all FFLs, composite, TF-FFL, miRNA-FFL) under several null models:
- the degree-preserving null, the reciprocity-preserving null and the bipartite-layer null;
- three models that start from the reciprocity-preserving null and restrict the swaps of the miRNA -> target layer to its evidence tier, those of
  TF -> target arcs to their annotated sign, or both;
- two waterfall steps beyond the bipartite-layer null, adding the tier split and then the sign split.

Each model is run on the analysed network and on the census graph, with the published TF typing and with the 22
co-regulators retyped as non-TF. The design is in `SETTINGS.md`. The values, with file, line and column, are in
`VALUES.md`.

| file | what it is |
|---|---|
| `motif_null_tier_sign.c` | `motif_null_models.c` with the tier- and sign-split models added; A, B and C are unchanged |
| `make_tier_sign_inputs.py` | Builds the graphs in both typings and the per-arc tier and sign attributes; checks the retyped analysed network against E5 |
| `run_tier_sign_nulls.sh` | Compiles and runs the 32 runs (1,000 replicates each) |
| `summarise_tier_sign_nulls.py` | Checks A, B and C against the published replicate files, then writes the summary and the waterfall |
| `write_values.py` | Writes `VALUES.md` |
| `null_model_summary.csv` | Observed, expected, fold, excess %, empirical p for every run and metric |
| `waterfall.csv` | Excess % from the degree-preserving null to the bipartite-layer null + tier + sign, with the step changes |
| `control_reproduction.csv` | The reproduction of the published replicate files |
| `arc_attribute_counts.csv` | Arcs per tier and per sign |
| `runs/` | The replicate files (one line per randomised network) |

**Run:** `FFL_ROOT=<repository checkout> bash run_tier_sign_nulls.sh`, then `python3 summarise_tier_sign_nulls.py`
and `python3 write_values.py`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `make_tier_sign_inputs.py` | Inputs of the tier- and sign-aware null models: graphs in both typings and one attribute line per arc. |
| `motif_null_tier_sign.c` | Three-node FFL motif test under the null models of motif_null_models.c (analyses/census_and_motif_nulls/nulls/) and under tier- and sign-restricted versions of them. |
| `run_tier_sign_nulls.sh` | Three-node FFL motif test under the degree-preserving null/B/C and the tier- and sign-restricted null models (see SETTINGS.md): 8 models x 2 graphs (analysed network, census graph) x 2 typings (published, 22 co-regulators retyped), 1,000 replicates … |
| `summarise_tier_sign_nulls.py` | Summary of the tier- and sign-aware null models (see SETTINGS.md). |
| `write_values.py` | Writes VALUES.md for the tier- and sign-aware null models from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
