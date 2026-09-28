# Tests with public perturbation data (Supplementary Note S1)

The settings were fixed in advance (`SETTINGS.md`, with every later change logged; the parties to the logged decisions
are named "the authors", and the log is otherwise as written). `DATASETS.tsv` lists every dataset found and the reason for
any exclusion; `P3a/p3a_triage_note.txt` explains how the 40 miRNA contrasts were chosen. The per-gene
differential-expression tables are not included; the scripts rebuild them from GEO. Besides public downloads, the
scripts read the network tables, `results/v3/seqreg_ext_occlusion_allsites.csv`, TRRUST v2 and, for `P2/p2_analyse.py`,
the six-node composite instances of `analyses/six_node_pattern/S1`.

| Folder | What it holds |
|---|---|
| `P0/` | The compute resources and the download cap |
| `P1/` | Scripts of P1a (TF knockdowns and the collagens) and P1b (miR-29 gain and loss in fibroblast-lineage and carcinoma cells); `p1_run_dataset.py` runs the differential expression of one contrast |
| `P1a/`, `P1b/` | Their summaries and pooled results |
| `P2/` | Co-transcribed miRNA clusters: whether shared targets are repressed convergently (human; mouse through orthologues) |
| `P3/`, `P3a/` | miRNA→target edges by evidence tier with miRNA gain and loss data (P3a), and TF→target edges with TF knockdowns (P3b) |
| `P4/` | Promoter occlusion with a second sequence model (Borzoi) |
| `P5/` | Mature miR-29 and miR-101 in purified cell types (two atlases); the files with `_as_delivered` are the first computation, before the correction logged in `SETTINGS.md` |
| `_*.py`, `_*.R` | Shared helpers: downloads, GEO and array processing, differential expression and random-effects pooling |
