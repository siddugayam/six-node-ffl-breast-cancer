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

## Scripts in this folder

| Script | What it does |
|---|---|
| `_agilent.R` | Agilent Feature Extraction files -> normalised matrices (SETTINGS.md: "normexp + quantile for Agilent"). usage: Rscript _agilent.R single <out.tsv> <files...> single-channel: normexp (offset 16), quantile; log2 E Rscript … |
| `_annodb.R` | Id -> SYMBOL from a Bioconductor annotation package (used when the GPL table has no gene symbols, e.g. |
| `_arr.py` | Array helpers for analyses/perturbation_tests: GPL probe -> gene symbol, and raw CEL download plus RMA (via _cel.R). |
| `_cel.R` | RMA of Affymetrix CEL files with oligo (SETTINGS.md: "RMA for Affymetrix CEL"). |
| `_datasets.py` | DATASETS.tsv: every dataset found, per part, with its status and reason (SETTINGS.md, "Every hit is listed"). |
| `_de2.R` | One-sample limma on two-colour log ratios (perturbed / control), for paired two-colour designs. usage: Rscript _de2.R <spec.json> <out.tsv> spec: {"M": path, "A": path or null, "gene_map": path, "orient": {"<GSM>": 1 \| -1}, … |
| `_de.R` | Differential expression for one dataset of analyses/perturbation_tests (SETTINGS.md, "Differential expression"). usage: Rscript _de.R <spec.json> <out.tsv> spec: {"kind": "counts" \| "log_matrix" \| "linear_matrix", "matrix": path … |
| `_design.py` | Group a series' samples by their title with replicate numbers removed, to read the design quickly. usage: _design.py <part> <regex\|-> <GSE> [<GSE> ...] (groups not matching the regex are only counted) |
| `_dl.py` | Download helper for analyses/perturbation_tests. |
| `_geo.py` | GEO helpers for analyses/perturbation_tests. |
| `_mouse.py` | Mouse -> human gene mapping by one-to-one orthologues (MGI HOM_MouseHumanSequence.rpt: a homology class with exactly one mouse and one human gene), as in P2/p2_mouse.py. |
| `_rma.R` | Random-effects pooling (metafor REML) and meta-regression for analyses/perturbation_tests. usage: Rscript _rma.R <in.tsv> <out.tsv> in.tsv columns: group, yi, sei [, mod] (one row per dataset). |
| `_screen.py` | Search and first-pass screen shared by P1a, P1b, P2 and P3a. search(part, queries): runs each GEO query (db=gds, GSE, Homo sapiens unless stated, expression profiling), logs the query, date and count to … |
| `_seed.py` | Per-dataset bootstrap seeds: rng_for(key) is a generator seeded from 20250908 and the key (dataset + test), so a bootstrap SE does not depend on which other datasets were processed before it. |
