# P1

## Scripts in this folder

| Script | What it does |
|---|---|
| `p1_batch.sh` | Run p1_run_dataset.py over a list file: part<TAB>arms<TAB>gse<TAB>tf<TAB>cell<TAB>tag<TAB>extra-args |
| `p1_knocktf_index.py` | P1a/P3b: parse the KnockTF 2.0 download table (one row per dataset) into knocktf_datasets.tsv. |
| `p1_run_dataset.py` | P1: differential expression for one contrast, from raw data where available (SETTINGS.md). usage: p1_run_dataset.py <part> <arms.tsv> <gse> <tf> <cell> [--route ncbi\|cel\|matrix] [--tag primary\|secondary] [--pseudo] Routes: ncbi … |
| `p1a_arms.py` | P1a: assign each candidate's samples to the perturbed and control arms from the regexes of p1a_triage.tsv (matched against the sample title, then characteristics), and write p1a_arms.tsv (one row per dataset and sample). |
| `p1a_encode_lincs.py` | P1a: (1) ENCODE portal search for knockdown/knockout RNA-seq of ETS1, NFKB1, RELA and SP1 (human); (2) LINCS L1000 landmark status of COL1A1 and COL3A1 (GSE92742 gene table, pr_is_lm). |
| `p1a_screen_col1a1.py` | P1a: COL1A1 expression in the controls of every candidate dataset (SETTINGS change 3: this check may use the NCBI TPM table or the series matrix before raw data are fetched). |
| `p1a_search.py` | P1a search: GEO (the analysis plan's terms per TF) plus the KnockTF 2.0 series for the four TFs, then the first-pass screen. |
| `p1a_summarise.py` | P1a summary: per contrast (P1a/de/*.tsv): TF, COL1A1 and COL3A1 log2FC with SE; knockdown verification (TF log2FC <= -0.74); COL1A1 expression in the controls (RNA-seq TPM >= 10; arrays percentile >= 50); positive control (TRRUST … |
| `p1b_search.py` | P1b search: GEO (the analysis plan's terms: miR-29 mimic / overexpression / transfection / inhibitor / antagomir / knockout; and miR-101 gain for the positive control), Homo sapiens, expression profiling; then the first-pass … |
| `p1b_summarise.py` | P1b summary (SETTINGS.md P1b and change 6). |
| `prep_counts.py` | Build a GSM-labelled count matrix from a submitter's GEO supplementary count file. usage: prep_counts.py <part> <gse> <suppl file name> <out.tsv> [--sep ,] [--title-key bracket\|title\|regex:<pattern>] Columns are matched to GSMs … |
| `ts_human_sites.py` | TargetScan 8.0 Summary Counts (all predictions) -> human rows only (Species ID 9606), one row per gene and miRNA family, written to the cache as ts80_human_family_sites.tsv.gz (gene, family, cons8, cons7m8, cons7a1, cons_total, … |
