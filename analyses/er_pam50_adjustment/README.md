# Oestrogen-receptor and PAM50 adjustment

**Part 1.** The six-node circuit survival test of `analyses/six_node_pattern/S8/s8_circuit_survival.R` in TCGA-BRCA
(overall survival; GSVA scores; age and stage; likelihood-ratio test of core vs core + extension; BH within family) is
repeated:
- with oestrogen-receptor status and PAM50 added;
- with age and stage alone on the same tumours.

Results are reported for all 2,383 six-node composite pattern circuits, 71 model-architecture circuits, and the 16
patterns that reached q < 0.05. METABRIC has no miRNA assay, so the circuits cannot be scored there.

**Part 2.** hsa-miR-29a, and the other ten miRNAs of Table S8, are re-fitted adjusted for ER and/or PAM50 in the five
miRNA cohorts that have ER or subtype annotation, then pooled.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `circuit_survival_er_pam50.R` | Part 1. It first reproduces the published S8 table |
| `mir29a_er_pam50.R` | Part 2 |
| `write_values.py` | Writes `VALUES.md` |
| `circuits_er_pam50.csv` | Every circuit × fit: HR of the six-node score (95 % CI), p, C-indices, LRT p and BH q |
| `sixteen_patterns_before_after.csv` | The 16 patterns, before and after |
| `circuit_summary.csv` | Counts per fit and family |
| `mir29a_er_pam50_percohort.csv`, `mir29a_er_pam50_pooled.csv` | Part 2 |
| `control_reproduction.csv` | The reproduction of the published S8 table |
| `run_circuits.log`, `run_mir29a.log` | Run outputs |

**Run:** `FFL_ROOT=<repository checkout> Rscript circuit_survival_er_pam50.R`, `Rscript mir29a_er_pam50.R`, then
`python3 write_values.py`. The scripts need the TCGA-BRCA data and the miRNA cohort tables, which the repository cites
but does not include.

## Scripts in this folder

| Script | What it does |
|---|---|
| `circuit_survival_er_pam50.R` | The six-node circuit survival test of analyses/six_node_pattern/S8/s8_circuit_survival.R (GSVA scores of each six-node set, its three-node core and its extension; Cox models for overall survival adjusted for age and AJCC stage; … |
| `mir29a_er_pam50.R` | hsa-miR-29a and the other ten miRNAs of Table S8: Cox models for each cohort's primary endpoint, unadjusted and adjusted for oestrogen-receptor status and/or PAM50 subtype, in the five miRNA cohorts that have ER or subtype … |
| `write_values.py` | Writes VALUES.md for the oestrogen-receptor and PAM50 adjustment from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
