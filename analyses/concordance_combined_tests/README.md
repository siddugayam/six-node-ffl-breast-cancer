# Calibrated combination of the stromal-adjusted concordance tests

The stromal-adjusted concordance tests average the two-sided p values over five matched-null draws. Here each test is
also computed in four other ways:
- (b) one test against the pooled null pairs of all five draws;
- (c) the Cauchy and the harmonic-mean (Vovk–Wang) combinations of the five p values;
- (d) one-sided, as in the unadjusted test;
- and twice the mean, a valid p value under any dependence.

These are given for the unadjusted test, CAF-A (the score behind the 58.1 % / 56.0 % TF -> target activation test),
CAF-B, ABSOLUTE purity, the ESTIMATE stromal score and each of the 42 distinct stromal estimates. Each covers every
evidence tier, with all edges and with the Farmer-correlated targets left out.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `combined_concordance_tests.R` | The analysis. It reproduces the published table of `04_caf_adjusted_vs_matched_null.R` first |
| `write_values.py` | Writes `VALUES.md` |
| `combined_tests_by_adjustment.csv` | Every adjustment × edge set × edge selection: concordance, null, all combined p values, per-draw p values |
| `combined_tests_summary_42_distinct.csv` | Summary over the 42 distinct estimates, in the layout of `concordance_summary_42_distinct.csv` with the combined-test columns added |
| `conclusion_changes.csv` | Significance at p < 0.05 by method, and the methods that change the call from the arithmetic mean |
| `control_reproduction.csv` | The reproduction of the published table |
| `run.log` | Output of the run |

**Run:** `FFL_ROOT=<repository checkout> Rscript combined_concordance_tests.R`, then `python3 write_values.py`. It needs
the TCGA-BRCA expression matrices (`data/brca_*.rds`), which the repository cites but does not include.

## Scripts in this folder

| Script | What it does |
|---|---|
| `combined_concordance_tests.R` | Sign concordance of network edges by evidence tier against the five-draw within-network matched null, combined over the five draws in several ways (arithmetic mean of the two-sided p values; one test against the pooled null pairs … |
| `write_values.py` | Writes VALUES.md for the combined concordance tests from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
