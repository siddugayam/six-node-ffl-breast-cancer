# Tier-stratified miRNA–target sign concordance after stromal adjustment

This analysis takes the sign concordance of network edges by evidence tier and tests it against the five-draw
within-network matched null of `scripts/07_expression_validation/edge_correlations/04_caf_adjusted_vs_matched_null.R`.
Three versions are computed:
- unadjusted;
- adjusted for each of the 45 stromal estimates of the mediation analysis, one at a time (ABSOLUTE purity and the
  ESTIMATE stromal score among them);
- with the targets correlated with the Farmer stroma-related score left out.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `stroma_adjusted_concordance.R` | The analysis. It first reproduces the table of 04, then runs every estimate and edge selection |
| `write_values.py` | Writes `VALUES.md` from the tables |
| `concordance_summary.csv` | Per edge set: unadjusted; ABSOLUTE alone; ESTIMATE alone; median and range over the 45 estimates |
| `concordance_summary_42_distinct.csv` | The same summary over the 42 distinct estimates (three of the 45 runs repeat another estimate exactly); written by `summarise_distinct_estimates.py`, which first checks that the three runs are exact duplicates |
| `concordance_adjusted_by_estimate.csv` | Every estimate × edge selection × edge set, adjusted, and unadjusted on the same tumours |
| `control_reproduction.csv` | The reproduction of 04's published table |
| `run.log` | Output of the run |

**Run:** `FFL_ROOT=<repository checkout> Rscript stroma_adjusted_concordance.R`, then `python3 write_values.py`.

The script needs the TCGA-BRCA expression matrices (`data/brca_*.rds`), which the repository cites but does not
include.

## Scripts in this folder

| Script | What it does |
|---|---|
| `stroma_adjusted_concordance.R` | Sign concordance of network edges by evidence tier, before and after adjustment for each of the 45 stromal estimates of the mediation analysis (ABSOLUTE purity and the ESTIMATE stromal score among them), against the five-draw … |
| `summarise_distinct_estimates.py` | The 45 stromal estimates contain three runs that reproduce other estimates exactly (MCPcounter_Fibroblasts_recomputed = MCPcounter_Fibroblasts, NNLS_Wu641_full_CAFs = NNLS_Wu_CAFs, QPROG_Wu641_full_CAFs = QPROG_Wu_CAFs). |
| `write_values.py` | Writes VALUES.md for the stroma-adjusted tier concordance from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
