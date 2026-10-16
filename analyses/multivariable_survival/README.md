# Multivariable survival for the miRNAs of the survival meta-analysis

Unadjusted and covariate-adjusted Cox models for the 11 miRNAs of Table S8, in the six miRNA cohorts of the
meta-analysis. The folder holds:
- Schoenfeld tests;
- fixed- and random-effects pooling of the adjusted HRs, with BH over the 11 miRNAs;
- an overall-survival-only pool;
- a pool without the stromal score.

The design is in `SETTINGS.md` and the values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `mirna_multivariable_cox.R` | All models and pools; it first checks that the unadjusted pool reproduces Table S8 |
| `write_values.py` | Writes `VALUES.md` from the tables |
| `percohort_cox.csv` | Every Cox fit: cohort, miRNA, endpoint, model, covariates, n, events, HR, CI, p, Schoenfeld p |
| `pooled_adjusted.csv` | Pooled unadjusted and adjusted HRs, primary endpoint (fixed and random effects, I², BH q) |
| `pooled_os_harmonised.csv` | The same on overall survival only (three cohorts) |
| `pooled_adjusted_without_stromal.csv` | The adjusted pool with the stromal score left out |
| `table_s8_reproduction.csv` | The unadjusted pool against Table S8 |
| `cohort_covariates.csv` | Number of patients with each covariate per cohort |
| `run.log` | Output of the run |

**Run:** `FFL_ROOT=<repository checkout> Rscript mirna_multivariable_cox.R`, then `python3 write_values.py`.

The script also reads three inputs that the repository cites but does not include: the TCGA-BRCA clinical matrix, the
miRNA cohort tables under `cache/`, and the GEO mRNA matrices. These are the inputs of
`scripts/10_survival_and_clinical/mirna_survival_meta_analysis/c11_mirna_cox_meta.R`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `mirna_multivariable_cox.R` | Unadjusted and covariate-adjusted Cox models for the 11 miRNAs of the survival meta-analysis (Table S8), Schoenfeld tests, fixed- and random-effects pooling, BH over the 11 miRNAs, the harmonised overall-survival sensitivity and … |
| `write_values.py` | Writes VALUES.md for the multivariable survival analysis from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
