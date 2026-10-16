# Stromal mediation on a common sample of tumours

This analysis repeats the mediation of the twelve paths of 46 for each of the 42 distinct stromal estimates: ETS1, NFKB1
and hsa-miR-29a/b/c with COL1A1 and COL3A1, and SP1 and RELA with COL1A1. Seven of them are the rows of Table 3. All runs use one common sample of tumours (both assays and every
estimate available), so each path's total effect is the same for every estimate. The model is `run1()` of
`scripts/08_stroma_and_cell_types/deconvolution_and_mediation/46_deconv_extended_mediation.R`, unchanged.

The design is in `SETTINGS.md`. The values, with file, line and column, are in `VALUES.md`.

| file | what it is |
|---|---|
| `mediation_common_sample.R` | Reproduces 46's published rows first, then runs the common sample |
| `write_values.py` | Writes `VALUES.md` |
| `mediation_common_sample_summary.csv` | Per path: total, ACME/ADE/proportion mediated (median and range), inconsistent counts, beside the published runs; Table 3 rows flagged |
| `mediation_common_sample_per_estimate.csv` | Every path × estimate |
| `control_reproduction.csv` | The reproduction of the published rows |
| `run.log` | Output of the run |

**Run:** `FFL_ROOT=<repository checkout> Rscript mediation_common_sample.R`, then `python3 write_values.py`. The script
reads the stromal estimate matrix that 46 reads (`cache/systems_and_dynamics/deconv_extended/fibroblast_estimates_extended.rds`).
The repository's CSV export of that matrix is rounded to 15 significant digits, which changes the ranks of the
MCPcounter columns.

## Scripts in this folder

| Script | What it does |
|---|---|
| `mediation_common_sample.R` | The stromal mediation of the twelve paths of 46 (ETS1, NFKB1 and hsa-miR-29a/b/c, each with COL1A1 and COL3A1; SP1 and RELA with COL1A1; Table 3 shows seven of them), repeated for each of the 42 distinct stromal estimates on one … |
| `write_values.py` | Writes VALUES.md for the common-sample mediation from the tables in this folder: every number with its file, line (header = line 1) and column. usage (from this folder): python3 write_values.py |
