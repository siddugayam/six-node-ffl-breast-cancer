# SETTINGS: stromal mediation on a common sample of tumours

Settings of the analysis.

## Model (unchanged)

The function `run1()` of `scripts/08_stroma_and_cell_types/deconvolution_and_mediation/46_deconv_extended_mediation.R`,
copied unchanged:
- normal-score transforms of exposure, mediator and outcome; m ~ x and y ~ x + m;
- ACME = a·b, ADE = c′, total = ACME + ADE;
- percentile-bootstrap CIs from 2,000 resamples, with seed 11 for every model as in 46;
- Sobel test;
- inconsistent = ACME and ADE of opposite sign.

## Paths

All twelve paths of 46: ETS1, NFKB1, hsa-miR-29a, hsa-miR-29b and hsa-miR-29c, each with COL1A1 and COL3A1; SP1 and
RELA with COL1A1.

Table 3 shows seven of them: ETS1, NFKB1, SP1 and RELA -> COL1A1; hsa-miR-29a -> COL3A1; hsa-miR-29a -> COL1A1;
hsa-miR-29b -> COL3A1. The published inconsistent counts of those seven rows (23, 32, 38, 22, 6, 6, 3) are checked by
the script.

## Stromal estimates

- The 45 columns of `cache/systems_and_dynamics/deconv_extended/fibroblast_estimates_extended.rds`, the mediators of 46.
- Three repeat another column exactly and are counted once: MCPcounter_Fibroblasts_recomputed,
  NNLS_Wu641_full_CAFs and QPROG_Wu641_full_CAFs are left out. The script checks that each equals a kept column.
- That leaves 42 distinct estimates.

## Common sample

- **Who is in it.** Every TCGA-BRCA primary tumour with both the gene and the miRNA assay (the miR-29 paths need
  both), finite values of all 42 estimates, and every exposure and both outcomes (COL1A1, COL3A1).
- **One sample for all twelve paths**, so the total effect of a path is identical across the 42 estimates; the script
  checks this.
- **n.** The script reports n. It is at most 767 (the methylation-based EpiDISH estimates).

## Control

The same function on each estimate's own available tumours reproduces the published rows of
`results/systems_and_dynamics/deconv_mediation_extended.csv` for the twelve paths. Checked: n, total, ACME, ADE and their
CIs, relative tolerance 1e-9.

## Reported

- **Per path × estimate on the common sample:** n, total, ACME (95 % CI), ADE (95 % CI), proportion mediated (95 % CI),
  Sobel p, inconsistent.
- **Per path,** beside the published counts on all available tumours:
  - the total effect;
  - the median and range over the 42 of ACME, ADE and proportion mediated;
  - the number of inconsistent estimates;
  - the number of inconsistent estimates whose ADE CI excludes 0;
  - the number with Sobel p < 0.05.

## Software

R (data.table, parallel). The bootstrap seed is 11 per model, as in 46. Run with 4 cores; mclapply does not change any
model's random stream, because each model sets its own seed.

## Outputs

- `mediation_common_sample_per_estimate.csv`
- `mediation_common_sample_summary.csv`
- `control_reproduction.csv`

## Input of the estimates

- The script reads the matrix of `cache/systems_and_dynamics/deconv_extended/fibroblast_estimates_extended.rds`. The file
  `results/systems_and_dynamics/deconv_fibroblast_estimates_extended.csv` holds the same values at 17 significant digits.
