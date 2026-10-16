# SETTINGS: multivariable survival for the miRNAs of the survival meta-analysis (Table S8)

The design of the analysis.

## miRNAs, expression values and endpoints

**miRNAs.** The 11 of Table S8: miR-21, miR-195, miR-204, miR-383, miR-124, miR-155, miR-429, miR-141, miR-34a,
miR-101 and miR-29a.
- The row per miRNA and cohort is chosen with `pick()` of
  `scripts/10_survival_and_clinical/mirna_survival_meta_analysis/c11_mirna_cox_meta.R`, unchanged.
- Each is z-scored within the cohort (within platform for GSE37405), so HR is per SD as in Table S8.

**Cohorts and endpoints.** As c11; each cohort's primary endpoint is used.

| cohort | endpoint |
|---|---|
| TCGA-BRCA | OS |
| GSE19783 | OS (the deposited time variable is the disease-free time, as c11 notes) |
| GSE22216 | DRFS |
| GSE37405 | RFS, with platform as a stratum |
| GSE59829 | DMFS |
| GSE78870 | TTP |

**METABRIC is not analysed.** Its public release (cBioPortal brca_metabric) has no miRNA profiles; the miRNA data are
controlled-access at EGA.

## Covariates per cohort ("maximal available" model)

| cohort | covariates |
|---|---|
| TCGA-BRCA | age; stage; ER; PAM50; stromal score |
| GSE19783 | ER; PAM50 (subtype); stromal score |
| GSE22216 | age; grade; size; node; ER; stromal score |
| GSE37405 | age; size; node (an ER-positive cohort by design; platform stratum) |
| GSE59829 | age; size; ER |
| GSE78870 | age; grade; ER |

Coding:
- **TCGA-BRCA:** age continuous.
  - Stage I / II / III / IV from AJCC (`stage` of the cohort table). Stage X and missing values are missing.
  - ER from `breast_carcinoma_estrogen_receptor_status` of `data/brca_clinicalMatrix.tsv` (Positive / Negative; other
    values missing).
  - PAM50 from `PAM50Call_RNAseq` (LumA as reference; Normal kept as a level).
  - Stromal score as below.
- **GSE19783:** ER (Positive / Negative); subtype (Lum A reference; "not classified" missing); stromal score.
- **GSE22216:** age, size and node as stored; grade 1 / 2 / 3 as a factor; ER 0 / 1; stromal score.
- **GSE37405:** age, size and node as in c11.
- **GSE59829:** age; size; ER from esr1_status.
- **GSE78870:** age (`age at start therapy`, as c11); grade 1 / 2 / 3 as a factor; ER from `estrogen receptor status`.
- No cohort has stage except TCGA. Grade, size and node are used where the cohort has them.

**Stromal score.** CAF-A, the stromal covariate of the single-covariate mediation analysis
(`scripts/08_stroma_and_cell_types/deconvolution_and_mediation/05_mediation.R`). It is the mean z over the 123 genes of
the ESTIMATE stromal signature, after removing every collagen and every network node.
- TCGA-BRCA: `meanZ_CAF_A_estimate` of `results/systems_and_dynamics/deconv_fibroblast_estimates_extended.csv`.
- GSE19783: the matched GPL6480 mRNA (the gene-level matrix in `gse19783.rds`).
- GSE22216: the matched GSE22219 mRNA (gene-level matrix of the cohort table, matched with
  `GSE22219_GSE22216_link.csv`).
- In both GEO cohorts it is computed with `score_of()` of `scripts/07_expression_validation/external_cohorts/N2_cohort_analysis.R` (the mean of row z over the
  genes present).
- The other cohorts have no matched mRNA, so they have no stromal score.

## Models

All models are Cox proportional hazards (survival::coxph), on complete cases with time > 0. A fit is made only if
n ≥ 25 and there are ≥ 8 events (the c11 rule).

For each cohort × miRNA:
1. **unadjusted**: miRNA z alone (platform stratum for GSE37405);
2. **adjusted**: miRNA z + the cohort's maximal covariate set above.

TCGA's adjusted model is the full model (age, stage, ER, PAM50, stromal score).

- **Reported per fit:**
  - HR per SD, Wald 95 % CI and p;
  - n and events;
  - the covariates used;
  - the Schoenfeld test (survival::cox.zph, transform "km"): the global p and the p of the miRNA term.
- **Pooling.** The adjusted log-HRs are pooled over the cohorts with an adjusted fit:
  - fixed effect (inverse variance) and random effects (DerSimonian–Laird, `metafor::rma`, as in c11);
  - with k, I², tau², Q and its p.
  - BH q over the 11 miRNAs is computed separately for the fixed-effect and the random-effects p.
- **The unadjusted primary-endpoint pool is recomputed in the same script.** It must reproduce Table S8 (HR, CI, p
  per miRNA, at the published rounding) before anything else is reported.
- **Harmonised endpoint (sensitivity).**
  - OS where a cohort has it: TCGA-BRCA OS, GSE19783 OS (deposited time variable, as above) and GSE37405 OS.
  - This drops GSE22216 (DRFS only), GSE59829 (DMFS only) and GSE78870 (TTP/DFI only).
  - Both unadjusted and adjusted models are pooled, fixed and random, with I² and BH over the 11 miRNAs.
- **No-stromal-score sensitivity.** The same adjusted pool, with the stromal score left out in the three cohorts that have
  it. This shows how much of any change comes from the stromal score.

**No decision rule.** The analysis reports HRs. "Remains associated" is read off as pooled adjusted p < 0.05 with
BH q < 0.05, and is stated as such whichever way it falls.

## Software

R with survival, metafor and data.table. No random numbers are used.

## Outputs

- `percohort_cox.csv`: every fit, with its covariates, n, events, HR, CI, p and the Schoenfeld p.
- `pooled_adjusted.csv`
- `pooled_os_harmonised.csv`
- `pooled_adjusted_without_stromal.csv`
- `table_s8_reproduction.csv`
- `cohort_covariates.csv`: n available per covariate per cohort.

No sample-level values are written.

