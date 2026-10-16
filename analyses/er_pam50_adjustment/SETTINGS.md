# SETTINGS: oestrogen-receptor and PAM50 adjustment of the six-node circuit survival test and of hsa-miR-29a

Settings of the analysis.

## Part 1: six-node circuits in TCGA-BRCA

**The published test** (`analyses/six_node_pattern/S8/s8_circuit_survival.R`, copied unchanged for data and scores):
- **Circuits.** The 2,498 node sets of the six-node composite pattern and the 71 sets of the six-node model architecture;
  2,383 and 71 have all six members measured.
- **Data.** TCGA-BRCA primary tumours with gene and miRNA expression and survival.
- **Scores.** GSVA (kcdf Gaussian, minSize 2), z-scaled, for the six-node set, its three-node core and its extension.
- **Model.** Overall survival; Cox with age and AJCC stage group.
  - The six-node model is ~ six + age + stage; the core model ~ core + age + stage; the nested model ~ core + extension
    + age + stage.
  - Test: likelihood-ratio test of core vs nested, with BH within family.
- **The 16 patterns** are the six-node composite circuits with BH q < 0.05 in this test.

**This analysis** runs on the same endpoint (overall survival) under the same multiple-testing rule (BH over the circuits
of each family).
- **Covariates added.** ER status from `breast_carcinoma_estrogen_receptor_status` of `data/brca_clinicalMatrix.tsv`
  (Positive / Negative; other values missing) and PAM50 from `PAM50Call_RNAseq` (LumA reference; Normal kept).
- **Three fits per circuit:**
  1. **published:** age + stage on the published sample. Its HR, p, LRT p and q must reproduce
     `six_node_pattern_survival_all_circuits.csv` (OS rows, relative tolerance 1e-6) before anything else is reported.
  2. **same tumours:** age + stage on the tumours that have ER and PAM50 (complete cases), to separate the effect of
     the smaller sample;
  3. **ER and PAM50:** age + stage + ER + PAM50 on those tumours.
- **Reported per circuit and fit:**
  - HR of the six-node score per SD, its 95 % Wald CI and p;
  - the LRT p (extension adds to the core) and its BH q.
- **For the 16 patterns:** before and after, with the number that keep nominal significance (LRT p < 0.05; and,
  separately, HR p < 0.05) and the number that keep q < 0.05.

## Part 2: hsa-miR-29a (and the other ten miRNAs of Table S8) adjusted for ER and/or PAM50

- **Method.** The cohort loading, `pick()` and z-scoring of `mirna_multivariable_cox.R` (from the multivariable survival
  analysis), with each cohort's primary endpoint as in Table S8.
- **Cohorts and covariates:**

  | cohort | covariates |
  |---|---|
  | TCGA-BRCA | ER (IHC) + PAM50 |
  | GSE19783 | ER + PAM50 subtype |
  | GSE22216 | ER |
  | GSE59829 | ER (esr1_status) |
  | GSE78870 | ER |

  GSE37405 is ER-positive by design and has no subtype, so it is not part of this pool.
- **Models.** Unadjusted and adjusted, both on the complete cases of the adjusted model.
- **Pooling.** Fixed effect and random effects (DerSimonian–Laird) over the five cohorts, with I² and BH over the 11
  miRNAs.
- **Schoenfeld test** per fit (cox.zph, transform "km").

## Software

R (survival, GSVA, metafor, data.table). The only seed is set.seed(1234), copied from S8; GSVA with a Gaussian kernel is
deterministic.

## Outputs

- `circuits_er_pam50.csv`: one row per circuit × fit.
- `sixteen_patterns_before_after.csv`
- `circuit_summary.csv`
- `mir29a_er_pam50_percohort.csv`
- `mir29a_er_pam50_pooled.csv`
- `control_reproduction.csv`
