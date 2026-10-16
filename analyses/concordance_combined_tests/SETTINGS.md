# SETTINGS: calibrated combination of the stromal-adjusted concordance tests

Settings of the analysis.

## Base design (unchanged)

- **The script.** `scripts/07_expression_validation/edge_correlations/04_caf_adjusted_vs_matched_null.R`, as extended in
  the stroma-adjusted concordance analysis (`stroma_adjusted_concordance.R`). Its data preparation and null draws are
  copied unchanged.
- **Samples.** TCGA-BRCA primary tumours with both assays, n = 1,066.
- **Correlation.** Spearman, as rank-z dot products. Adjusted correlations are first-order partial Spearman.
- **Matched null.** For every edge, 10 pairs of network nodes per draw (K = 10) of the same kind (gene or miRNA) and the
  same mean-expression decile; self-pairs and real edges excluded. Five draws, seeds 501–505.
- **Concordance.** The fraction of edges whose correlation sign equals the predicted sign.
- **Control.** The published table of 04 is reproduced first (relative tolerance 1e-9).

## Adjustments covered

- unadjusted;
- CAF-A, the stromal score behind the 58.1 % / 56.0 % row: the ESTIMATE stromal signature without collagens and
  network nodes, 123 genes, mean z, from `results/expression_validation/caf_scores.rds`;
- CAF-B (DCN, LUM, FAP, THY1);
- ABSOLUTE purity alone (1 − purity; 1,023 tumours);
- the ESTIMATE stromal score alone;
- each of the 42 distinct stromal estimates of the mediation analysis.
  - Of the 45 columns of `deconv_fibroblast_estimates_extended.csv`, three repeat another exactly and are counted once
    (dropped): MCPcounter_Fibroblasts_recomputed, NNLS_Wu641_full_CAFs and QPROG_Wu641_full_CAFs.
  - That they are exact repeats is checked by the script.

## Edge sets

- TF -> target activation (the 58.1 % row), TF -> target repression, TF -> miRNA activation and repression;
- miRNA -> target: all, strong, weak, prediction-only;
- each with all edges, and with the Farmer-correlated targets left out (rho > 0.5 with the Farmer stroma-related score,
  as in the stroma-adjusted concordance analysis).

## Combined tests per adjustment × edge set × edge selection

Let p1…p5 be the per-draw p values, each an exact binomial test of x concordant of n edges against that draw's null
concordance p0_d (10 n pairs).

- **(a) Arithmetic mean** of the five two-sided p values (the published method).
  - Also shown: twice the mean, capped at 1. This is a valid p value under any dependence (Rüschendorf; Vovk and Wang
    2020).
- **(b) Pooled reference distribution.** One exact two-sided binomial test of x of n against the concordance of all 50 n
  null pairs of the five draws pooled.
- **(c) Combinations of the five two-sided p values:**
  - the Cauchy combination (Liu and Xie 2020): T = mean of tan((0.5 − p_i)π), p = 0.5 − arctan(T)/π;
  - the harmonic mean with the Vovk–Wang factor for arbitrary dependence: p = min(1, e · ln 5 · HM(p1…p5)).
- **(d) One-sided, matching the unadjusted test of Fig. S6** (alternative: concordance greater than the null):
  - the pooled-reference one-sided exact binomial (primary for (d));
  - the mean and the Cauchy combination of the five one-sided p values.

**Conclusions.** Significance is read at p < 0.05 for each method. The table shows, for each adjustment × set, whether
any method changes the call relative to (a).

## Summaries

- The 42 distinct estimates are summarised by: the median and range of p per method; the number of estimates with
  p < 0.05 per method; the number above the null.
- **Layout.** The `concordance_summary_42_distinct.csv` of `analyses/stroma_adjusted_concordance/`: edges, set, n_estimates, p_median, p_min, p_max,
  n_p_below_0.05, n_above_null, excess_median, excess_min, excess_max. Here p_* are the arithmetic-mean protocol (a).
  - The combined-test columns are added after them, for each method: <method>_p_median, _p_min, _p_max and
    _n_below_0.05.
  - The reference row (all edges, miRNA_target_weak: 42; 0.05152…; 20; 40) must be reproduced.

## Software

R (stats, matrixStats). Seeds 501–505 are from 04. No other random numbers.

## Outputs

- `combined_tests_by_adjustment.csv`: every adjustment × set × edge selection × method.
- `combined_tests_summary_42_distinct.csv`
- `conclusion_changes.csv`
- `control_reproduction.csv`
