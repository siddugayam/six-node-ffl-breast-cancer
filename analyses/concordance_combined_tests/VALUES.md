# VALUES: calibrated combination of the stromal-adjusted concordance tests

Lines count the header as line 1. Tables are written by `combined_concordance_tests.R` (run log `run.log`).

- **Matched pairs per edge per draw:** 10, five draws, seeds 501–505 (`combined_tests_by_adjustment.csv`, column
  pairs_per_edge_per_draw).
- **The 58.1 % and 56.0 % of the TF -> target activation test** come from the CAF-A stromal score: the ESTIMATE stromal
  signature without collagens and network nodes, 123 genes, mean z. They are the CAF_A row below.
- **Columns:**
  - conc, null_pooled;
  - p_mean_two_sided (a), p_twice_mean_two_sided (2 × a, valid under any dependence);
  - p_pooled_two_sided (b);
  - p_cauchy_two_sided, p_hmp_vovk_wang_two_sided (c);
  - p_pooled_one_sided (d), p_mean_one_sided, p_cauchy_one_sided.
- **One-sided:** the alternative is concordance greater than the null, as in the unadjusted test of Fig. S6.

## 1. Controls (`control_reproduction.csv`)

The published table of 04 is reproduced. The unadjusted, CAF-A and CAF-B rows equal its conc and p columns; the script
asserts this.

## 2. Key tests (`combined_tests_by_adjustment.csv`)

| line | adjustment | edges | set | n | conc | null (pooled) | (a) mean | 2 × mean | (b) pooled | (c) Cauchy | (c) HM Vovk–Wang | (d) pooled one-sided | mean one-sided | Cauchy one-sided |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | unadjusted | all edges | TF_target_Activation | 301 | 0.7110 | 0.6352 | 0.00705 | 0.0141 | 0.00587 | 0.00577 | 0.0253 | 0.00337 | 0.00367 | 0.00299 |
| 6 | unadjusted | all edges | miRNA_target_all | 4628 | 0.4786 | 0.4735 | 0.499 | 0.998 | 0.489 | 0.515 | 1 | 0.248 | 0.252 | 0.241 |
| 7 | unadjusted | all edges | miRNA_target_strong | 701 | 0.4950 | 0.4882 | 0.73 | 1 | 0.734 | 0.783 | 1 | 0.373 | 0.374 | 0.371 |
| 8 | unadjusted | all edges | miRNA_target_weak | 2294 | 0.5266 | 0.4850 | 0.000155 | 0.00031 | 7.17e-05 | 3.06e-05 | 0.000134 | 3.72e-05 | 8.03e-05 | 1.63e-05 |
| 9 | unadjusted | all edges | miRNA_target_predicted_only | 1633 | 0.4042 | 0.4510 | 0.000183 | 0.000367 | 0.000141 | 0.000104 | 0.000454 | 1 | 1 | 1 |
| 10 | unadjusted | Farmer-correlated targets left out | TF_target_Activation | 197 | 0.7107 | 0.6350 | 0.0307 | 0.0614 | 0.0316 | 0.0278 | 0.122 | 0.0154 | 0.0161 | 0.0146 |
| 14 | unadjusted | Farmer-correlated targets left out | miRNA_target_all | 3536 | 0.4542 | 0.4691 | 0.0838 | 0.168 | 0.0769 | 0.0704 | 0.307 | 0.964 | 0.96 | 0.967 |
| 15 | unadjusted | Farmer-correlated targets left out | miRNA_target_strong | 496 | 0.4456 | 0.4808 | 0.125 | 0.251 | 0.126 | 0.113 | 0.493 | 0.947 | 0.944 | 0.951 |
| 16 | unadjusted | Farmer-correlated targets left out | miRNA_target_weak | 1780 | 0.5017 | 0.4823 | 0.119 | 0.238 | 0.102 | 0.0791 | 0.342 | 0.0532 | 0.0611 | 0.04 |
| 17 | unadjusted | Farmer-correlated targets left out | miRNA_target_predicted_only | 1260 | 0.3905 | 0.4459 | 8.44e-05 | 0.000169 | 7.19e-05 | 6.87e-05 | 0.0003 | 1 | 1 | 1 |
| 18 | CAF-A | all edges | TF_target_Activation | 301 | 0.5814 | 0.5603 | 0.491 | 0.981 | 0.486 | 0.491 | 1 | 0.249 | 0.252 | 0.245 |
| 22 | CAF-A | all edges | miRNA_target_all | 4628 | 0.4697 | 0.4642 | 0.472 | 0.945 | 0.452 | 0.507 | 1 | 0.229 | 0.238 | 0.214 |
| 23 | CAF-A | all edges | miRNA_target_strong | 701 | 0.5107 | 0.4724 | 0.0493 | 0.0986 | 0.0449 | 0.0378 | 0.165 | 0.0231 | 0.026 | 0.0198 |
| 24 | CAF-A | all edges | miRNA_target_weak | 2294 | 0.4843 | 0.4715 | 0.256 | 0.513 | 0.225 | 0.195 | 0.797 | 0.113 | 0.13 | 0.0938 |
| 25 | CAF-A | all edges | miRNA_target_predicted_only | 1633 | 0.4317 | 0.4505 | 0.133 | 0.267 | 0.129 | 0.126 | 0.551 | 0.94 | 0.938 | 0.942 |
| 26 | CAF-A | Farmer-correlated targets left out | TF_target_Activation | 197 | 0.6142 | 0.5552 | 0.103 | 0.207 | 0.0995 | 0.101 | 0.44 | 0.0549 | 0.0557 | 0.0537 |
| 30 | CAF-A | Farmer-correlated targets left out | miRNA_target_all | 3536 | 0.4505 | 0.4590 | 0.333 | 0.667 | 0.319 | 0.302 | 1 | 0.847 | 0.84 | 0.861 |
| 31 | CAF-A | Farmer-correlated targets left out | miRNA_target_strong | 496 | 0.4738 | 0.4636 | 0.664 | 1 | 0.653 | 0.697 | 1 | 0.34 | 0.343 | 0.336 |
| 32 | CAF-A | Farmer-correlated targets left out | miRNA_target_weak | 1780 | 0.4691 | 0.4677 | 0.801 | 1 | 0.905 | 0.942 | 1 | 0.462 | 0.463 | 0.459 |
| 33 | CAF-A | Farmer-correlated targets left out | miRNA_target_predicted_only | 1260 | 0.4151 | 0.4448 | 0.0356 | 0.0712 | 0.0359 | 0.0329 | 0.144 | 0.984 | 0.984 | 0.985 |
| 34 | CAF-B | all edges | TF_target_Activation | 301 | 0.6545 | 0.6006 | 0.0668 | 0.134 | 0.0596 | 0.0511 | 0.223 | 0.0314 | 0.0351 | 0.0269 |
| 38 | CAF-B | all edges | miRNA_target_all | 4628 | 0.4475 | 0.4526 | 0.51 | 1 | 0.488 | 0.493 | 1 | 0.761 | 0.753 | 0.782 |
| 39 | CAF-B | all edges | miRNA_target_strong | 701 | 0.4679 | 0.4647 | 0.874 | 1 | 0.88 | 1 | 1 | 0.447 | 0.447 | 0.446 |
| 40 | CAF-B | all edges | miRNA_target_weak | 2294 | 0.4690 | 0.4606 | 0.447 | 0.895 | 0.426 | 0.541 | 1 | 0.213 | 0.226 | 0.192 |
| 41 | CAF-B | all edges | miRNA_target_predicted_only | 1633 | 0.4085 | 0.4362 | 0.0256 | 0.0511 | 0.0247 | 0.023 | 0.101 | 0.989 | 0.988 | 0.99 |
| 42 | CAF-B | Farmer-correlated targets left out | TF_target_Activation | 197 | 0.6954 | 0.6054 | 0.0124 | 0.0249 | 0.0106 | 0.00815 | 0.0356 | 0.00542 | 0.0066 | 0.00438 |
| 46 | CAF-B | Farmer-correlated targets left out | miRNA_target_all | 3536 | 0.4403 | 0.4477 | 0.399 | 0.798 | 0.388 | 0.381 | 1 | 0.814 | 0.808 | 0.825 |
| 47 | CAF-B | Farmer-correlated targets left out | miRNA_target_strong | 496 | 0.4556 | 0.4558 | 0.879 | 1 | 1 | 0.913 | 1 | 0.52 | 0.519 | 0.52 |
| 48 | CAF-B | Farmer-correlated targets left out | miRNA_target_weak | 1780 | 0.4663 | 0.4580 | 0.508 | 1 | 0.49 | 0.523 | 1 | 0.249 | 0.257 | 0.231 |
| 49 | CAF-B | Farmer-correlated targets left out | miRNA_target_predicted_only | 1260 | 0.3976 | 0.4298 | 0.0215 | 0.043 | 0.0212 | 0.0208 | 0.0908 | 0.991 | 0.99 | 0.991 |
| 50 | ESTIMATE stromal score | all edges | TF_target_Activation | 301 | 0.6047 | 0.5668 | 0.206 | 0.413 | 0.201 | 0.193 | 0.836 | 0.102 | 0.105 | 0.097 |
| 54 | ESTIMATE stromal score | all edges | miRNA_target_all | 4628 | 0.4477 | 0.4511 | 0.664 | 1 | 0.647 | 0.754 | 1 | 0.686 | 0.676 | 0.71 |
| 55 | ESTIMATE stromal score | all edges | miRNA_target_strong | 701 | 0.4750 | 0.4608 | 0.473 | 0.946 | 0.449 | 0.456 | 1 | 0.236 | 0.244 | 0.217 |
| 56 | ESTIMATE stromal score | all edges | miRNA_target_weak | 2294 | 0.4656 | 0.4578 | 0.484 | 0.969 | 0.463 | 0.622 | 1 | 0.235 | 0.245 | 0.217 |
| 57 | ESTIMATE stromal score | all edges | miRNA_target_predicted_only | 1633 | 0.4109 | 0.4376 | 0.0311 | 0.0622 | 0.03 | 0.0294 | 0.128 | 0.986 | 0.986 | 0.987 |
| 58 | ESTIMATE stromal score | Farmer-correlated targets left out | TF_target_Activation | 197 | 0.6548 | 0.5677 | 0.0157 | 0.0315 | 0.0143 | 0.013 | 0.0571 | 0.00783 | 0.00851 | 0.00716 |
| 62 | ESTIMATE stromal score | Farmer-correlated targets left out | miRNA_target_all | 3536 | 0.4400 | 0.4460 | 0.496 | 0.992 | 0.477 | 0.494 | 1 | 0.769 | 0.76 | 0.786 |
| 63 | ESTIMATE stromal score | Farmer-correlated targets left out | miRNA_target_strong | 496 | 0.4556 | 0.4524 | 0.858 | 1 | 0.892 | 1 | 1 | 0.46 | 0.46 | 0.459 |
| 64 | ESTIMATE stromal score | Farmer-correlated targets left out | miRNA_target_weak | 1780 | 0.4601 | 0.4545 | 0.654 | 1 | 0.634 | 1 | 1 | 0.326 | 0.332 | 0.313 |
| 65 | ESTIMATE stromal score | Farmer-correlated targets left out | miRNA_target_predicted_only | 1260 | 0.4056 | 0.4316 | 0.0644 | 0.129 | 0.0645 | 0.0629 | 0.275 | 0.971 | 0.971 | 0.971 |
| 626 | ABSOLUTE purity | all edges | TF_target_Activation | 301 | 0.6844 | 0.6070 | 0.00656 | 0.0131 | 0.00654 | 0.00532 | 0.0233 | 0.00326 | 0.00364 | 0.00291 |
| 630 | ABSOLUTE purity | all edges | miRNA_target_all | 4628 | 0.5028 | 0.4922 | 0.156 | 0.312 | 0.154 | 0.146 | 0.637 | 0.0774 | 0.0796 | 0.0746 |
| 631 | ABSOLUTE purity | all edges | miRNA_target_strong | 701 | 0.5264 | 0.5028 | 0.237 | 0.474 | 0.213 | 0.197 | 0.832 | 0.113 | 0.122 | 0.0988 |
| 632 | ABSOLUTE purity | all edges | miRNA_target_weak | 2294 | 0.5318 | 0.5006 | 0.0039 | 0.0078 | 0.00283 | 0.00206 | 0.00901 | 0.00149 | 0.00201 | 0.00107 |
| 633 | ABSOLUTE purity | all edges | miRNA_target_predicted_only | 1633 | 0.4519 | 0.4760 | 0.0544 | 0.109 | 0.0533 | 0.0518 | 0.226 | 0.976 | 0.975 | 0.976 |
| 634 | ABSOLUTE purity | Farmer-correlated targets left out | TF_target_Activation | 197 | 0.6853 | 0.6016 | 0.0187 | 0.0374 | 0.0163 | 0.0161 | 0.0702 | 0.00932 | 0.00992 | 0.00857 |
| 638 | ABSOLUTE purity | Farmer-correlated targets left out | miRNA_target_all | 3536 | 0.4745 | 0.4878 | 0.12 | 0.24 | 0.118 | 0.111 | 0.485 | 0.945 | 0.943 | 0.947 |
| 639 | ABSOLUTE purity | Farmer-correlated targets left out | miRNA_target_strong | 496 | 0.4758 | 0.4964 | 0.387 | 0.773 | 0.369 | 0.362 | 1 | 0.832 | 0.826 | 0.846 |
| 640 | ABSOLUTE purity | Farmer-correlated targets left out | miRNA_target_weak | 1780 | 0.5039 | 0.4973 | 0.592 | 1 | 0.586 | 0.616 | 1 | 0.297 | 0.3 | 0.29 |
| 641 | ABSOLUTE purity | Farmer-correlated targets left out | miRNA_target_predicted_only | 1260 | 0.4325 | 0.4710 | 0.00645 | 0.0129 | 0.00619 | 0.00634 | 0.0278 | 0.997 | 0.997 | 0.997 |

## 3. The 42 distinct estimates (`combined_tests_summary_42_distinct.csv`)

The layout is that of `concordance_summary_42_distinct.csv` (columns edges … excess_max, from the arithmetic-mean protocol), with the combined-test
columns added: <method>_p_median, _p_min, _p_max, _n_below_0.05. The reference row for all edges, miRNA_target_weak
is reproduced to 15 significant digits.

| line | edges | set | p median (a) | n p < 0.05: (a) | 2 × mean | (b) pooled | (c) Cauchy | (c) HM Vovk–Wang | (d) one-sided pooled | above null |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | all edges | miRNA_target_all | 0.555 | 0 | 0 | 0 | 1 | 0 | 2 | 23 |
| 3 | all edges | miRNA_target_predicted_only | 0.0129 | 28 | 24 | 28 | 29 | 21 | 0 | 0 |
| 4 | all edges | miRNA_target_strong | 0.55 | 0 | 0 | 1 | 1 | 0 | 3 | 35 |
| 5 | all edges | miRNA_target_weak | 0.0515 | 20 | 19 | 22 | 24 | 19 | 26 | 40 |
| 6 | all edges | TF_miRNA_Activation | 0.14 | 18 | 6 | 18 | 18 | 4 | 19 | 42 |
| 7 | all edges | TF_miRNA_Repression | 0.167 | 5 | 2 | 9 | 11 | 2 | 0 | 0 |
| 8 | all edges | TF_target_Activation | 0.0668 | 18 | 13 | 19 | 19 | 9 | 24 | 42 |
| 9 | all edges | TF_target_Repression | 0.5 | 0 | 0 | 0 | 0 | 0 | 0 | 16 |
| 10 | Farmer-correlated targets left out | miRNA_target_all | 0.558 | 0 | 0 | 0 | 0 | 0 | 0 | 7 |
| 11 | Farmer-correlated targets left out | miRNA_target_predicted_only | 0.0225 | 29 | 23 | 29 | 29 | 18 | 0 | 0 |
| 12 | Farmer-correlated targets left out | miRNA_target_strong | 0.808 | 0 | 0 | 0 | 0 | 0 | 0 | 20 |
| 13 | Farmer-correlated targets left out | miRNA_target_weak | 0.259 | 3 | 1 | 4 | 4 | 1 | 12 | 41 |
| 14 | Farmer-correlated targets left out | TF_miRNA_Activation | 0.14 | 18 | 6 | 18 | 18 | 4 | 19 | 42 |
| 15 | Farmer-correlated targets left out | TF_miRNA_Repression | 0.167 | 5 | 2 | 9 | 11 | 2 | 0 | 0 |
| 16 | Farmer-correlated targets left out | TF_target_Activation | 0.0198 | 30 | 26 | 30 | 32 | 12 | 36 | 42 |
| 17 | Farmer-correlated targets left out | TF_target_Repression | 0.485 | 0 | 0 | 0 | 0 | 0 | 0 | 19 |

## 4. Does any conclusion change? (`conclusion_changes.csv`, columns sig_<method> and changes_vs_mean)

Significance is read at p < 0.05. Relative to the arithmetic mean (a), for the key tests:

- line 9: unadjusted, all edges, miRNA_target_predicted_only: (a) significant = TRUE; differs for: pooled_one_sided;mean_one_sided;cauchy_one_sided
- line 10: unadjusted, Farmer-correlated targets left out, TF_target_Activation: (a) significant = TRUE; differs for: twice_mean;harmonic_mean_vovk_wang
- line 16: unadjusted, Farmer-correlated targets left out, miRNA_target_weak: (a) significant = FALSE; differs for: cauchy_one_sided
- line 17: unadjusted, Farmer-correlated targets left out, miRNA_target_predicted_only: (a) significant = TRUE; differs for: pooled_one_sided;mean_one_sided;cauchy_one_sided
- line 23: CAF-A, all edges, miRNA_target_strong: (a) significant = TRUE; differs for: twice_mean;harmonic_mean_vovk_wang
- line 33: CAF-A, Farmer-correlated targets left out, miRNA_target_predicted_only: (a) significant = TRUE; differs for: twice_mean;harmonic_mean_vovk_wang;pooled_one_sided;mean_one_sided;cauchy_one_sided
- line 34: CAF-B, all edges, TF_target_Activation: (a) significant = FALSE; differs for: pooled_one_sided;mean_one_sided;cauchy_one_sided
- line 41: CAF-B, all edges, miRNA_target_predicted_only: (a) significant = TRUE; differs for: twice_mean;harmonic_mean_vovk_wang;pooled_one_sided;mean_one_sided;cauchy_one_sided
- line 49: CAF-B, Farmer-correlated targets left out, miRNA_target_predicted_only: (a) significant = TRUE; differs for: harmonic_mean_vovk_wang;pooled_one_sided;mean_one_sided;cauchy_one_sided
- line 57: ESTIMATE stromal score, all edges, miRNA_target_predicted_only: (a) significant = TRUE; differs for: twice_mean;harmonic_mean_vovk_wang;pooled_one_sided;mean_one_sided;cauchy_one_sided
- line 58: ESTIMATE stromal score, Farmer-correlated targets left out, TF_target_Activation: (a) significant = TRUE; differs for: harmonic_mean_vovk_wang
- line 634: ABSOLUTE purity, Farmer-correlated targets left out, TF_target_Activation: (a) significant = TRUE; differs for: harmonic_mean_vovk_wang
- line 641: ABSOLUTE purity, Farmer-correlated targets left out, miRNA_target_predicted_only: (a) significant = TRUE; differs for: pooled_one_sided;mean_one_sided;cauchy_one_sided

Over all 720 rows (45 adjustments: unadjusted, CAF-A, CAF-B and the 42 estimates; × 8 sets × 2 edge selections): the call changes
under at least one two-sided method in 110 rows. Under the one-sided tests it also changes wherever the concordance
lies below the null, because there the one-sided p is large.

Per method, over the same rows (gained = significant only under the method; lost = significant only under (a)):

| method | rows that change | gained | lost |
|---|---|---|---|
| twice_mean | 56 | 0 | 56 |
| pooled_two_sided | 13 | 13 | 0 |
| cauchy_two_sided | 23 | 23 | 0 |
| harmonic_mean_vovk_wang | 87 | 0 | 87 |
| pooled_one_sided | 109 | 35 | 74 |
| mean_one_sided | 107 | 33 | 74 |
| cauchy_one_sided | 116 | 42 | 74 |

