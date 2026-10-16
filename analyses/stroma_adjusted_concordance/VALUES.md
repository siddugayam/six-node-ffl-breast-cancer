# VALUES: tier-stratified miRNA–target sign concordance after stromal adjustment

Lines count the header as line 1. Tables are written by `stroma_adjusted_concordance.R` (run log `run.log`).
Concordance = fraction of edges with the predicted (negative) correlation sign; null = mean null concordance over the
five matched-null draws; p = mean over the draws of the exact two-sided binomial p (p_max = the largest).

## 1. Control (`control_reproduction.csv`)

- line 2: n: max relative difference 0, reproduced TRUE
- line 3: conc_raw: max relative difference 9.77e-16, reproduced TRUE
- line 4: null_raw: max relative difference 1.35e-15, reproduced TRUE
- line 5: p_raw: max relative difference 3.51e-15, reproduced TRUE
- line 6: conc_cafA: max relative difference 1.16e-15, reproduced TRUE
- line 7: null_cafA: max relative difference 1.08e-15, reproduced TRUE
- line 8: p_cafA: max relative difference 1.87e-15, reproduced TRUE
- line 9: p_cafA_max: max relative difference 1.21e-15, reproduced TRUE
- line 10: conc_cafB: max relative difference 9.49e-16, reproduced TRUE
- line 11: null_cafB: max relative difference 1.02e-15, reproduced TRUE
- line 12: p_cafB: max relative difference 1.98e-15, reproduced TRUE
- line 13: p_cafB_max: max relative difference 3.9e-15, reproduced TRUE

The published table of 04 (raw, CAF-A, CAF-B for all 8 edge sets) is reproduced.

## 2. Per tier: unadjusted, ABSOLUTE alone, ESTIMATE alone (`concordance_summary.csv`)

Columns: n_edges; unadj_conc, unadj_null, unadj_p; ABSOLUTE_tumours, ABSOLUTE_unadj_conc, ABSOLUTE_unadj_p (unadjusted,
same 1,023 tumours); ABSOLUTE_conc, ABSOLUTE_null, ABSOLUTE_p; ESTIMATE_conc, ESTIMATE_null, ESTIMATE_p.

| line | edges | set | edges n | unadjusted conc / null (p) | same tumours as ABSOLUTE: conc (p) | ABSOLUTE-adjusted conc / null (p) | ESTIMATE-adjusted conc / null (p) |
|---|---|---|---|---|---|---|---|
| 6 | all edges | miRNA_target_all | 4628 | 0.479 / 0.474 (0.499) | 0.477 (0.549) | 0.503 / 0.492 (0.156) | 0.448 / 0.451 (0.664) |
| 7 | all edges | miRNA_target_strong | 701 | 0.495 / 0.488 (0.73) | 0.501 (0.505) | 0.526 / 0.503 (0.237) | 0.475 / 0.461 (0.473) |
| 8 | all edges | miRNA_target_weak | 2294 | 0.527 / 0.485 (0.000155) | 0.522 (0.000492) | 0.532 / 0.501 (0.0039) | 0.466 / 0.458 (0.484) |
| 9 | all edges | miRNA_target_predicted_only | 1633 | 0.404 / 0.451 (0.000183) | 0.404 (0.00022) | 0.452 / 0.476 (0.0544) | 0.411 / 0.438 (0.0311) |
| 14 | Farmer-correlated targets left out | miRNA_target_all | 3536 | 0.454 / 0.469 (0.0838) | 0.452 (0.0718) | 0.475 / 0.488 (0.12) | 0.440 / 0.446 (0.496) |
| 15 | Farmer-correlated targets left out | miRNA_target_strong | 496 | 0.446 / 0.481 (0.125) | 0.452 (0.202) | 0.476 / 0.496 (0.387) | 0.456 / 0.452 (0.858) |
| 16 | Farmer-correlated targets left out | miRNA_target_weak | 1780 | 0.502 / 0.482 (0.119) | 0.498 (0.161) | 0.504 / 0.497 (0.592) | 0.460 / 0.455 (0.654) |
| 17 | Farmer-correlated targets left out | miRNA_target_predicted_only | 1260 | 0.390 / 0.446 (8.44e-05) | 0.388 (6.82e-05) | 0.433 / 0.471 (0.00645) | 0.406 / 0.432 (0.0644) |

## 3. Per tier: the 45 stromal estimates, each alone (`concordance_summary.csv`)

Columns all45_conc_median/min/max, all45_excess_median/min/max (concordance minus null), all45_p_median/min/max,
all45_n_p_below_0.05, all45_n_above_null.

| line | edges | set | adjusted conc median (range) | excess over null median (range) | p median (range) | estimates with p < 0.05 | estimates above null |
|---|---|---|---|---|---|---|---|
| 6 | all edges | miRNA_target_all | 0.450 (0.411–0.511) | +0.0005 (-0.0121 to +0.0143) | 0.648 (0.0574–0.888) | 0 of 45 | 25 of 45 |
| 7 | all edges | miRNA_target_strong | 0.475 (0.435–0.539) | +0.0116 (-0.0086 to +0.0383) | 0.553 (0.0508–0.898) | 0 of 45 | 37 of 45 |
| 8 | all edges | miRNA_target_weak | 0.473 (0.423–0.546) | +0.0216 (-0.0056 to +0.0434) | 0.0501 (8.33e-05–0.786) | 22 of 45 | 43 of 45 |
| 9 | all edges | miRNA_target_predicted_only | 0.407 (0.376–0.452) | -0.0317 (-0.0502 to -0.0086) | 0.0105 (5.27e-05–0.497) | 30 of 45 | 0 of 45 |
| 14 | Farmer-correlated targets left out | miRNA_target_all | 0.440 (0.423–0.485) | -0.0049 (-0.0141 to +0.0034) | 0.575 (0.101–0.893) | 0 of 45 | 8 of 45 |
| 15 | Farmer-correlated targets left out | miRNA_target_strong | 0.454 (0.435–0.486) | -0.0001 (-0.0353 to +0.0141) | 0.811 (0.126–0.914) | 0 of 45 | 22 of 45 |
| 16 | Farmer-correlated targets left out | miRNA_target_weak | 0.468 (0.439–0.520) | +0.0140 (-0.0001 to +0.0290) | 0.254 (0.0184–0.825) | 4 of 45 | 44 of 45 |
| 17 | Farmer-correlated targets left out | miRNA_target_predicted_only | 0.398 (0.377–0.435) | -0.0323 (-0.0517 to -0.0116) | 0.0215 (0.000233–0.413) | 31 of 45 | 0 of 45 |

## 4. TF edge sets (same table, for completeness)

| line | edges | set | n | unadjusted conc / null (p) | ABSOLUTE conc / null (p) | ESTIMATE conc / null (p) |
|---|---|---|---|---|---|---|
| 2 | all edges | TF_target_Activation | 301 | 0.711 / 0.635 (0.00705) | 0.684 / 0.607 (0.00656) | 0.605 / 0.567 (0.206) |
| 3 | all edges | TF_target_Repression | 136 | 0.426 / 0.371 (0.204) | 0.449 / 0.406 (0.363) | 0.382 / 0.441 (0.224) |
| 4 | all edges | TF_miRNA_Activation | 708 | 0.578 / 0.512 (0.000568) | 0.547 / 0.495 (0.00629) | 0.558 / 0.537 (0.275) |
| 5 | all edges | TF_miRNA_Repression | 374 | 0.455 / 0.496 (0.134) | 0.455 / 0.507 (0.0692) | 0.422 / 0.466 (0.113) |
| 10 | Farmer-correlated targets left out | TF_target_Activation | 197 | 0.711 / 0.635 (0.0307) | 0.685 / 0.602 (0.0187) | 0.655 / 0.568 (0.0157) |
| 11 | Farmer-correlated targets left out | TF_target_Repression | 96 | 0.448 / 0.372 (0.146) | 0.469 / 0.411 (0.297) | 0.396 / 0.434 (0.497) |
| 12 | Farmer-correlated targets left out | TF_miRNA_Activation | 708 | 0.578 / 0.512 (0.000568) | 0.547 / 0.495 (0.00629) | 0.558 / 0.537 (0.275) |
| 13 | Farmer-correlated targets left out | TF_miRNA_Repression | 374 | 0.455 / 0.496 (0.134) | 0.455 / 0.507 (0.0692) | 0.422 / 0.466 (0.113) |

## 5. Every estimate

`concordance_adjusted_by_estimate.csv` holds 1440 rows: 45 estimates × 2 edge selections × 8 sets × 2 analyses
(adjusted; unadjusted on the same tumours). Columns: estimate, n_tumours, edges, set, n_edges, analysis, conc, null,
p_mean, p_max.

The Farmer-correlated set: genes with Spearman rho > 0.5 against the Farmer stroma-related ssGSEA score (698 genes,
67 of them network miRNA targets; counts printed in `run.log`).

