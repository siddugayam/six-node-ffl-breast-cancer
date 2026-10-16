# VALUES: stromal mediation of the twelve paths (seven in Table 3) on a common sample of tumours

Lines count the header as line 1. Tables are written by `mediation_common_sample.R` (run log `run.log`).

## 1. Control (`control_reproduction.csv`)

- line 2: 504 of 504 models (12 paths × 42 estimates, each on its own tumours) reproduce the published rows: n equal TRUE; largest relative difference in total 4.6e-15, ACME 4.7e-15, ADE 4.5e-15, CI bounds 4.7e-15; reproduced TRUE.

## 2. Per path: common sample against the published runs (`mediation_common_sample_summary.csv`)

Common sample: 734 tumours with both assays and all 42 distinct estimates (column n_common). The
total effect is the same for every estimate (checked by the script).

| line | path | in Table 3 (published inconsistent there) | total (common) | ACME median (range) | ADE median (range) | proportion mediated median (range) | inconsistent of 42 | of which ADE CI excludes 0 | Sobel p < 0.05 | published: n range | published total range | published inconsistent | published inconsistent, ADE CI excl. 0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | ETS1 -> COL1A1 | yes (23) | +0.275 | +0.277 (-0.108 to +0.692) | -0.002 (-0.416 to +0.383) | 1.007 (-0.391 to 2.510) | 23 | 20 | 40 | 783-1097 | 0.2797 to 0.3012 | 23 | 20 |
| 3 | NFKB1 -> COL1A1 | yes (32) | +0.112 | +0.139 (-0.021 to +0.263) | -0.028 (-0.151 to +0.132) | 1.248 (-0.184 to 2.354) | 33 | 20 | 38 | 783-1097 | 0.1003 to 0.124 | 32 | 19 |
| 4 | ETS1 -> COL3A1 | no | +0.395 | +0.262 (-0.148 to +0.659) | +0.132 (-0.264 to +0.543) | 0.664 (-0.375 to 1.669) | 11 | 8 | 39 | 783-1097 | 0.3944 to 0.4101 | 12 | 9 |
| 5 | NFKB1 -> COL3A1 | no | +0.218 | +0.143 (-0.021 to +0.272) | +0.075 (-0.054 to +0.239) | 0.655 (-0.097 to 1.246) | 10 | 7 | 38 | 783-1097 | 0.1973 to 0.2255 | 9 | 5 |
| 6 | SP1 -> COL1A1 | yes (38) | +0.011 | +0.114 (-0.018 to +0.251) | -0.103 (-0.241 to +0.028) | 10.832 (-1.701 to 23.904) | 39 | 31 | 32 | 783-1097 | 0.02088 to 0.03866 | 38 | 31 |
| 7 | RELA -> COL1A1 | yes (22) | +0.035 | -0.017 (-0.109 to +0.045) | +0.051 (-0.011 to +0.143) | -0.475 (-3.120 to 1.304) | 33 | 26 | 11 | 783-1097 | 0.03503 to 0.05168 | 22 | 19 |
| 8 | hsa-miR-29a -> COL1A1 | yes (6) | -0.238 | -0.145 (-0.238 to +0.068) | -0.093 (-0.306 to +0.000) | 0.609 (-0.285 to 1.001) | 6 | 5 | 36 | 767-1066 | -0.2481 to -0.2208 | 6 | 6 |
| 9 | hsa-miR-29a -> COL3A1 | yes (6) | -0.251 | -0.149 (-0.234 to +0.072) | -0.102 (-0.323 to -0.017) | 0.593 (-0.289 to 0.934) | 5 | 5 | 36 | 767-1066 | -0.2615 to -0.2266 | 6 | 6 |
| 10 | hsa-miR-29b -> COL1A1 | no | -0.290 | -0.168 (-0.313 to +0.002) | -0.122 (-0.292 to +0.023) | 0.581 (-0.005 to 1.079) | 2 | 1 | 37 | 767-1066 | -0.2996 to -0.2955 | 3 | 2 |
| 11 | hsa-miR-29b -> COL3A1 | yes (3) | -0.278 | -0.178 (-0.311 to +0.002) | -0.099 (-0.280 to +0.033) | 0.643 (-0.008 to 1.120) | 2 | 2 | 37 | 767-1066 | -0.2858 to -0.2818 | 3 | 3 |
| 12 | hsa-miR-29c -> COL1A1 | no | -0.262 | -0.130 (-0.283 to +0.021) | -0.132 (-0.283 to +0.021) | 0.496 (-0.081 to 1.081) | 5 | 4 | 35 | 767-1066 | -0.269 to -0.2503 | 7 | 6 |
| 13 | hsa-miR-29c -> COL3A1 | no | -0.270 | -0.136 (-0.279 to +0.021) | -0.133 (-0.291 to +0.010) | 0.505 (-0.080 to 1.036) | 5 | 4 | 35 | 767-1066 | -0.2786 to -0.2607 | 7 | 6 |

Columns: total_common; ACME_median, ACME_min, ACME_max; ADE_*; propmed_*; n_inconsistent;
n_inconsistent_ADE_CI_excl0; n_sobel_p05; and published_* (the same 42 estimates in
results/systems_and_dynamics/deconv_mediation_extended.csv).

The proportion mediated is ACME / total. It is unstable where the total is near 0.

The script checks that the published inconsistent counts of Table 3's seven rows are 23, 32, 38, 22, 6, 6 and 3
(columns published_n_inconsistent and withexemplar_inconsistent).

## 3. Every path × estimate (`mediation_common_sample_per_estimate.csv`)

504 rows, with columns x, y, mediator, n, total, ACME (acme_lo, acme_hi, acme_p),
ADE (ade_lo, ade_hi, ade_p), prop_mediated (prop_lo, prop_hi), sobel_p, inconsistent and ADE_CI_excludes_zero.

The three exact repeats left out are MCPcounter_Fibroblasts_recomputed (= MCPcounter_Fibroblasts), NNLS_Wu641_full_CAFs
(= NNLS_Wu_CAFs) and QPROG_Wu641_full_CAFs (= QPROG_Wu_CAFs); equality is checked in `run.log`.

