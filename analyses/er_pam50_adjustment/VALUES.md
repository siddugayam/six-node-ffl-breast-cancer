# VALUES: oestrogen-receptor and PAM50 adjustment

Lines count the header as line 1. Part 1 tables are written by `circuit_survival_er_pam50.R` (log `run_circuits.log`),
part 2 tables by `mir29a_er_pam50.R` (log `run_mir29a.log`).

## Part 1. Six-node circuits, TCGA-BRCA, overall survival

### Control (`control_reproduction.csv`)

- line 2: 2454 of 2454 published rows compared; n and events equal TRUE; largest relative differences HR 5.0e-15, p 4.7e-15, LRT p 4.7e-15, q 4.2e-15; reproduced TRUE.

### Summary per fit and family (`circuit_summary.csv`)

Fits: published = age + stage on the published sample; same_tumours = age + stage on the tumours with ER and PAM50;
er_pam50 = age + stage + ER + PAM50 on those tumours. LRT = core vs core + extension; BH within fit and family.

| line | fit | family | circuits | n | events | LRT p < 0.05 | q < 0.05 | HR p < 0.05 | the 16: LRT p < 0.05 | the 16: q < 0.05 | the 16: HR p < 0.05 | median ΔC |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | published | six-node composite pattern | 2383 | 1028 | 137 | 288 | 16 | 467 | 16 | 16 | 14 | +0.00043 |
| 3 | published | six-node model architecture | 71 | 1028 | 137 | 9 | 0 | 10 | 0 | 0 | 0 | +0.00036 |
| 4 | same_tumours | six-node composite pattern | 2383 | 747 | 103 | 376 | 0 | 426 | 16 | 0 | 8 | +0.00106 |
| 5 | same_tumours | six-node model architecture | 71 | 747 | 103 | 0 | 0 | 1 | 0 | 0 | 0 | +0.00048 |
| 6 | er_pam50 | six-node composite pattern | 2383 | 747 | 103 | 386 | 2 | 372 | 16 | 0 | 7 | +0.00048 |
| 7 | er_pam50 | six-node model architecture | 71 | 747 | 103 | 0 | 0 | 0 | 0 | 0 | 0 | -0.00032 |

In the six-node model architecture family, the columns for the 16 count none, because the 16 are all six-node
composite pattern circuits. Median ΔC is the median over circuits of C_six − C_core (column median_deltaC).

### Circuits with q < 0.05 in the same-tumour or ER + PAM50 fits (`circuits_er_pam50.csv`)

| line | fit | family | members | core | extension | HR (95 % CI) | p | LRT p | q |
|---|---|---|---|---|---|---|---|---|---|
| 5837 | er_pam50 | six-node composite pattern | COL1A1;HIF1A;hsa-miR-29b;hsa-miR-29c;SPARC;VHL | HIF1A;hsa-miR-29c;COL1A1 | VHL;hsa-miR-29b;SPARC | 1.351 (1.071–1.705) | 0.0112 | 2.6e-05 | 0.037 |
| 5838 | er_pam50 | six-node composite pattern | COL3A1;HIF1A;hsa-miR-29b;hsa-miR-29c;SPARC;VHL | HIF1A;hsa-miR-29c;COL3A1 | VHL;hsa-miR-29b;SPARC | 1.365 (1.084–1.718) | 0.00817 | 3.11e-05 | 0.037 |

Of these 2, 0 are among the 16 patterns.

### The 16 patterns before and after (`sixteen_patterns_before_after.csv`)

HR of the six-node score per SD (95 % CI), its p, and the LRT p (q) of each fit.

| line | members | published: HR (CI), p; LRT p (q) | same tumours: HR (CI), p; LRT p (q) | ER + PAM50: HR (CI), p; LRT p (q) |
|---|---|---|---|---|
| 2 | ADA;HIF1A;hsa-miR-29b;hsa-miR-29c;NT5E;STAT3 | 1.100 (0.922–1.314), 0.29; 0.000158 (0.0387) | 1.094 (0.888–1.348), 0.397; 0.0001 (0.0734) | 1.076 (0.869–1.334), 0.501; 0.000235 (0.103) |
| 3 | ARNT;F3;HIF1A;hsa-miR-106b;hsa-miR-93;SERPINE1 | 1.274 (1.064–1.525), 0.00825; 8.19e-05 (0.0387) | 1.192 (0.967–1.469), 0.1; 0.000956 (0.104) | 1.159 (0.925–1.452), 0.2; 0.00241 (0.153) |
| 4 | BCL2L1;GATA3;HRK;hsa-miR-29b;hsa-miR-29c;RELA | 0.788 (0.655–0.948), 0.0114; 0.000161 (0.0387) | 0.875 (0.710–1.077), 0.208; 0.0237 (0.259) | 0.890 (0.710–1.116), 0.314; 0.0273 (0.233) |
| 5 | E2F1;FN1;hsa-miR-200b;hsa-miR-429;THBS1;TP73 | 1.131 (0.951–1.345), 0.165; 0.000103 (0.0387) | 1.145 (0.933–1.406), 0.195; 0.000565 (0.0938) | 1.186 (0.948–1.484), 0.136; 0.000783 (0.118) |
| 6 | F3;HIF1A;hsa-miR-17;hsa-miR-19a;SERPINE1;STAT3 | 1.327 (1.106–1.591), 0.00227; 0.000304 (0.0474) | 1.270 (1.028–1.568), 0.0264; 0.000944 (0.104) | 1.248 (0.995–1.566), 0.0549; 0.00365 (0.157) |
| 7 | F3;HIF1A;hsa-miR-17;hsa-miR-19b;SERPINE1;STAT3 | 1.328 (1.109–1.590), 0.00202; 0.000186 (0.0403) | 1.280 (1.038–1.579), 0.0212; 0.000646 (0.0938) | 1.272 (1.017–1.590), 0.0351; 0.00203 (0.145) |
| 8 | F3;HIF1A;hsa-miR-18a;hsa-miR-19a;SERPINE1;STAT3 | 1.253 (1.045–1.501), 0.0147; 6.76e-05 (0.0387) | 1.213 (0.983–1.498), 0.0721; 0.000318 (0.0938) | 1.176 (0.938–1.475), 0.161; 0.00125 (0.136) |
| 9 | F3;HIF1A;hsa-miR-18a;hsa-miR-19b;SERPINE1;STAT3 | 1.270 (1.060–1.521), 0.00965; 2.88e-05 (0.0387) | 1.235 (1.000–1.525), 0.0502; 0.000197 (0.0925) | 1.205 (0.963–1.508), 0.103; 0.000588 (0.118) |
| 10 | F3;HIF1A;hsa-miR-18a;hsa-miR-20a;SERPINE1;STAT3 | 1.262 (1.054–1.512), 0.0115; 0.000163 (0.0387) | 1.222 (0.988–1.512), 0.0647; 0.00102 (0.106) | 1.202 (0.958–1.507), 0.112; 0.0021 (0.145) |
| 11 | F3;HIF1A;hsa-miR-19a;hsa-miR-19b;SERPINE1;STAT3 | 1.266 (1.057–1.517), 0.0104; 4.53e-05 (0.0387) | 1.282 (1.038–1.583), 0.0212; 0.000714 (0.0938) | 1.266 (1.011–1.585), 0.0395; 0.00205 (0.145) |
| 12 | F3;HIF1A;hsa-miR-19a;hsa-miR-20a;SERPINE1;STAT3 | 1.259 (1.052–1.506), 0.0119; 0.000231 (0.0424) | 1.264 (1.024–1.560), 0.0295; 0.00294 (0.171) | 1.261 (1.007–1.578), 0.0433; 0.00589 (0.157) |
| 13 | F3;HIF1A;hsa-miR-19b;hsa-miR-20a;SERPINE1;STAT3 | 1.280 (1.070–1.531), 0.00684; 4.91e-05 (0.0387) | 1.285 (1.041–1.585), 0.0197; 0.000748 (0.0938) | 1.291 (1.034–1.613), 0.0244; 0.00294 (0.157) |
| 14 | HDAC2;hsa-miR-15b;hsa-miR-16;MYC;SMAD7;TGFBR2 | 1.284 (1.080–1.527), 0.00463; 0.000159 (0.0387) | 1.270 (1.038–1.553), 0.0201; 0.00641 (0.222) | 1.313 (1.040–1.657), 0.0218; 0.0212 (0.214) |
| 15 | HDAC2;hsa-miR-25;hsa-miR-93;MYC;SMAD7;TGFBR2 | 1.247 (1.050–1.482), 0.0119; 0.000291 (0.0474) | 1.182 (0.972–1.437), 0.093; 0.000734 (0.0938) | 1.192 (0.949–1.497), 0.131; 0.00307 (0.157) |
| 16 | HIF1A;hsa-miR-18a;hsa-miR-19b;MMP2;STAT3;THBS1 | 1.380 (1.138–1.673), 0.00105; 0.000217 (0.0424) | 1.414 (1.122–1.784), 0.00339; 0.00177 (0.132) | 1.408 (1.102–1.799), 0.0062; 0.00569 (0.157) |
| 17 | HIF1A;hsa-miR-19a;hsa-miR-19b;MMP2;STAT3;THBS1 | 1.366 (1.128–1.655), 0.00141; 0.000318 (0.0474) | 1.451 (1.154–1.825), 0.00144; 0.00481 (0.198) | 1.452 (1.141–1.847), 0.0024; 0.0143 (0.183) |

Every circuit and fit: `circuits_er_pam50.csv` (columns fit, covariates, family, members, core, extension, n, events, HR_six, HR_lo, HR_hi,
p_six, C_core, C_six, LRT_p_extension_adds, q_BH). METABRIC: NOT POSSIBLE (no miRNA assay; the circuits cannot be scored).

## Part 2. hsa-miR-29a and the other miRNAs, adjusted for ER and/or PAM50

### Pooled (`mir29a_er_pam50_pooled.csv`)

| line | analysis | miRNA | k | n | events | HR FE (95 % CI) | p FE | HR RE (95 % CI) | p RE | q RE | I² % |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | unadjusted, same cases | miR-21 | 5 | 1288 | 370 | 0.973 (0.874–1.083) | 0.616 | 0.973 (0.874–1.083) | 0.616 | 0.646 | 0.0 |
| 3 | unadjusted, same cases | miR-195 | 5 | 1288 | 370 | 0.826 (0.746–0.914) | 0.000217 | 0.826 (0.746–0.914) | 0.000217 | 0.00119 | 0.0 |
| 4 | unadjusted, same cases | miR-204 | 5 | 1288 | 370 | 0.840 (0.752–0.939) | 0.00209 | 0.850 (0.720–1.003) | 0.0537 | 0.197 | 53.8 |
| 5 | unadjusted, same cases | miR-383 | 3 | 1101 | 250 | 0.870 (0.751–1.009) | 0.0647 | 0.847 (0.673–1.065) | 0.155 | 0.283 | 55.9 |
| 6 | unadjusted, same cases | miR-124 | 2 | 891 | 171 | 1.099 (0.976–1.237) | 0.118 | 1.099 (0.976–1.237) | 0.118 | 0.259 | 0.0 |
| 7 | unadjusted, same cases | miR-155 | 4 | 1193 | 292 | 0.902 (0.804–1.013) | 0.0809 | 0.902 (0.804–1.013) | 0.0809 | 0.222 | 0.0 |
| 8 | unadjusted, same cases | miR-429 | 5 | 1288 | 370 | 1.100 (0.986–1.228) | 0.088 | 1.102 (0.878–1.382) | 0.403 | 0.506 | 76.3 |
| 9 | unadjusted, same cases | miR-141 | 5 | 1285 | 367 | 0.981 (0.877–1.098) | 0.744 | 0.964 (0.827–1.125) | 0.646 | 0.646 | 44.8 |
| 10 | unadjusted, same cases | miR-34a | 5 | 1288 | 370 | 1.087 (0.978–1.208) | 0.123 | 1.068 (0.912–1.250) | 0.414 | 0.506 | 53.2 |
| 11 | unadjusted, same cases | miR-101 | 5 | 1285 | 367 | 0.873 (0.782–0.974) | 0.0153 | 0.892 (0.719–1.106) | 0.297 | 0.467 | 73.0 |
| 12 | unadjusted, same cases | miR-29a | 5 | 1288 | 370 | 0.807 (0.733–0.889) | 1.44e-05 | 0.807 (0.733–0.889) | 1.44e-05 | 0.000158 | 0.0 |
| 13 | adjusted for ER and/or PAM50 | miR-21 | 5 | 1288 | 370 | 0.963 (0.865–1.074) | 0.501 | 0.963 (0.865–1.074) | 0.501 | 0.551 | 0.0 |
| 14 | adjusted for ER and/or PAM50 | miR-195 | 5 | 1288 | 370 | 0.835 (0.751–0.930) | 0.000999 | 0.835 (0.751–0.930) | 0.000999 | 0.0055 | 0.0 |
| 15 | adjusted for ER and/or PAM50 | miR-204 | 5 | 1288 | 370 | 0.829 (0.740–0.929) | 0.00128 | 0.835 (0.725–0.962) | 0.0125 | 0.046 | 33.5 |
| 16 | adjusted for ER and/or PAM50 | miR-383 | 3 | 1101 | 250 | 0.879 (0.758–1.019) | 0.0881 | 0.854 (0.679–1.073) | 0.175 | 0.321 | 54.9 |
| 17 | adjusted for ER and/or PAM50 | miR-124 | 2 | 891 | 171 | 1.090 (0.965–1.231) | 0.165 | 1.090 (0.965–1.231) | 0.165 | 0.321 | 0.0 |
| 18 | adjusted for ER and/or PAM50 | miR-155 | 4 | 1193 | 292 | 0.868 (0.767–0.983) | 0.0252 | 0.868 (0.767–0.983) | 0.0252 | 0.0694 | 0.0 |
| 19 | adjusted for ER and/or PAM50 | miR-429 | 5 | 1288 | 370 | 1.101 (0.983–1.233) | 0.0953 | 1.099 (0.870–1.390) | 0.428 | 0.523 | 76.4 |
| 20 | adjusted for ER and/or PAM50 | miR-141 | 5 | 1285 | 367 | 1.001 (0.893–1.124) | 0.981 | 0.995 (0.874–1.133) | 0.943 | 0.943 | 19.5 |
| 21 | adjusted for ER and/or PAM50 | miR-34a | 5 | 1288 | 370 | 1.095 (0.982–1.221) | 0.103 | 1.080 (0.941–1.239) | 0.273 | 0.429 | 34.5 |
| 22 | adjusted for ER and/or PAM50 | miR-101 | 5 | 1285 | 367 | 0.898 (0.800–1.007) | 0.0653 | 0.913 (0.735–1.134) | 0.409 | 0.523 | 70.7 |
| 23 | adjusted for ER and/or PAM50 | miR-29a | 5 | 1288 | 370 | 0.811 (0.736–0.895) | 2.9e-05 | 0.811 (0.736–0.895) | 2.9e-05 | 0.000319 | 0.0 |

### Per cohort, hsa-miR-29a (`mir29a_er_pam50_percohort.csv`)

| line | cohort | endpoint | model | covariates | n | events | HR (95 % CI) | p | Schoenfeld miRNA p |
|---|---|---|---|---|---|---|---|---|---|
| 22 | TCGA_BRCA | OS | unadjusted_same_cases | none (same cases as the adjusted model) | 768 | 112 | 0.816 (0.684–0.974) | 0.0242 | 0.248 |
| 23 | TCGA_BRCA | OS | adjusted | ER (IHC); PAM50 | 768 | 112 | 0.824 (0.686–0.988) | 0.0371 | 0.202 |
| 40 | GSE19783 | OS | unadjusted_same_cases | none (same cases as the adjusted model) | 92 | 42 | 0.904 (0.665–1.229) | 0.52 | 0.165 |
| 41 | GSE19783 | OS | adjusted | ER; PAM50 subtype | 92 | 42 | 0.932 (0.678–1.282) | 0.667 | 0.195 |
| 60 | GSE22216 | DRFS | unadjusted_same_cases | none (same cases as the adjusted model) | 210 | 79 | 0.780 (0.656–0.927) | 0.00486 | 0.737 |
| 61 | GSE22216 | DRFS | adjusted | ER | 210 | 79 | 0.765 (0.643–0.912) | 0.00272 | 0.744 |
| 82 | GSE59829 | DMFS | unadjusted_same_cases | none (same cases as the adjusted model) | 123 | 59 | 0.875 (0.683–1.122) | 0.292 | 0.832 |
| 83 | GSE59829 | DMFS | adjusted | ER | 123 | 59 | 0.878 (0.684–1.127) | 0.307 | 0.825 |
| 98 | GSE78870 | TTP | unadjusted_same_cases | none (same cases as the adjusted model) | 95 | 78 | 0.721 (0.556–0.934) | 0.0134 | 0.28 |
| 99 | GSE78870 | TTP | adjusted | ER | 95 | 78 | 0.750 (0.581–0.969) | 0.0276 | 0.413 |

GSE37405 is ER-positive by design and has no subtype annotation; it is not in this pool.

