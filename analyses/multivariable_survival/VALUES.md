# VALUES: multivariable survival for the miRNAs of the survival meta-analysis

Lines count the header as line 1. Every table is written by `mirna_multivariable_cox.R` (run log `run.log`).
HR is per SD of the miRNA (z-scored within cohort, within platform for GSE37405).

## 1. Control: Table S8 reproduced (`table_s8_reproduction.csv`)

The unadjusted primary-endpoint pool equals Table S8 for all 11 miRNAs (HR, CI, p, I² and q within relative
tolerance 1e-9; column `identical_within_1e9`).

| line | miRNA | k | HR (RE) | 95 % CI | p | identical |
|---|---|---|---|---|---|---|
| 2 | miR-101 | 6 | 0.903 | 0.777–1.050 | 0.186 | TRUE |
| 3 | miR-124 | 3 | 0.994 | 0.844–1.170 | 0.938 | TRUE |
| 4 | miR-141 | 6 | 0.939 | 0.853–1.034 | 0.203 | TRUE |
| 5 | miR-155 | 5 | 0.933 | 0.844–1.033 | 0.181 | TRUE |
| 6 | miR-195 | 6 | 0.833 | 0.753–0.923 | 0.000444 | TRUE |
| 7 | miR-204 | 6 | 0.795 | 0.663–0.953 | 0.0131 | TRUE |
| 8 | miR-21 | 6 | 1.002 | 0.912–1.101 | 0.964 | TRUE |
| 9 | miR-29a | 6 | 0.826 | 0.758–0.900 | 1.36e-05 | TRUE |
| 10 | miR-34a | 6 | 1.026 | 0.898–1.173 | 0.704 | TRUE |
| 11 | miR-383 | 4 | 0.873 | 0.753–1.013 | 0.0734 | TRUE |
| 12 | miR-429 | 6 | 1.066 | 0.901–1.261 | 0.455 | TRUE |

## 2. Covariates available per cohort (`cohort_covariates.csv`; columns n_samples, n_available, levels)

| line | cohort | covariate | samples | available | levels |
|---|---|---|---|---|---|
| 2 | TCGA_BRCA | age | 1068 | 1066 |  |
| 3 | TCGA_BRCA | stage4 | 1068 | 1044 | I 178; II 604; III 243; IV 19 |
| 4 | TCGA_BRCA | ER_ihc | 1068 | 1017 | Positive 786; Negative 231 |
| 5 | TCGA_BRCA | PAM50 | 1068 | 817 | LumA 414; LumB 182; Her2 66; Basal 132; Normal 23 |
| 6 | TCGA_BRCA | stromal | 1068 | 1066 |  |
| 7 | GSE19783 | ER2 | 101 | 99 | Positive 61; Negative 38 |
| 8 | GSE19783 | subtype5 | 101 | 95 | Lum A 41; Lum B 12; ERBB2 17; Basal 15; Normal-like 10 |
| 9 | GSE19783 | stromal | 101 | 101 |  |
| 10 | GSE22216 | age | 210 | 210 |  |
| 11 | GSE22216 | grade3 | 210 | 186 | 1 42; 2 81; 3 63 |
| 12 | GSE22216 | size | 210 | 210 |  |
| 13 | GSE22216 | node | 210 | 210 |  |
| 14 | GSE22216 | ER | 210 | 210 |  |
| 15 | GSE22216 | stromal | 210 | 205 |  |
| 16 | GSE37405 | age | 153 | 152 |  |
| 17 | GSE37405 | size | 153 | 151 |  |
| 18 | GSE37405 | node | 153 | 151 |  |
| 19 | GSE59829 | age | 123 | 123 |  |
| 20 | GSE59829 | size | 123 | 122 |  |
| 21 | GSE59829 | ER2 | 123 | 123 | 0 23; 1 100 |
| 22 | GSE78870 | age | 106 | 105 |  |
| 23 | GSE78870 | grade3 | 106 | 92 | 1 11; 2 45; 3 36 |
| 24 | GSE78870 | ER2 | 106 | 103 | positive 95; negative 8 |
| 25 | GSE19783 | stromal score genes present |  | 118 | of 123 |
| 26 | GSE22216 | stromal score genes present |  | 102 | of 123 |

METABRIC: NOT POSSIBLE (no miRNA profiles in its public release). No cohort other than TCGA-BRCA has stage.

## 3. Pooled adjusted HR, primary endpoint (`pooled_adjusted.csv`)

Columns HR_FE, lo_FE, hi_FE, p_FE, q_BH_FE (fixed effect) and HR_RE, lo_RE, hi_RE, p_RE, q_BH_RE (DerSimonian–Laird),
I2, k. BH is over the 11 miRNAs, separately for FE and RE.

| line | analysis | miRNA | k | HR FE (95 % CI) | p FE | q FE | HR RE (95 % CI) | p RE | q RE | I² % |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | primary endpoint, unadjusted | miR-21 | 6 | 1.002 (0.912–1.101) | 0.964 | 0.964 | 1.002 (0.912–1.101) | 0.964 | 0.964 | 0.0 |
| 3 | primary endpoint, unadjusted | miR-195 | 6 | 0.833 (0.762–0.910) | 5.05e-05 | 0.000185 | 0.833 (0.753–0.923) | 0.000444 | 0.00244 | 20.7 |
| 4 | primary endpoint, unadjusted | miR-204 | 6 | 0.798 (0.723–0.880) | 7.09e-06 | 7.49e-05 | 0.795 (0.663–0.953) | 0.0131 | 0.048 | 68.7 |
| 5 | primary endpoint, unadjusted | miR-383 | 4 | 0.881 (0.783–0.993) | 0.0374 | 0.0824 | 0.873 (0.753–1.013) | 0.0734 | 0.202 | 32.1 |
| 6 | primary endpoint, unadjusted | miR-124 | 3 | 1.027 (0.926–1.139) | 0.616 | 0.678 | 0.994 (0.844–1.170) | 0.938 | 0.964 | 44.8 |
| 7 | primary endpoint, unadjusted | miR-155 | 5 | 0.933 (0.844–1.033) | 0.181 | 0.279 | 0.933 (0.844–1.033) | 0.181 | 0.318 | 0.0 |
| 8 | primary endpoint, unadjusted | miR-429 | 6 | 1.066 (0.968–1.174) | 0.192 | 0.279 | 1.066 (0.901–1.261) | 0.455 | 0.625 | 66.3 |
| 9 | primary endpoint, unadjusted | miR-141 | 6 | 0.939 (0.853–1.034) | 0.203 | 0.279 | 0.939 (0.853–1.034) | 0.203 | 0.318 | 0.0 |
| 10 | primary endpoint, unadjusted | miR-34a | 6 | 1.043 (0.951–1.144) | 0.373 | 0.456 | 1.026 (0.898–1.173) | 0.704 | 0.861 | 50.1 |
| 11 | primary endpoint, unadjusted | miR-101 | 6 | 0.880 (0.801–0.965) | 0.00692 | 0.019 | 0.903 (0.777–1.050) | 0.186 | 0.318 | 59.6 |
| 12 | primary endpoint, unadjusted | miR-29a | 6 | 0.826 (0.758–0.900) | 1.36e-05 | 7.49e-05 | 0.826 (0.758–0.900) | 1.36e-05 | 0.00015 | 0.0 |
| 13 | primary endpoint, adjusted (maximal available) | miR-21 | 6 | 0.970 (0.874–1.076) | 0.564 | 0.69 | 0.970 (0.874–1.076) | 0.564 | 0.776 | 0.0 |
| 14 | primary endpoint, adjusted (maximal available) | miR-195 | 6 | 0.842 (0.758–0.936) | 0.00136 | 0.005 | 0.842 (0.758–0.936) | 0.00136 | 0.0075 | 0.0 |
| 15 | primary endpoint, adjusted (maximal available) | miR-204 | 6 | 0.797 (0.711–0.892) | 8.14e-05 | 0.000448 | 0.789 (0.641–0.971) | 0.025 | 0.0688 | 68.9 |
| 16 | primary endpoint, adjusted (maximal available) | miR-383 | 4 | 0.868 (0.750–1.004) | 0.0562 | 0.103 | 0.864 (0.731–1.021) | 0.0854 | 0.188 | 23.2 |
| 17 | primary endpoint, adjusted (maximal available) | miR-124 | 3 | 1.075 (0.959–1.204) | 0.214 | 0.336 | 1.006 (0.803–1.259) | 0.961 | 0.961 | 66.5 |
| 18 | primary endpoint, adjusted (maximal available) | miR-155 | 5 | 0.856 (0.757–0.967) | 0.0126 | 0.0276 | 0.856 (0.748–0.980) | 0.024 | 0.0688 | 16.1 |
| 19 | primary endpoint, adjusted (maximal available) | miR-429 | 6 | 1.052 (0.944–1.173) | 0.361 | 0.496 | 1.067 (0.873–1.304) | 0.527 | 0.776 | 70.4 |
| 20 | primary endpoint, adjusted (maximal available) | miR-141 | 6 | 1.015 (0.903–1.140) | 0.808 | 0.808 | 0.992 (0.826–1.191) | 0.931 | 0.961 | 57.2 |
| 21 | primary endpoint, adjusted (maximal available) | miR-34a | 6 | 0.979 (0.881–1.089) | 0.7 | 0.77 | 0.973 (0.846–1.119) | 0.703 | 0.859 | 41.2 |
| 22 | primary endpoint, adjusted (maximal available) | miR-101 | 6 | 0.855 (0.766–0.953) | 0.00486 | 0.0134 | 0.866 (0.728–1.031) | 0.106 | 0.195 | 59.0 |
| 23 | primary endpoint, adjusted (maximal available) | miR-29a | 6 | 0.811 (0.735–0.895) | 2.94e-05 | 0.000324 | 0.811 (0.735–0.895) | 2.94e-05 | 0.000324 | 0.0 |

Adjusted, random effects, p < 0.05 and q < 0.05: miR-195, miR-29a (lines 14, 23).

## 4. Harmonised endpoint: overall survival only (`pooled_os_harmonised.csv`)

| line | analysis | miRNA | cohorts | HR FE (95 % CI) | p FE | HR RE (95 % CI) | p RE | q RE | I² % |
|---|---|---|---|---|---|---|---|---|---|
| 2 | OS only, unadjusted | miR-21 | TCGA_BRCA;GSE19783;GSE37405 | 1.058 (0.928–1.206) | 0.4 | 1.058 (0.928–1.206) | 0.4 | 0.489 | 0.0 |
| 3 | OS only, unadjusted | miR-195 | TCGA_BRCA;GSE19783;GSE37405 | 0.892 (0.791–1.007) | 0.0649 | 0.908 (0.775–1.063) | 0.228 | 0.489 | 32.5 |
| 4 | OS only, unadjusted | miR-204 | TCGA_BRCA;GSE19783;GSE37405 | 0.846 (0.737–0.973) | 0.0186 | 0.833 (0.556–1.246) | 0.373 | 0.489 | 85.2 |
| 5 | OS only, unadjusted | miR-383 | TCGA_BRCA;GSE37405 | 0.844 (0.721–0.987) | 0.0337 | 0.844 (0.721–0.987) | 0.0337 | 0.165 | 0.0 |
| 6 | OS only, unadjusted | miR-124 | TCGA_BRCA;GSE37405 | 1.017 (0.903–1.145) | 0.782 | 0.847 (0.505–1.422) | 0.531 | 0.584 | 85.4 |
| 7 | OS only, unadjusted | miR-155 | TCGA_BRCA;GSE19783;GSE37405 | 0.905 (0.799–1.026) | 0.12 | 0.905 (0.799–1.026) | 0.12 | 0.331 | 0.0 |
| 8 | OS only, unadjusted | miR-429 | TCGA_BRCA;GSE19783;GSE37405 | 1.039 (0.915–1.179) | 0.556 | 1.030 (0.862–1.231) | 0.743 | 0.743 | 47.2 |
| 9 | OS only, unadjusted | miR-141 | TCGA_BRCA;GSE19783;GSE37405 | 0.931 (0.818–1.060) | 0.279 | 0.925 (0.795–1.077) | 0.317 | 0.489 | 23.3 |
| 10 | OS only, unadjusted | miR-34a | TCGA_BRCA;GSE19783;GSE37405 | 1.109 (0.977–1.258) | 0.109 | 1.088 (0.923–1.283) | 0.317 | 0.489 | 35.6 |
| 11 | OS only, unadjusted | miR-101 | TCGA_BRCA;GSE19783;GSE37405 | 0.816 (0.721–0.923) | 0.00126 | 0.816 (0.721–0.923) | 0.00126 | 0.0139 | 0.0 |
| 12 | OS only, unadjusted | miR-29a | TCGA_BRCA;GSE19783;GSE37405 | 0.883 (0.782–0.997) | 0.0451 | 0.883 (0.782–0.997) | 0.0451 | 0.165 | 0.0 |
| 13 | OS only, adjusted (maximal available) | miR-21 | TCGA_BRCA;GSE19783;GSE37405 | 0.987 (0.852–1.143) | 0.861 | 0.983 (0.833–1.160) | 0.838 | 0.838 | 19.1 |
| 14 | OS only, adjusted (maximal available) | miR-195 | TCGA_BRCA;GSE19783;GSE37405 | 0.846 (0.731–0.979) | 0.0247 | 0.846 (0.731–0.979) | 0.0247 | 0.0751 | 0.0 |
| 15 | OS only, adjusted (maximal available) | miR-204 | TCGA_BRCA;GSE19783;GSE37405 | 0.805 (0.683–0.950) | 0.01 | 0.768 (0.466–1.264) | 0.299 | 0.547 | 87.1 |
| 16 | OS only, adjusted (maximal available) | miR-383 | TCGA_BRCA;GSE37405 | 0.797 (0.652–0.975) | 0.0273 | 0.797 (0.652–0.975) | 0.0273 | 0.0751 | 0.0 |
| 17 | OS only, adjusted (maximal available) | miR-124 | TCGA_BRCA;GSE37405 | 1.078 (0.942–1.233) | 0.276 | 0.879 (0.484–1.593) | 0.67 | 0.838 | 88.0 |
| 18 | OS only, adjusted (maximal available) | miR-155 | TCGA_BRCA;GSE19783;GSE37405 | 0.829 (0.707–0.971) | 0.0204 | 0.829 (0.707–0.971) | 0.0204 | 0.0751 | 0.0 |
| 19 | OS only, adjusted (maximal available) | miR-429 | TCGA_BRCA;GSE19783;GSE37405 | 1.058 (0.912–1.228) | 0.456 | 1.057 (0.766–1.457) | 0.738 | 0.838 | 78.6 |
| 20 | OS only, adjusted (maximal available) | miR-141 | TCGA_BRCA;GSE19783;GSE37405 | 1.007 (0.849–1.194) | 0.94 | 0.920 (0.653–1.296) | 0.633 | 0.838 | 70.7 |
| 21 | OS only, adjusted (maximal available) | miR-34a | TCGA_BRCA;GSE19783;GSE37405 | 1.055 (0.907–1.227) | 0.486 | 1.033 (0.824–1.294) | 0.78 | 0.838 | 51.9 |
| 22 | OS only, adjusted (maximal available) | miR-101 | TCGA_BRCA;GSE19783;GSE37405 | 0.763 (0.655–0.888) | 0.000497 | 0.763 (0.655–0.888) | 0.000497 | 0.00546 | 0.0 |
| 23 | OS only, adjusted (maximal available) | miR-29a | TCGA_BRCA;GSE19783;GSE37405 | 0.861 (0.743–0.997) | 0.0449 | 0.861 (0.743–0.997) | 0.0449 | 0.0987 | 0.0 |

## 5. Adjusted pool without the stromal score (`pooled_adjusted_without_stromal.csv`)

| line | analysis | miRNA | cohorts | HR FE (95 % CI) | p FE | HR RE (95 % CI) | p RE | q RE | I² % |
|---|---|---|---|---|---|---|---|---|---|
| 2 | primary endpoint, adjusted without the stromal score | miR-21 | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 0.972 (0.877–1.078) | 0.596 | 0.972 (0.877–1.078) | 0.596 | 0.841 | 0.0 |
| 3 | primary endpoint, adjusted without the stromal score | miR-195 | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 0.848 (0.763–0.941) | 0.00191 | 0.848 (0.763–0.941) | 0.00191 | 0.0105 | 0.0 |
| 4 | primary endpoint, adjusted without the stromal score | miR-204 | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 0.813 (0.727–0.909) | 0.000267 | 0.808 (0.658–0.993) | 0.0423 | 0.116 | 69.4 |
| 5 | primary endpoint, adjusted without the stromal score | miR-383 | TCGA_BRCA;GSE22216;GSE37405;GSE59829 | 0.901 (0.781–1.038) | 0.149 | 0.884 (0.718–1.089) | 0.246 | 0.452 | 52.1 |
| 6 | primary endpoint, adjusted without the stromal score | miR-124 | TCGA_BRCA;GSE37405;GSE59829 | 1.070 (0.956–1.197) | 0.243 | 1.004 (0.806–1.250) | 0.974 | 0.974 | 65.1 |
| 7 | primary endpoint, adjusted without the stromal score | miR-155 | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829 | 0.868 (0.771–0.978) | 0.0201 | 0.868 (0.763–0.989) | 0.0333 | 0.116 | 14.0 |
| 8 | primary endpoint, adjusted without the stromal score | miR-429 | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 1.040 (0.935–1.157) | 0.468 | 1.049 (0.871–1.264) | 0.612 | 0.841 | 67.0 |
| 9 | primary endpoint, adjusted without the stromal score | miR-141 | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 0.996 (0.893–1.111) | 0.939 | 0.987 (0.848–1.148) | 0.864 | 0.951 | 46.4 |
| 10 | primary endpoint, adjusted without the stromal score | miR-34a | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 0.982 (0.884–1.092) | 0.74 | 0.975 (0.848–1.122) | 0.727 | 0.888 | 42.1 |
| 11 | primary endpoint, adjusted without the stromal score | miR-101 | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 0.855 (0.766–0.953) | 0.0047 | 0.866 (0.729–1.029) | 0.103 | 0.226 | 58.9 |
| 12 | primary endpoint, adjusted without the stromal score | miR-29a | TCGA_BRCA;GSE19783;GSE22216;GSE37405;GSE59829;GSE78870 | 0.820 (0.744–0.905) | 7.5e-05 | 0.820 (0.744–0.905) | 7.5e-05 | 0.000825 | 0.0 |

The OS-only analysis keeps TCGA-BRCA (OS), GSE19783 (OS; the deposited time variable is the disease-free time) and
GSE37405 (OS). It drops GSE22216 (DRFS), GSE59829 (DMFS) and GSE78870 (TTP).

## 6. Per-cohort fits for hsa-miR-29a (`percohort_cox.csv`)

Columns n, nevent, HR, lo, hi, p, schoenfeld_global_p, schoenfeld_mirna_p, covariate_set.

| line | cohort | endpoint | model | covariates | n | events | HR (95 % CI) | p | Schoenfeld global p | Schoenfeld miRNA p |
|---|---|---|---|---|---|---|---|---|---|---|
| 62 | TCGA_BRCA | OS (primary) | unadjusted | none | 1053 | 147 | 0.856 (0.733–1.000) | 0.0499 | 0.588 | 0.588 |
| 63 | TCGA_BRCA | OS (primary) | adjusted | age; stage I-IV; ER (IHC); PAM50; stromal score | 747 | 103 | 0.806 (0.651–0.998) | 0.0474 | 0.0693 | 0.126 |
| 64 | TCGA_BRCA | OS (primary) | adjusted_without_stromal | age; stage I-IV; ER (IHC); PAM50 | 747 | 103 | 0.820 (0.666–1.009) | 0.0605 | 0.0515 | 0.145 |
| 116 | GSE19783 | OS (primary) | unadjusted | none | 99 | 43 | 0.960 (0.713–1.294) | 0.791 | 0.229 | 0.229 |
| 117 | GSE19783 | OS (primary) | adjusted | ER; PAM50 subtype; stromal score | 92 | 42 | 0.935 (0.679–1.288) | 0.681 | 0.0117 | 0.196 |
| 118 | GSE19783 | OS (primary) | adjusted_without_stromal | ER; PAM50 subtype | 92 | 42 | 0.932 (0.678–1.282) | 0.667 | 0.0159 | 0.195 |
| 149 | GSE22216 | DRFS (primary) | unadjusted | none | 210 | 79 | 0.780 (0.656–0.927) | 0.00486 | 0.737 | 0.737 |
| 150 | GSE22216 | DRFS (primary) | adjusted | age; grade; size; nodes; ER; stromal score | 181 | 68 | 0.720 (0.589–0.881) | 0.00139 | 0.0136 | 0.976 |
| 151 | GSE22216 | DRFS (primary) | adjusted_without_stromal | age; grade; size; nodes; ER | 186 | 70 | 0.744 (0.608–0.912) | 0.00436 | 0.00867 | 0.907 |
| 192 | GSE37405 | RFS (primary) | unadjusted | none | 150 | 75 | 0.833 (0.652–1.064) | 0.144 | 0.875 | 0.875 |
| 193 | GSE37405 | RFS (primary) | adjusted | age; size; nodes (ER-positive cohort; platform stratum) | 149 | 75 | 0.811 (0.631–1.041) | 0.1 | 0.69 | 0.978 |
| 194 | GSE37405 | OS (OS) | unadjusted | none | 129 | 70 | 0.904 (0.697–1.171) | 0.443 | 0.153 | 0.153 |
| 195 | GSE37405 | OS (OS) | adjusted | age; size; nodes (ER-positive cohort; platform stratum) | 128 | 70 | 0.898 (0.693–1.166) | 0.42 | 0.483 | 0.131 |
| 216 | GSE59829 | DMFS (primary) | unadjusted | none | 123 | 59 | 0.875 (0.683–1.122) | 0.292 | 0.832 | 0.832 |
| 217 | GSE59829 | DMFS (primary) | adjusted | age; size; ER | 122 | 58 | 0.903 (0.702–1.161) | 0.425 | 0.974 | 0.804 |
| 232 | GSE78870 | TTP (primary) | unadjusted | none | 97 | 80 | 0.710 (0.551–0.914) | 0.00799 | 0.296 | 0.296 |
| 233 | GSE78870 | TTP (primary) | adjusted | age; grade; ER | 84 | 70 | 0.809 (0.624–1.048) | 0.108 | 0.337 | 0.715 |

For TCGA-BRCA and GSE19783 the primary endpoint is OS, so their OS rows repeat the primary rows.

## 7. Proportional hazards (Schoenfeld, `percohort_cox.csv` columns schoenfeld_mirna_p and schoenfeld_global_p)

- Distinct fits (the duplicated OS rows of TCGA-BRCA and GSE19783 counted once): 172.
- miRNA term p < 0.05: 6. Global p < 0.05: 46.
- The global test includes the clinical covariates, which carry most of the departures.
- miRNA-term departures: TCGA_BRCA miR-429 OS adjusted_without_stromal (0.0431); GSE19783 miR-101 OS unadjusted (0.0328); GSE19783 miR-101 OS adjusted (0.0451); GSE19783 miR-101 OS adjusted_without_stromal (0.0453); GSE37405 miR-101 OS unadjusted (0.0439); GSE37405 miR-101 OS adjusted (0.034).

