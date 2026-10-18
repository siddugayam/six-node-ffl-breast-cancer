# VALUES: dynamics of the higher-order modules (Table 4, Fig. 5)

Lines count the header as line 1. Written by `behaviour_intervals_and_prcc.py` (a, b, d), `run_excursions.py` and
`summarise_excursions.py` (c). n = 16,384 Sobol points per module.

## (a) Exact 95 % Clopper–Pearson intervals (`behaviour_frequencies_clopper_pearson.csv`)

Columns k, n, pct, cp95_lo_pct, cp95_hi_pct. Flag 'uncertified (as stored)' is the per-point flag of the sweep
(Table 4); 'certified limit cycle' rows use the 1,000-tau certification.

| line | family | module | behaviour | k | % | 95 % CI (%) | flag |
|---|---|---|---|---|---|---|---|
| 2 | I1_miRNA_FFL | n3 | is_memory | 5970 | 36.438 | 35.700–37.180 | uncertified (as stored) |
| 3 | I1_miRNA_FFL | n3 | is_pulse | 2279 | 13.910 | 13.383–14.449 | uncertified (as stored) |
| 4 | I1_miRNA_FFL | n3 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 5 | I1_miRNA_FFL | n3 | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 6 | I1_miRNA_FFL | n3 | is_bistable | 78 | 0.476 | 0.376–0.594 | uncertified (as stored) |
| 7 | I1_miRNA_FFL | n3 | bistability_unresolved | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 8 | I1_miRNA_FFL | n3 | is_ultrasensitive | 14236 | 86.890 | 86.363–87.403 | uncertified (as stored) |
| 9 | I1_miRNA_FFL | n3 | is_very_ultrasensitive | 5888 | 35.938 | 35.202–36.678 | uncertified (as stored) |
| 10 | I1_miRNA_FFL | n3 | is_noise_rejecting | 13008 | 79.395 | 78.767–80.012 | uncertified (as stored) |
| 11 | I1_miRNA_FFL | n3 | nonmonotone_dose | 7098 | 43.323 | 42.562–44.086 | uncertified (as stored) |
| 12 | I1_miRNA_FFL | n4 | is_memory | 11219 | 68.475 | 67.758–69.186 | uncertified (as stored) |
| 13 | I1_miRNA_FFL | n4 | is_pulse | 1110 | 6.775 | 6.395–7.171 | uncertified (as stored) |
| 14 | I1_miRNA_FFL | n4 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 15 | I1_miRNA_FFL | n4 | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 16 | I1_miRNA_FFL | n4 | is_bistable | 392 | 2.393 | 2.164–2.638 | uncertified (as stored) |
| 17 | I1_miRNA_FFL | n4 | bistability_unresolved | 986 | 6.018 | 5.659–6.393 | uncertified (as stored) |
| 18 | I1_miRNA_FFL | n4 | is_ultrasensitive | 11790 | 71.960 | 71.266–72.647 | uncertified (as stored) |
| 19 | I1_miRNA_FFL | n4 | is_very_ultrasensitive | 3932 | 23.999 | 23.347–24.661 | uncertified (as stored) |
| 20 | I1_miRNA_FFL | n4 | is_noise_rejecting | 13331 | 81.366 | 80.761–81.960 | uncertified (as stored) |
| 21 | I1_miRNA_FFL | n4 | nonmonotone_dose | 4558 | 27.820 | 27.135–28.513 | uncertified (as stored) |
| 22 | I1_miRNA_FFL | n5 | is_memory | 10210 | 62.317 | 61.570–63.060 | uncertified (as stored) |
| 23 | I1_miRNA_FFL | n5 | is_pulse | 1844 | 11.255 | 10.775–11.749 | uncertified (as stored) |
| 24 | I1_miRNA_FFL | n5 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 25 | I1_miRNA_FFL | n5 | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 26 | I1_miRNA_FFL | n5 | is_bistable | 393 | 2.399 | 2.170–2.645 | uncertified (as stored) |
| 27 | I1_miRNA_FFL | n5 | bistability_unresolved | 906 | 5.530 | 5.185–5.891 | uncertified (as stored) |
| 28 | I1_miRNA_FFL | n5 | is_ultrasensitive | 12383 | 75.580 | 74.914–76.236 | uncertified (as stored) |
| 29 | I1_miRNA_FFL | n5 | is_very_ultrasensitive | 5116 | 31.226 | 30.516–31.942 | uncertified (as stored) |
| 30 | I1_miRNA_FFL | n5 | is_noise_rejecting | 12570 | 76.721 | 76.066–77.367 | uncertified (as stored) |
| 31 | I1_miRNA_FFL | n5 | nonmonotone_dose | 6592 | 40.234 | 39.482–40.990 | uncertified (as stored) |
| 32 | I1_miRNA_FFL | n6 | is_memory | 9722 | 59.338 | 58.582–60.092 | uncertified (as stored) |
| 33 | I1_miRNA_FFL | n6 | is_pulse | 2675 | 16.327 | 15.764–16.902 | uncertified (as stored) |
| 34 | I1_miRNA_FFL | n6 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 35 | I1_miRNA_FFL | n6 | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 36 | I1_miRNA_FFL | n6 | is_bistable | 369 | 2.252 | 2.030–2.491 | uncertified (as stored) |
| 37 | I1_miRNA_FFL | n6 | bistability_unresolved | 891 | 5.438 | 5.096–5.797 | uncertified (as stored) |
| 38 | I1_miRNA_FFL | n6 | is_ultrasensitive | 12675 | 77.362 | 76.713–78.001 | uncertified (as stored) |
| 39 | I1_miRNA_FFL | n6 | is_very_ultrasensitive | 5287 | 32.269 | 31.554–32.991 | uncertified (as stored) |
| 40 | I1_miRNA_FFL | n6 | is_noise_rejecting | 13536 | 82.617 | 82.028–83.195 | uncertified (as stored) |
| 41 | I1_miRNA_FFL | n6 | nonmonotone_dose | 5408 | 33.008 | 32.288–33.734 | uncertified (as stored) |
| 42 | I1_miRNA_FFL | n4d | is_memory | 6049 | 36.920 | 36.181–37.664 | uncertified (as stored) |
| 43 | I1_miRNA_FFL | n4d | is_pulse | 1614 | 9.851 | 9.399–10.318 | uncertified (as stored) |
| 44 | I1_miRNA_FFL | n4d | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 45 | I1_miRNA_FFL | n4d | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 46 | I1_miRNA_FFL | n4d | is_bistable | 643 | 3.925 | 3.632–4.233 | uncertified (as stored) |
| 47 | I1_miRNA_FFL | n4d | bistability_unresolved | 24 | 0.146 | 0.094–0.218 | uncertified (as stored) |
| 48 | I1_miRNA_FFL | n4d | is_ultrasensitive | 14369 | 87.701 | 87.189–88.201 | uncertified (as stored) |
| 49 | I1_miRNA_FFL | n4d | is_very_ultrasensitive | 5668 | 34.595 | 33.866–35.329 | uncertified (as stored) |
| 50 | I1_miRNA_FFL | n4d | is_noise_rejecting | 13618 | 83.118 | 82.535–83.688 | uncertified (as stored) |
| 51 | I1_miRNA_FFL | n4d | nonmonotone_dose | 5739 | 35.028 | 34.297–35.764 | uncertified (as stored) |
| 52 | COMP_C2_toggle | n3 | is_memory | 4061 | 24.786 | 24.127–25.455 | uncertified (as stored) |
| 53 | COMP_C2_toggle | n3 | is_pulse | 3 | 0.018 | 0.004–0.054 | uncertified (as stored) |
| 54 | COMP_C2_toggle | n3 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 55 | COMP_C2_toggle | n3 | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 56 | COMP_C2_toggle | n3 | is_bistable | 4130 | 25.208 | 24.544–25.880 | uncertified (as stored) |
| 57 | COMP_C2_toggle | n3 | bistability_unresolved | 8 | 0.049 | 0.021–0.096 | uncertified (as stored) |
| 58 | COMP_C2_toggle | n3 | is_ultrasensitive | 9980 | 60.913 | 60.161–61.661 | uncertified (as stored) |
| 59 | COMP_C2_toggle | n3 | is_very_ultrasensitive | 3207 | 19.574 | 18.969–20.190 | uncertified (as stored) |
| 60 | COMP_C2_toggle | n3 | is_noise_rejecting | 13962 | 85.217 | 84.665–85.758 | uncertified (as stored) |
| 61 | COMP_C2_toggle | n3 | nonmonotone_dose | 907 | 5.536 | 5.191–5.897 | uncertified (as stored) |
| 62 | COMP_C2_toggle | n4 | is_memory | 9684 | 59.106 | 58.349–59.861 | uncertified (as stored) |
| 63 | COMP_C2_toggle | n4 | is_pulse | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 64 | COMP_C2_toggle | n4 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 65 | COMP_C2_toggle | n4 | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 66 | COMP_C2_toggle | n4 | is_bistable | 3741 | 22.833 | 22.192–23.484 | uncertified (as stored) |
| 67 | COMP_C2_toggle | n4 | bistability_unresolved | 637 | 3.888 | 3.597–4.195 | uncertified (as stored) |
| 68 | COMP_C2_toggle | n4 | is_ultrasensitive | 8460 | 51.636 | 50.867–52.404 | uncertified (as stored) |
| 69 | COMP_C2_toggle | n4 | is_very_ultrasensitive | 1627 | 9.930 | 9.477–10.399 | uncertified (as stored) |
| 70 | COMP_C2_toggle | n4 | is_noise_rejecting | 13182 | 80.457 | 79.841–81.061 | uncertified (as stored) |
| 71 | COMP_C2_toggle | n4 | nonmonotone_dose | 885 | 5.402 | 5.060–5.759 | uncertified (as stored) |
| 72 | COMP_C2_toggle | n5 | is_memory | 7954 | 48.547 | 47.779–49.316 | uncertified (as stored) |
| 73 | COMP_C2_toggle | n5 | is_pulse | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 74 | COMP_C2_toggle | n5 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 75 | COMP_C2_toggle | n5 | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 76 | COMP_C2_toggle | n5 | is_bistable | 4272 | 26.074 | 25.403–26.754 | uncertified (as stored) |
| 77 | COMP_C2_toggle | n5 | bistability_unresolved | 513 | 3.131 | 2.870–3.409 | uncertified (as stored) |
| 78 | COMP_C2_toggle | n5 | is_ultrasensitive | 6981 | 42.609 | 41.850–43.370 | uncertified (as stored) |
| 79 | COMP_C2_toggle | n5 | is_very_ultrasensitive | 1220 | 7.446 | 7.049–7.859 | uncertified (as stored) |
| 80 | COMP_C2_toggle | n5 | is_noise_rejecting | 11903 | 72.650 | 71.961–73.332 | uncertified (as stored) |
| 81 | COMP_C2_toggle | n5 | nonmonotone_dose | 1033 | 6.305 | 5.937–6.688 | uncertified (as stored) |
| 82 | COMP_C2_toggle | n6 | is_memory | 9593 | 58.551 | 57.792–59.307 | uncertified (as stored) |
| 83 | COMP_C2_toggle | n6 | is_pulse | 4116 | 25.122 | 24.459–25.794 | uncertified (as stored) |
| 84 | COMP_C2_toggle | n6 | is_sustained_osc | 108 | 0.659 | 0.541–0.795 | uncertified (as stored) |
| 85 | COMP_C2_toggle | n6 | is_damped_osc | 862 | 5.261 | 4.924–5.614 | uncertified (as stored) |
| 86 | COMP_C2_toggle | n6 | is_bistable | 2046 | 12.488 | 11.985–13.004 | uncertified (as stored) |
| 87 | COMP_C2_toggle | n6 | bistability_unresolved | 1419 | 8.661 | 8.235–9.102 | uncertified (as stored) |
| 88 | COMP_C2_toggle | n6 | is_ultrasensitive | 10838 | 66.150 | 65.419–66.875 | uncertified (as stored) |
| 89 | COMP_C2_toggle | n6 | is_very_ultrasensitive | 2854 | 17.419 | 16.841–18.009 | uncertified (as stored) |
| 90 | COMP_C2_toggle | n6 | is_noise_rejecting | 13144 | 80.225 | 79.606–80.832 | uncertified (as stored) |
| 91 | COMP_C2_toggle | n6 | nonmonotone_dose | 2657 | 16.217 | 15.656–16.790 | uncertified (as stored) |
| 92 | COMP_C2_toggle | n4d | is_memory | 4266 | 26.038 | 25.367–26.717 | uncertified (as stored) |
| 93 | COMP_C2_toggle | n4d | is_pulse | 5 | 0.031 | 0.010–0.071 | uncertified (as stored) |
| 94 | COMP_C2_toggle | n4d | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 95 | COMP_C2_toggle | n4d | is_damped_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 96 | COMP_C2_toggle | n4d | is_bistable | 4532 | 27.661 | 26.977–28.353 | uncertified (as stored) |
| 97 | COMP_C2_toggle | n4d | bistability_unresolved | 11 | 0.067 | 0.034–0.120 | uncertified (as stored) |
| 98 | COMP_C2_toggle | n4d | is_ultrasensitive | 10089 | 61.578 | 60.829–62.324 | uncertified (as stored) |
| 99 | COMP_C2_toggle | n4d | is_very_ultrasensitive | 3543 | 21.625 | 20.996–22.263 | uncertified (as stored) |
| 100 | COMP_C2_toggle | n4d | is_noise_rejecting | 13958 | 85.193 | 84.640–85.733 | uncertified (as stored) |
| 101 | COMP_C2_toggle | n4d | nonmonotone_dose | 908 | 5.542 | 5.197–5.903 | uncertified (as stored) |
| 102 | COMP_I1_negfeedback | n3 | is_memory | 4833 | 29.498 | 28.801–30.203 | uncertified (as stored) |
| 103 | COMP_I1_negfeedback | n3 | is_pulse | 4505 | 27.496 | 26.814–28.187 | uncertified (as stored) |
| 104 | COMP_I1_negfeedback | n3 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 105 | COMP_I1_negfeedback | n3 | is_damped_osc | 418 | 2.551 | 2.315–2.804 | uncertified (as stored) |
| 106 | COMP_I1_negfeedback | n3 | is_bistable | 246 | 1.501 | 1.321–1.700 | uncertified (as stored) |
| 107 | COMP_I1_negfeedback | n3 | bistability_unresolved | 1 | 0.006 | 0.000–0.034 | uncertified (as stored) |
| 108 | COMP_I1_negfeedback | n3 | is_ultrasensitive | 13591 | 82.953 | 82.368–83.526 | uncertified (as stored) |
| 109 | COMP_I1_negfeedback | n3 | is_very_ultrasensitive | 5525 | 33.722 | 32.998–34.452 | uncertified (as stored) |
| 110 | COMP_I1_negfeedback | n3 | is_noise_rejecting | 11358 | 69.324 | 68.611–70.029 | uncertified (as stored) |
| 111 | COMP_I1_negfeedback | n3 | nonmonotone_dose | 6946 | 42.395 | 41.637–43.156 | uncertified (as stored) |
| 112 | COMP_I1_negfeedback | n4 | is_memory | 9874 | 60.266 | 59.512–61.017 | uncertified (as stored) |
| 113 | COMP_I1_negfeedback | n4 | is_pulse | 3251 | 19.843 | 19.234–20.462 | uncertified (as stored) |
| 114 | COMP_I1_negfeedback | n4 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 115 | COMP_I1_negfeedback | n4 | is_damped_osc | 164 | 1.001 | 0.854–1.165 | uncertified (as stored) |
| 116 | COMP_I1_negfeedback | n4 | is_bistable | 565 | 3.448 | 3.174–3.739 | uncertified (as stored) |
| 117 | COMP_I1_negfeedback | n4 | bistability_unresolved | 831 | 5.072 | 4.741–5.419 | uncertified (as stored) |
| 118 | COMP_I1_negfeedback | n4 | is_ultrasensitive | 11649 | 71.100 | 70.399–71.793 | uncertified (as stored) |
| 119 | COMP_I1_negfeedback | n4 | is_very_ultrasensitive | 3917 | 23.907 | 23.256–24.568 | uncertified (as stored) |
| 120 | COMP_I1_negfeedback | n4 | is_noise_rejecting | 11827 | 72.186 | 71.493–72.871 | uncertified (as stored) |
| 121 | COMP_I1_negfeedback | n4 | nonmonotone_dose | 4678 | 28.552 | 27.861–29.251 | uncertified (as stored) |
| 122 | COMP_I1_negfeedback | n5 | is_memory | 8840 | 53.955 | 53.188–54.721 | uncertified (as stored) |
| 123 | COMP_I1_negfeedback | n5 | is_pulse | 4204 | 25.659 | 24.992–26.335 | uncertified (as stored) |
| 124 | COMP_I1_negfeedback | n5 | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 125 | COMP_I1_negfeedback | n5 | is_damped_osc | 41 | 0.250 | 0.180–0.339 | uncertified (as stored) |
| 126 | COMP_I1_negfeedback | n5 | is_bistable | 503 | 3.070 | 2.811–3.346 | uncertified (as stored) |
| 127 | COMP_I1_negfeedback | n5 | bistability_unresolved | 782 | 4.773 | 4.452–5.111 | uncertified (as stored) |
| 128 | COMP_I1_negfeedback | n5 | is_ultrasensitive | 11918 | 72.742 | 72.053–73.423 | uncertified (as stored) |
| 129 | COMP_I1_negfeedback | n5 | is_very_ultrasensitive | 4652 | 28.394 | 27.704–29.091 | uncertified (as stored) |
| 130 | COMP_I1_negfeedback | n5 | is_noise_rejecting | 10423 | 63.617 | 62.875–64.354 | uncertified (as stored) |
| 131 | COMP_I1_negfeedback | n5 | nonmonotone_dose | 6353 | 38.776 | 38.029–39.527 | uncertified (as stored) |
| 132 | COMP_I1_negfeedback | n6 | is_memory | 8278 | 50.525 | 49.756–51.293 | uncertified (as stored) |
| 133 | COMP_I1_negfeedback | n6 | is_pulse | 5040 | 30.762 | 30.055–31.475 | uncertified (as stored) |
| 134 | COMP_I1_negfeedback | n6 | is_sustained_osc | 37 | 0.226 | 0.159–0.311 | uncertified (as stored) |
| 135 | COMP_I1_negfeedback | n6 | is_damped_osc | 635 | 3.876 | 3.585–4.183 | uncertified (as stored) |
| 136 | COMP_I1_negfeedback | n6 | is_bistable | 426 | 2.600 | 2.362–2.855 | uncertified (as stored) |
| 137 | COMP_I1_negfeedback | n6 | bistability_unresolved | 824 | 5.029 | 4.700–5.375 | uncertified (as stored) |
| 138 | COMP_I1_negfeedback | n6 | is_ultrasensitive | 12528 | 76.465 | 75.808–77.113 | uncertified (as stored) |
| 139 | COMP_I1_negfeedback | n6 | is_very_ultrasensitive | 5276 | 32.202 | 31.487–32.924 | uncertified (as stored) |
| 140 | COMP_I1_negfeedback | n6 | is_noise_rejecting | 11977 | 73.102 | 72.416–73.780 | uncertified (as stored) |
| 141 | COMP_I1_negfeedback | n6 | nonmonotone_dose | 5210 | 31.799 | 31.086–32.519 | uncertified (as stored) |
| 142 | COMP_I1_negfeedback | n4d | is_memory | 5151 | 31.439 | 30.729–32.156 | uncertified (as stored) |
| 143 | COMP_I1_negfeedback | n4d | is_pulse | 3709 | 22.638 | 21.999–23.287 | uncertified (as stored) |
| 144 | COMP_I1_negfeedback | n4d | is_sustained_osc | 0 | 0.000 | 0.000–0.023 | uncertified (as stored) |
| 145 | COMP_I1_negfeedback | n4d | is_damped_osc | 163 | 0.995 | 0.849–1.159 | uncertified (as stored) |
| 146 | COMP_I1_negfeedback | n4d | is_bistable | 1097 | 6.696 | 6.317–7.089 | uncertified (as stored) |
| 147 | COMP_I1_negfeedback | n4d | bistability_unresolved | 20 | 0.122 | 0.075–0.188 | uncertified (as stored) |
| 148 | COMP_I1_negfeedback | n4d | is_ultrasensitive | 13730 | 83.801 | 83.228–84.362 | uncertified (as stored) |
| 149 | COMP_I1_negfeedback | n4d | is_very_ultrasensitive | 5396 | 32.935 | 32.215–33.660 | uncertified (as stored) |
| 150 | COMP_I1_negfeedback | n4d | is_noise_rejecting | 12052 | 73.560 | 72.877–74.234 | uncertified (as stored) |
| 151 | COMP_I1_negfeedback | n4d | nonmonotone_dose | 5854 | 35.730 | 34.996–36.469 | uncertified (as stored) |
| 152 | COMP_C2_toggle | n6 | is_sustained_osc | 83 | 0.507 | 0.404–0.628 | certified limit cycle |
| 153 | COMP_I1_negfeedback | n6 | is_sustained_osc | 31 | 0.189 | 0.129–0.268 | certified limit cycle |

## (b) PRCC, top 5 per module and behaviour (`prcc_top5.csv`; all parameters in `prcc_binary.csv`)

Columns k (points with the behaviour), rank, parameter, prcc, p, q_bh. p values below about 1e-300 print as 0.

| line | family | module | behaviour | k | rank | parameter | PRCC | p | q (BH) |
|---|---|---|---|---|---|---|---|---|---|
| 2 | I1_miRNA_FFL | n3 | memory | 5970 | 1 | gp | -0.611 | 0 | 0 |
| 3 | I1_miRNA_FFL | n3 | memory | 5970 | 2 | gx | -0.276 | 5.95e-285 | 4.76e-284 |
| 4 | I1_miRNA_FFL | n3 | memory | 5970 | 3 | Kxy | -0.255 | 1.88e-240 | 1e-239 |
| 5 | I1_miRNA_FFL | n3 | memory | 5970 | 4 | Kxm | +0.086 | 3.42e-28 | 1.37e-27 |
| 6 | I1_miRNA_FFL | n3 | memory | 5970 | 5 | Kmp | -0.060 | 1.09e-14 | 3.49e-14 |
| 7 | I1_miRNA_FFL | n3 | pulse | 2279 | 1 | KS | +0.254 | 1.66e-238 | 2.66e-237 |
| 8 | I1_miRNA_FFL | n3 | pulse | 2279 | 2 | theta | -0.190 | 2.04e-133 | 1.63e-132 |
| 9 | I1_miRNA_FFL | n3 | pulse | 2279 | 3 | lam_y | +0.188 | 2.66e-130 | 1.42e-129 |
| 10 | I1_miRNA_FFL | n3 | pulse | 2279 | 4 | gp | +0.179 | 5.54e-118 | 2.22e-117 |
| 11 | I1_miRNA_FFL | n3 | pulse | 2279 | 5 | Kmy | -0.143 | 7.33e-76 | 2.34e-75 |
| 12 | I1_miRNA_FFL | n4 | memory | 11219 | 1 | gp | -0.458 | 0 | 0 |
| 13 | I1_miRNA_FFL | n4 | memory | 11219 | 2 | kappa | +0.321 | 0 | 0 |
| 14 | I1_miRNA_FFL | n4 | memory | 11219 | 3 | Kxy | -0.295 | 0 | 0 |
| 15 | I1_miRNA_FFL | n4 | memory | 11219 | 4 | Kg12 | -0.154 | 1.11e-87 | 5.82e-87 |
| 16 | I1_miRNA_FFL | n4 | memory | 11219 | 5 | gx | -0.147 | 2.19e-79 | 9.22e-79 |
| 17 | I1_miRNA_FFL | n4 | pulse | 1110 | 1 | lam_y | +0.185 | 1.12e-125 | 2.35e-124 |
| 18 | I1_miRNA_FFL | n4 | pulse | 1110 | 2 | KS | +0.174 | 1.28e-111 | 1.35e-110 |
| 19 | I1_miRNA_FFL | n4 | pulse | 1110 | 3 | theta | -0.171 | 1.47e-107 | 1.03e-106 |
| 20 | I1_miRNA_FFL | n4 | pulse | 1110 | 4 | Kmy | -0.126 | 1.11e-58 | 5.82e-58 |
| 21 | I1_miRNA_FFL | n4 | pulse | 1110 | 5 | Kmp | -0.108 | 1.54e-43 | 6.46e-43 |
| 22 | I1_miRNA_FFL | n5 | memory | 10210 | 1 | gp | -0.449 | 0 | 0 |
| 23 | I1_miRNA_FFL | n5 | memory | 10210 | 2 | Kxy | -0.316 | 0 | 0 |
| 24 | I1_miRNA_FFL | n5 | memory | 10210 | 3 | kappa | +0.251 | 7.94e-233 | 6.36e-232 |
| 25 | I1_miRNA_FFL | n5 | memory | 10210 | 4 | Kxm | +0.181 | 3.66e-120 | 2.2e-119 |
| 26 | I1_miRNA_FFL | n5 | memory | 10210 | 5 | gx | -0.164 | 8.08e-99 | 3.88e-98 |
| 27 | I1_miRNA_FFL | n5 | pulse | 1844 | 1 | KS | +0.258 | 5.91e-248 | 1.42e-246 |
| 28 | I1_miRNA_FFL | n5 | pulse | 1844 | 2 | lam_y | +0.152 | 5.62e-85 | 6.74e-84 |
| 29 | I1_miRNA_FFL | n5 | pulse | 1844 | 3 | theta | -0.147 | 3.05e-79 | 2.44e-78 |
| 30 | I1_miRNA_FFL | n5 | pulse | 1844 | 4 | Kc | +0.134 | 2.09e-66 | 1.25e-65 |
| 31 | I1_miRNA_FFL | n5 | pulse | 1844 | 5 | gp | +0.130 | 4.83e-63 | 2.32e-62 |
| 32 | I1_miRNA_FFL | n6 | memory | 9722 | 1 | gp | -0.470 | 0 | 0 |
| 33 | I1_miRNA_FFL | n6 | memory | 9722 | 2 | Kxy | -0.293 | 0 | 0 |
| 34 | I1_miRNA_FFL | n6 | memory | 9722 | 3 | kappa | +0.241 | 5.5e-215 | 5.68e-214 |
| 35 | I1_miRNA_FFL | n6 | memory | 9722 | 4 | Kx2y | -0.189 | 5.9e-132 | 4.58e-131 |
| 36 | I1_miRNA_FFL | n6 | memory | 9722 | 5 | Kg12 | -0.173 | 1.04e-109 | 6.46e-109 |
| 37 | I1_miRNA_FFL | n6 | pulse | 2675 | 1 | KS | +0.236 | 2.26e-205 | 7e-204 |
| 38 | I1_miRNA_FFL | n6 | pulse | 2675 | 2 | lam_y | +0.177 | 1.04e-114 | 1.62e-113 |
| 39 | I1_miRNA_FFL | n6 | pulse | 2675 | 3 | Kc | +0.176 | 1.61e-113 | 1.66e-112 |
| 40 | I1_miRNA_FFL | n6 | pulse | 2675 | 4 | Kmp | -0.160 | 1.46e-94 | 1.13e-93 |
| 41 | I1_miRNA_FFL | n6 | pulse | 2675 | 5 | theta | -0.146 | 6.23e-79 | 3.86e-78 |
| 42 | I1_miRNA_FFL | n4d | memory | 6049 | 1 | gp | -0.633 | 0 | 0 |
| 43 | I1_miRNA_FFL | n4d | memory | 6049 | 2 | gx | -0.280 | 9.78e-293 | 9.29e-292 |
| 44 | I1_miRNA_FFL | n4d | memory | 6049 | 3 | Kxy | -0.238 | 2.18e-210 | 1.38e-209 |
| 45 | I1_miRNA_FFL | n4d | memory | 6049 | 4 | Kxm | +0.071 | 7.09e-20 | 3.37e-19 |
| 46 | I1_miRNA_FFL | n4d | memory | 6049 | 5 | Kmp | -0.057 | 3.72e-13 | 1.41e-12 |
| 47 | I1_miRNA_FFL | n4d | pulse | 1614 | 1 | KS | +0.208 | 7.13e-160 | 1.35e-158 |
| 48 | I1_miRNA_FFL | n4d | pulse | 1614 | 2 | theta | -0.185 | 4.86e-126 | 4.61e-125 |
| 49 | I1_miRNA_FFL | n4d | pulse | 1614 | 3 | lam_y | +0.171 | 4.74e-108 | 3e-107 |
| 50 | I1_miRNA_FFL | n4d | pulse | 1614 | 4 | Kmp | -0.148 | 1.05e-80 | 4.99e-80 |
| 51 | I1_miRNA_FFL | n4d | pulse | 1614 | 5 | gp | +0.124 | 7.14e-57 | 2.71e-56 |
| 52 | COMP_C2_toggle | n3 | memory | 4061 | 1 | gp | -0.654 | 0 | 0 |
| 53 | COMP_C2_toggle | n3 | memory | 4061 | 2 | gx | -0.182 | 3.76e-122 | 3.57e-121 |
| 54 | COMP_C2_toggle | n3 | memory | 4061 | 3 | Kxy | -0.157 | 1.09e-90 | 6.89e-90 |
| 55 | COMP_C2_toggle | n3 | memory | 4061 | 4 | gm | -0.127 | 6.16e-60 | 2.92e-59 |
| 56 | COMP_C2_toggle | n3 | memory | 4061 | 5 | Kxm | -0.123 | 1.28e-56 | 4.85e-56 |
| 57 | COMP_C2_toggle | n4 | memory | 9684 | 1 | gp | -0.447 | 0 | 0 |
| 58 | COMP_C2_toggle | n4 | memory | 9684 | 2 | Kmx | +0.310 | 0 | 0 |
| 59 | COMP_C2_toggle | n4 | memory | 9684 | 3 | Kxm | -0.296 | 0 | 0 |
| 60 | COMP_C2_toggle | n4 | memory | 9684 | 4 | kappa | +0.295 | 0 | 0 |
| 61 | COMP_C2_toggle | n4 | memory | 9684 | 5 | gm | -0.247 | 1.08e-225 | 5.16e-225 |
| 62 | COMP_C2_toggle | n5 | memory | 7954 | 1 | gp | -0.460 | 0 | 0 |
| 63 | COMP_C2_toggle | n5 | memory | 7954 | 2 | Kmx | +0.364 | 0 | 0 |
| 64 | COMP_C2_toggle | n5 | memory | 7954 | 3 | Kxm | -0.338 | 0 | 0 |
| 65 | COMP_C2_toggle | n5 | memory | 7954 | 4 | kappa | +0.252 | 1.64e-234 | 1.11e-233 |
| 66 | COMP_C2_toggle | n5 | memory | 7954 | 5 | Kxy | -0.200 | 6.01e-148 | 3.25e-147 |
| 67 | COMP_C2_toggle | n6 | memory | 9593 | 1 | gp | -0.463 | 0 | 0 |
| 68 | COMP_C2_toggle | n6 | memory | 9593 | 2 | Kxm | -0.313 | 0 | 0 |
| 69 | COMP_C2_toggle | n6 | memory | 9593 | 3 | kappa | +0.297 | 0 | 0 |
| 70 | COMP_C2_toggle | n6 | memory | 9593 | 4 | Kxy | -0.187 | 5.62e-129 | 4.78e-128 |
| 71 | COMP_C2_toggle | n6 | memory | 9593 | 5 | Kx2y | -0.175 | 7.98e-113 | 5.43e-112 |
| 72 | COMP_C2_toggle | n6 | pulse | 4116 | 1 | Kxm | +0.406 | 0 | 0 |
| 73 | COMP_C2_toggle | n6 | pulse | 4116 | 2 | Kmx | -0.309 | 0 | 0 |
| 74 | COMP_C2_toggle | n6 | pulse | 4116 | 3 | Kx2m | -0.205 | 2.61e-155 | 2.96e-154 |
| 75 | COMP_C2_toggle | n6 | pulse | 4116 | 4 | Kxy | +0.137 | 4.97e-69 | 4.22e-68 |
| 76 | COMP_C2_toggle | n6 | pulse | 4116 | 5 | Kx2y | +0.112 | 1.46e-46 | 9.9e-46 |
| 77 | COMP_C2_toggle | n6 | oscillation | 108 | 1 | Kmx | -0.071 | 7.31e-20 | 2.49e-18 |
| 78 | COMP_C2_toggle | n6 | oscillation | 108 | 2 | nmx | +0.055 | 1.34e-12 | 2.28e-11 |
| 79 | COMP_C2_toggle | n6 | oscillation | 108 | 3 | gm | +0.053 | 1.22e-11 | 1.38e-10 |
| 80 | COMP_C2_toggle | n6 | oscillation | 108 | 4 | Kxm | +0.049 | 3.12e-10 | 2.65e-09 |
| 81 | COMP_C2_toggle | n6 | oscillation | 108 | 5 | gx | -0.047 | 1.54e-09 | 9.21e-09 |
| 82 | COMP_C2_toggle | n4d | memory | 4266 | 1 | gp | -0.659 | 0 | 0 |
| 83 | COMP_C2_toggle | n4d | memory | 4266 | 2 | gx | -0.198 | 1.75e-144 | 1.92e-143 |
| 84 | COMP_C2_toggle | n4d | memory | 4266 | 3 | Kxy | -0.177 | 4.44e-115 | 3.26e-114 |
| 85 | COMP_C2_toggle | n4d | memory | 4266 | 4 | Kmx | +0.122 | 3.55e-55 | 1.95e-54 |
| 86 | COMP_C2_toggle | n4d | memory | 4266 | 5 | gm | -0.120 | 1.38e-53 | 6.06e-53 |
| 87 | COMP_I1_negfeedback | n3 | memory | 4833 | 1 | gp | -0.556 | 0 | 0 |
| 88 | COMP_I1_negfeedback | n3 | memory | 4833 | 2 | gx | -0.174 | 6.64e-111 | 6.31e-110 |
| 89 | COMP_I1_negfeedback | n3 | memory | 4833 | 3 | Kxy | -0.164 | 1.86e-98 | 1.18e-97 |
| 90 | COMP_I1_negfeedback | n3 | memory | 4833 | 4 | Kmp | -0.095 | 2.31e-34 | 1.1e-33 |
| 91 | COMP_I1_negfeedback | n3 | memory | 4833 | 5 | Kmx | +0.094 | 3.28e-33 | 1.25e-32 |
| 92 | COMP_I1_negfeedback | n3 | pulse | 4505 | 1 | Kxy | +0.368 | 0 | 0 |
| 93 | COMP_I1_negfeedback | n3 | pulse | 4505 | 2 | Kxm | -0.254 | 1.73e-239 | 1.65e-238 |
| 94 | COMP_I1_negfeedback | n3 | pulse | 4505 | 3 | gp | +0.207 | 3.35e-158 | 2.12e-157 |
| 95 | COMP_I1_negfeedback | n3 | pulse | 4505 | 4 | KS | +0.156 | 5.85e-90 | 2.78e-89 |
| 96 | COMP_I1_negfeedback | n3 | pulse | 4505 | 5 | theta | -0.147 | 9.12e-80 | 3.46e-79 |
| 97 | COMP_I1_negfeedback | n4 | memory | 9874 | 1 | gp | -0.429 | 0 | 0 |
| 98 | COMP_I1_negfeedback | n4 | memory | 9874 | 2 | Kxy | -0.347 | 0 | 0 |
| 99 | COMP_I1_negfeedback | n4 | memory | 9874 | 3 | kappa | +0.282 | 1.76e-297 | 1.41e-296 |
| 100 | COMP_I1_negfeedback | n4 | memory | 9874 | 4 | Kxm | +0.219 | 1.63e-176 | 9.81e-176 |
| 101 | COMP_I1_negfeedback | n4 | memory | 9874 | 5 | gm | -0.178 | 3.78e-116 | 1.82e-115 |
| 102 | COMP_I1_negfeedback | n4 | pulse | 3251 | 1 | Kxy | +0.370 | 0 | 0 |
| 103 | COMP_I1_negfeedback | n4 | pulse | 3251 | 2 | Kxm | -0.272 | 2.68e-275 | 3.22e-274 |
| 104 | COMP_I1_negfeedback | n4 | pulse | 3251 | 3 | Kmx | -0.165 | 1.05e-99 | 8.4e-99 |
| 105 | COMP_I1_negfeedback | n4 | pulse | 3251 | 4 | nxy | +0.162 | 6.9e-97 | 4.14e-96 |
| 106 | COMP_I1_negfeedback | n4 | pulse | 3251 | 5 | gp | +0.152 | 4.8e-85 | 2.3e-84 |
| 107 | COMP_I1_negfeedback | n5 | memory | 8840 | 1 | gp | -0.405 | 0 | 0 |
| 108 | COMP_I1_negfeedback | n5 | memory | 8840 | 2 | Kxy | -0.334 | 0 | 0 |
| 109 | COMP_I1_negfeedback | n5 | memory | 8840 | 3 | Kxm | +0.224 | 2.42e-184 | 2.17e-183 |
| 110 | COMP_I1_negfeedback | n5 | memory | 8840 | 4 | kappa | +0.223 | 7.28e-183 | 4.92e-182 |
| 111 | COMP_I1_negfeedback | n5 | memory | 8840 | 5 | rho_g2 | +0.121 | 2.52e-54 | 1.36e-53 |
| 112 | COMP_I1_negfeedback | n5 | pulse | 4204 | 1 | Kxy | +0.366 | 0 | 0 |
| 113 | COMP_I1_negfeedback | n5 | pulse | 4204 | 2 | Kxm | -0.265 | 1.06e-261 | 1.43e-260 |
| 114 | COMP_I1_negfeedback | n5 | pulse | 4204 | 3 | gp | +0.183 | 1.4e-123 | 1.26e-122 |
| 115 | COMP_I1_negfeedback | n5 | pulse | 4204 | 4 | nxy | +0.178 | 2.49e-116 | 1.68e-115 |
| 116 | COMP_I1_negfeedback | n5 | pulse | 4204 | 5 | KS | +0.152 | 8.12e-85 | 4.39e-84 |
| 117 | COMP_I1_negfeedback | n6 | memory | 8278 | 1 | gp | -0.442 | 0 | 0 |
| 118 | COMP_I1_negfeedback | n6 | memory | 8278 | 2 | Kxy | -0.287 | 2.84e-307 | 4.83e-306 |
| 119 | COMP_I1_negfeedback | n6 | memory | 8278 | 3 | Kx2y | -0.220 | 5.74e-179 | 6.5e-178 |
| 120 | COMP_I1_negfeedback | n6 | memory | 8278 | 4 | kappa | +0.204 | 1.27e-152 | 1.08e-151 |
| 121 | COMP_I1_negfeedback | n6 | memory | 8278 | 5 | Kxm | +0.203 | 2.43e-152 | 1.66e-151 |
| 122 | COMP_I1_negfeedback | n6 | pulse | 5040 | 1 | Kxy | +0.303 | 0 | 0 |
| 123 | COMP_I1_negfeedback | n6 | pulse | 5040 | 2 | Kxm | -0.274 | 1.4e-279 | 2.39e-278 |
| 124 | COMP_I1_negfeedback | n6 | pulse | 5040 | 3 | Kx2y | +0.188 | 5.81e-130 | 6.58e-129 |
| 125 | COMP_I1_negfeedback | n6 | pulse | 5040 | 4 | gp | +0.182 | 2.33e-121 | 1.98e-120 |
| 126 | COMP_I1_negfeedback | n6 | pulse | 5040 | 5 | Kx2m | -0.157 | 4.53e-91 | 3.08e-90 |
| 127 | COMP_I1_negfeedback | n6 | oscillation | 37 | 1 | Kxm | -0.048 | 1.19e-09 | 4.06e-08 |
| 128 | COMP_I1_negfeedback | n6 | oscillation | 37 | 2 | Kmx | -0.046 | 3.22e-09 | 5.47e-08 |
| 129 | COMP_I1_negfeedback | n6 | oscillation | 37 | 3 | nx12 | +0.043 | 3.18e-08 | 3.61e-07 |
| 130 | COMP_I1_negfeedback | n6 | oscillation | 37 | 4 | gx | -0.039 | 6.93e-07 | 5.89e-06 |
| 131 | COMP_I1_negfeedback | n6 | oscillation | 37 | 5 | nmx | +0.038 | 8.77e-07 | 5.96e-06 |
| 132 | COMP_I1_negfeedback | n4d | memory | 5151 | 1 | gp | -0.585 | 0 | 0 |
| 133 | COMP_I1_negfeedback | n4d | memory | 5151 | 2 | gx | -0.196 | 8.95e-142 | 9.85e-141 |
| 134 | COMP_I1_negfeedback | n4d | memory | 5151 | 3 | Kxy | -0.176 | 1.22e-113 | 8.96e-113 |
| 135 | COMP_I1_negfeedback | n4d | memory | 5151 | 4 | Kmx | +0.083 | 3.75e-26 | 2.06e-25 |
| 136 | COMP_I1_negfeedback | n4d | memory | 5151 | 5 | Kmp | -0.081 | 2.18e-25 | 9.58e-25 |
| 137 | COMP_I1_negfeedback | n4d | pulse | 3709 | 1 | Kxy | +0.372 | 0 | 0 |
| 138 | COMP_I1_negfeedback | n4d | pulse | 3709 | 2 | Kxm | -0.272 | 3.09e-276 | 3.4e-275 |
| 139 | COMP_I1_negfeedback | n4d | pulse | 3709 | 3 | gp | +0.163 | 6.93e-98 | 5.08e-97 |
| 140 | COMP_I1_negfeedback | n4d | pulse | 3709 | 4 | nxy | +0.155 | 7.48e-89 | 4.11e-88 |
| 141 | COMP_I1_negfeedback | n4d | pulse | 3709 | 5 | KS | +0.139 | 4.43e-71 | 1.95e-70 |

Not estimable (fewer than 10 points with or without the behaviour; `prcc_binary.csv` column note):

- line 34: I1_miRNA_FFL n3 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 77: I1_miRNA_FFL n4 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 126: I1_miRNA_FFL n5 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 189: I1_miRNA_FFL n6 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 228: I1_miRNA_FFL n4d oscillation: not estimable (0 of 16384 points with the behaviour)
- line 248: COMP_C2_toggle n3 pulse: not estimable (3 of 16384 points with the behaviour)
- line 249: COMP_C2_toggle n3 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 274: COMP_C2_toggle n4 pulse: not estimable (0 of 16384 points with the behaviour)
- line 275: COMP_C2_toggle n4 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 303: COMP_C2_toggle n5 pulse: not estimable (0 of 16384 points with the behaviour)
- line 304: COMP_C2_toggle n5 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 429: COMP_C2_toggle n4d pulse: not estimable (5 of 16384 points with the behaviour)
- line 430: COMP_C2_toggle n4d oscillation: not estimable (0 of 16384 points with the behaviour)
- line 469: COMP_I1_negfeedback n3 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 518: COMP_I1_negfeedback n4 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 573: COMP_I1_negfeedback n5 oscillation: not estimable (0 of 16384 points with the behaviour)
- line 720: COMP_I1_negfeedback n4d oscillation: not estimable (0 of 16384 points with the behaviour)

Parameters that act in each module: `parameters_acting_by_module.csv` (column n_acting):

- line 2: I1_miRNA_FFL n3: 16 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;Kmx;nmx;lam_x;b;rho;phi;rho_g2;kappa;w2;Kg12;ng12
- line 3: I1_miRNA_FFL n4: 21 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;Kmx;nmx;lam_x;b;rho;phi
- line 4: I1_miRNA_FFL n5: 24 acting; not acting: gx2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;Kmx;nmx;lam_x;b
- line 5: I1_miRNA_FFL n6: 31 acting; not acting: Kmx;nmx;lam_x;b
- line 6: I1_miRNA_FFL n4d: 19 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;Kmx;nmx;lam_x;b;rho;phi;rho_g2;kappa
- line 7: COMP_C2_toggle n3: 19 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b;rho;phi;rho_g2;kappa;w2;Kg12;ng12
- line 8: COMP_C2_toggle n4: 24 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b;rho;phi
- line 9: COMP_C2_toggle n5: 27 acting; not acting: gx2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b
- line 10: COMP_C2_toggle n6: 34 acting; not acting: b
- line 11: COMP_C2_toggle n4d: 22 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b;rho;phi;rho_g2;kappa
- line 12: COMP_I1_negfeedback n3: 19 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b;rho;phi;rho_g2;kappa;w2;Kg12;ng12
- line 13: COMP_I1_negfeedback n4: 24 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b;rho;phi
- line 14: COMP_I1_negfeedback n5: 27 acting; not acting: gx2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b
- line 15: COMP_I1_negfeedback n6: 34 acting; not acting: b
- line 16: COMP_I1_negfeedback n4d: 22 acting; not acting: gx2;gm2;Kx2y;nx2y;Kx2m;nx2m;Kx12;nx12;b;rho;phi;rho_g2;kappa

The continuous-proxy PRCCs (secondary) are in `prcc_continuous.csv` (columns proxy, n_points, parameter, prcc, p,
q_bh, rank_abs).

## (c) One-at-a-time excursions (`excursion_frequencies.csv`)

Columns baseline_k, baseline_pct, k, pct, cp95_lo_pct, cp95_hi_pct, gained, lost, mcnemar_p, status. Rows for memory,
pulse and sustained oscillation are listed here; every behaviour is in the file.

| line | family | module | parameter | factor | behaviour | baseline % | excursion % (95 % CI) | gained / lost | McNemar p | status |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | COMP_C2_toggle | n6 | gx | x0.5 | is_memory | 58.55 | 61.07 (60.32–61.82) | 731 / 318 | 4.93e-38 | run |
| 3 | COMP_C2_toggle | n6 | gx | x0.5 | is_pulse | 25.12 | 26.35 (25.68–27.03) | 334 / 133 | 5.21e-21 | run |
| 4 | COMP_C2_toggle | n6 | gx | x0.5 | is_sustained_osc | 0.66 | 1.18 (1.02–1.36) | 153 / 67 | 6.48e-09 | run |
| 12 | COMP_C2_toggle | n6 | gx | x2.0 | is_memory | 58.55 | 56.70 (55.93–57.46) | 178 / 482 | 2.88e-33 | run |
| 13 | COMP_C2_toggle | n6 | gx | x2.0 | is_pulse | 25.12 | 24.43 (23.78–25.10) | 85 / 198 | 1.5e-11 | run |
| 14 | COMP_C2_toggle | n6 | gx | x2.0 | is_sustained_osc | 0.66 | 0.37 (0.28–0.48) | 51 / 98 | 0.000146 | run |
| 22 | COMP_C2_toggle | n6 | gm | x0.5 | is_memory | 58.55 | 61.52 (60.77–62.27) | 622 / 135 | 1.77e-75 | run |
| 23 | COMP_C2_toggle | n6 | gm | x0.5 | is_pulse | 25.12 | 23.19 (22.55–23.85) | 197 / 513 | 2.49e-33 | run |
| 24 | COMP_C2_toggle | n6 | gm | x0.5 | is_sustained_osc | 0.66 | 0.40 (0.31–0.51) | 50 / 92 | 0.000534 | run |
| 32 | COMP_C2_toggle | n6 | gm | x2.0 | is_memory | 58.55 | 55.41 (54.65–56.18) | 138 / 652 | 1.25e-80 | run |
| 33 | COMP_C2_toggle | n6 | gm | x2.0 | is_pulse | 25.12 | 26.51 (25.83–27.19) | 511 / 284 | 7.02e-16 | run |
| 34 | COMP_C2_toggle | n6 | gm | x2.0 | is_sustained_osc | 0.66 | 1.15 (1.00–1.33) | 148 / 67 | 3.39e-08 | run |
| 42 | COMP_C2_toggle | n6 | gp | x0.5 | is_memory | 58.55 | 73.51 (72.83–74.19) | 2640 / 189 | 0 | run |
| 43 | COMP_C2_toggle | n6 | gp | x0.5 | is_pulse | 25.12 | 22.25 (21.62–22.90) | 12 / 482 | 1.55e-125 | run |
| 44 | COMP_C2_toggle | n6 | gp | x0.5 | is_sustained_osc | 0.66 | 0.63 (0.52–0.77) | 3 / 7 | 0.344 | run |
| 52 | COMP_C2_toggle | n6 | gp | x2.0 | is_memory | 58.55 | 40.63 (39.88–41.39) | 116 / 3052 | 0 | run |
| 53 | COMP_C2_toggle | n6 | gp | x2.0 | is_pulse | 25.12 | 27.33 (26.64–28.01) | 373 / 12 | 4.88e-94 | run |
| 54 | COMP_C2_toggle | n6 | gp | x2.0 | is_sustained_osc | 0.66 | 0.72 (0.60–0.86) | 10 / 0 | 0.00195 | run |
| 62 | COMP_C2_toggle | n6 | Kxm | x0.5 | is_memory | 58.55 | 65.83 (65.10–66.56) | 1379 / 186 | 3.73e-225 | run |
| 63 | COMP_C2_toggle | n6 | Kxm | x0.5 | is_pulse | 25.12 | 16.23 (15.67–16.80) | 140 / 1597 | 2.48e-313 | run |
| 64 | COMP_C2_toggle | n6 | Kxm | x0.5 | is_sustained_osc | 0.66 | 0.52 (0.41–0.64) | 60 / 83 | 0.0654 | run |
| 72 | COMP_C2_toggle | n6 | Kxm | x2.0 | is_memory | 58.55 | 51.67 (50.90–52.44) | 191 / 1318 | 3.31e-207 | run |
| 73 | COMP_C2_toggle | n6 | Kxm | x2.0 | is_pulse | 25.12 | 33.59 (32.87–34.32) | 1569 / 181 | 6.19e-276 | run |
| 74 | COMP_C2_toggle | n6 | Kxm | x2.0 | is_sustained_osc | 0.66 | 0.93 (0.79–1.09) | 102 / 58 | 0.000629 | run |
| 82 | COMP_C2_toggle | n6 | Kxy | x0.5 | is_memory | 58.55 | 63.27 (62.53–64.01) | 936 / 163 | 2.74e-132 | run |
| 83 | COMP_C2_toggle | n6 | Kxy | x0.5 | is_pulse | 25.12 | 21.60 (20.97–22.24) | 194 / 771 | 5.68e-82 | run |
| 84 | COMP_C2_toggle | n6 | Kxy | x0.5 | is_sustained_osc | 0.66 | 0.70 (0.58–0.84) | 49 / 42 | 0.53 | run |
| 92 | COMP_C2_toggle | n6 | Kxy | x2.0 | is_memory | 58.55 | 51.72 (50.95–52.49) | 136 / 1255 | 3.83e-227 | run |
| 93 | COMP_C2_toggle | n6 | Kxy | x2.0 | is_pulse | 25.12 | 28.85 (28.15–29.55) | 870 / 260 | 2.78e-77 | run |
| 94 | COMP_C2_toggle | n6 | Kxy | x2.0 | is_sustained_osc | 0.66 | 0.80 (0.67–0.95) | 57 / 34 | 0.0206 | run |
| 102 | COMP_C2_toggle | n6 | Kx2y | x0.5 | is_memory | 58.55 | 62.67 (61.93–63.41) | 853 / 178 | 3.57e-106 | run |
| 103 | COMP_C2_toggle | n6 | Kx2y | x0.5 | is_pulse | 25.12 | 22.37 (21.73–23.02) | 195 / 646 | 3.89e-57 | run |
| 104 | COMP_C2_toggle | n6 | Kx2y | x0.5 | is_sustained_osc | 0.66 | 0.66 (0.54–0.80) | 41 / 41 | 1 | run |
| 112 | COMP_C2_toggle | n6 | Kx2y | x2.0 | is_memory | 58.55 | 51.76 (50.99–52.53) | 134 / 1247 | 1.65e-226 | run |
| 113 | COMP_C2_toggle | n6 | Kx2y | x2.0 | is_pulse | 25.12 | 28.55 (27.86–29.25) | 816 / 254 | 3.1e-69 | run |
| 114 | COMP_C2_toggle | n6 | Kx2y | x2.0 | is_sustained_osc | 0.66 | 0.82 (0.69–0.97) | 57 / 30 | 0.00501 | run |
| 122 | COMP_C2_toggle | n6 | Kx2m | x0.5 | is_memory | 58.55 | 54.95 (54.18–55.71) | 208 / 798 | 5.57e-82 | run |
| 123 | COMP_C2_toggle | n6 | Kx2m | x0.5 | is_pulse | 25.12 | 29.94 (29.24–30.65) | 1192 / 403 | 1.75e-90 | run |
| 124 | COMP_C2_toggle | n6 | Kx2m | x0.5 | is_sustained_osc | 0.66 | 0.99 (0.84–1.15) | 128 / 74 | 0.000177 | run |
| 132 | COMP_C2_toggle | n6 | Kx2m | x2.0 | is_memory | 58.55 | 63.28 (62.54–64.02) | 975 / 200 | 9.77e-123 | run |
| 133 | COMP_C2_toggle | n6 | Kx2m | x2.0 | is_pulse | 25.12 | 19.07 (18.47–19.68) | 304 / 1296 | 1.09e-145 | run |
| 134 | COMP_C2_toggle | n6 | Kx2m | x2.0 | is_sustained_osc | 0.66 | 0.46 (0.37–0.58) | 55 / 87 | 0.00904 | run |
| 142 | COMP_C2_toggle | n6 | Kmx | x0.5 | is_memory | 58.55 | 53.58 (52.82–54.35) | 240 / 1054 | 7.14e-122 | run |
| 143 | COMP_C2_toggle | n6 | Kmx | x0.5 | is_pulse | 25.12 | 32.93 (32.21–33.66) | 1593 / 313 | 5.06e-206 | run |
| 144 | COMP_C2_toggle | n6 | Kmx | x0.5 | is_sustained_osc | 0.66 | 1.35 (1.18–1.54) | 189 / 76 | 2.8e-12 | run |
| 152 | COMP_C2_toggle | n6 | Kmx | x2.0 | is_memory | 58.55 | 62.59 (61.85–63.33) | 857 / 195 | 1.86e-99 | run |
| 153 | COMP_C2_toggle | n6 | Kmx | x2.0 | is_pulse | 25.12 | 17.72 (17.14–18.32) | 229 / 1441 | 1.02e-214 | run |
| 154 | COMP_C2_toggle | n6 | Kmx | x2.0 | is_sustained_osc | 0.66 | 0.34 (0.25–0.44) | 38 / 91 | 3.46e-06 | run |
| 162 | COMP_C2_toggle | n6 | nmx | x0.5 | is_memory | 58.55 | 54.47 (53.71–55.24) | 348 / 1016 | 4.84e-76 | run |
| 163 | COMP_C2_toggle | n6 | nmx | x0.5 | is_pulse | 25.12 | 30.40 (29.70–31.11) | 1279 / 414 | 1.48e-102 | run |
| 164 | COMP_C2_toggle | n6 | nmx | x0.5 | is_sustained_osc | 0.66 | 0.08 (0.04–0.14) | 12 / 107 | 3.21e-20 | run |
| 172 | COMP_C2_toggle | n6 | nmx | x2.0 | is_memory | 58.55 | 61.02 (60.27–61.76) | 644 / 240 | 2.3e-43 | run |
| 173 | COMP_C2_toggle | n6 | nmx | x2.0 | is_pulse | 25.12 | 21.28 (20.65–21.91) | 264 / 894 | 1.91e-80 | run |
| 174 | COMP_C2_toggle | n6 | nmx | x2.0 | is_sustained_osc | 0.66 | 1.81 (1.61–2.02) | 268 / 80 | 7.63e-25 | run |
| 182 | COMP_C2_toggle | n6 | kappa | x0.5 | is_memory | 58.55 | 48.36 (47.60–49.13) | 90 / 1759 | 0 | run |
| 183 | COMP_C2_toggle | n6 | kappa | x0.5 | is_pulse | 25.12 | 25.38 (24.72–26.06) | 45 / 2 | 1.6e-11 | run |
| 184 | COMP_C2_toggle | n6 | kappa | x0.5 | is_sustained_osc | 0.66 | 0.67 (0.55–0.80) | 2 / 1 | 1 | run |
| 192 | COMP_C2_toggle | n6 | kappa | x2.0 | is_memory | 58.55 | 66.38 (65.65–67.10) | 1428 / 146 | 1.6e-264 | run |
| 193 | COMP_C2_toggle | n6 | kappa | x2.0 | is_pulse | 25.12 | 24.64 (23.98–25.31) | 3 / 82 | 5.3e-21 | run |
| 194 | COMP_C2_toggle | n6 | kappa | x2.0 | is_sustained_osc | 0.66 | 0.66 (0.54–0.80) | 0 / 0 | 1 | run |
| 202 | COMP_I1_negfeedback | n6 | gx | x0.5 | is_memory | 50.52 | 54.00 (53.23–54.76) | 652 / 83 | 2.08e-110 | run |
| 203 | COMP_I1_negfeedback | n6 | gx | x0.5 | is_pulse | 30.76 | 31.27 (30.56–31.99) | 202 / 118 | 3.08e-06 | run |
| 204 | COMP_I1_negfeedback | n6 | gx | x0.5 | is_sustained_osc | 0.23 | 0.28 (0.21–0.37) | 38 / 29 | 0.328 | run |
| 212 | COMP_I1_negfeedback | n6 | gx | x2.0 | is_memory | 50.52 | 48.09 (47.32–48.86) | 40 / 439 | 5.29e-86 | run |
| 213 | COMP_I1_negfeedback | n6 | gx | x2.0 | is_pulse | 30.76 | 30.29 (29.59–31.00) | 74 / 151 | 3.14e-07 | run |
| 214 | COMP_I1_negfeedback | n6 | gx | x2.0 | is_sustained_osc | 0.23 | 0.05 (0.03–0.10) | 5 / 33 | 4.26e-06 | run |
| 222 | COMP_I1_negfeedback | n6 | gp | x0.5 | is_memory | 50.52 | 66.31 (65.58–67.03) | 2703 / 117 | 0 | run |
| 223 | COMP_I1_negfeedback | n6 | gp | x0.5 | is_pulse | 30.76 | 24.59 (23.93–25.26) | 51 / 1062 | 8.94e-247 | run |
| 224 | COMP_I1_negfeedback | n6 | gp | x0.5 | is_sustained_osc | 0.23 | 0.18 (0.12–0.25) | 0 / 8 | 0.00781 | run |
| 232 | COMP_I1_negfeedback | n6 | gp | x2.0 | is_memory | 50.52 | 34.41 (33.68–35.14) | 87 / 2728 | 0 | run |
| 233 | COMP_I1_negfeedback | n6 | gp | x2.0 | is_pulse | 30.76 | 36.20 (35.46–36.94) | 927 / 36 | 9.52e-225 | run |
| 234 | COMP_I1_negfeedback | n6 | gp | x2.0 | is_sustained_osc | 0.23 | 0.27 (0.20–0.37) | 8 / 0 | 0.00781 | run |
| 242 | COMP_I1_negfeedback | n6 | Kxm | x0.5 | is_memory | 50.52 | 46.01 (45.24–46.78) | 406 / 1146 | 1.07e-81 | run |
| 243 | COMP_I1_negfeedback | n6 | Kxm | x0.5 | is_pulse | 30.76 | 37.28 (36.54–38.03) | 1334 / 266 | 6.21e-171 | run |
| 244 | COMP_I1_negfeedback | n6 | Kxm | x0.5 | is_sustained_osc | 0.23 | 0.35 (0.26–0.45) | 40 / 20 | 0.0135 | run |
| 252 | COMP_I1_negfeedback | n6 | Kxm | x2.0 | is_memory | 50.52 | 56.48 (55.72–57.24) | 1373 / 397 | 1.42e-125 | run |
| 253 | COMP_I1_negfeedback | n6 | Kxm | x2.0 | is_pulse | 30.76 | 22.99 (22.34–23.64) | 256 / 1530 | 8.79e-221 | run |
| 254 | COMP_I1_negfeedback | n6 | Kxm | x2.0 | is_sustained_osc | 0.23 | 0.12 (0.07–0.18) | 13 / 31 | 0.00956 | run |
| 262 | COMP_I1_negfeedback | n6 | Kxy | x0.5 | is_memory | 50.52 | 57.23 (56.47–57.99) | 1399 / 300 | 2.11e-169 | run |
| 263 | COMP_I1_negfeedback | n6 | Kxy | x0.5 | is_pulse | 30.76 | 23.51 (22.86–24.17) | 173 / 1361 | 6.04e-229 | run |
| 264 | COMP_I1_negfeedback | n6 | Kxy | x0.5 | is_sustained_osc | 0.23 | 0.20 (0.14–0.28) | 12 / 16 | 0.572 | run |
| 272 | COMP_I1_negfeedback | n6 | Kxy | x2.0 | is_memory | 50.52 | 42.18 (41.42–42.94) | 295 / 1663 | 6.28e-231 | run |
| 273 | COMP_I1_negfeedback | n6 | Kxy | x2.0 | is_pulse | 30.76 | 39.39 (38.64–40.14) | 1688 / 275 | 1.78e-247 | run |
| 274 | COMP_I1_negfeedback | n6 | Kxy | x2.0 | is_sustained_osc | 0.23 | 0.19 (0.13–0.27) | 14 / 20 | 0.392 | run |
| 282 | COMP_I1_negfeedback | n6 | Kx2y | x0.5 | is_memory | 50.52 | 55.94 (55.18–56.71) | 1099 / 211 | 4.64e-145 | run |
| 283 | COMP_I1_negfeedback | n6 | Kx2y | x0.5 | is_pulse | 30.76 | 26.60 (25.92–27.28) | 258 / 940 | 2.03e-91 | run |
| 284 | COMP_I1_negfeedback | n6 | Kx2y | x0.5 | is_sustained_osc | 0.23 | 0.20 (0.14–0.28) | 14 / 18 | 0.597 | run |
| 292 | COMP_I1_negfeedback | n6 | Kx2y | x2.0 | is_memory | 50.52 | 43.50 (42.74–44.26) | 252 / 1403 | 1.79e-193 | run |
| 293 | COMP_I1_negfeedback | n6 | Kx2y | x2.0 | is_pulse | 30.76 | 35.86 (35.12–36.60) | 1173 / 338 | 6.19e-108 | run |
| 294 | COMP_I1_negfeedback | n6 | Kx2y | x2.0 | is_sustained_osc | 0.23 | 0.24 (0.17–0.33) | 21 / 19 | 0.875 | run |
| 302 | COMP_I1_negfeedback | n6 | Kx2m | x0.5 | is_memory | 50.52 | 46.23 (45.46–47.00) | 265 / 969 | 1.53e-94 | run |
| 303 | COMP_I1_negfeedback | n6 | Kx2m | x0.5 | is_pulse | 30.76 | 34.46 (33.73–35.19) | 1031 / 425 | 2.63e-58 | run |
| 304 | COMP_I1_negfeedback | n6 | Kx2m | x0.5 | is_sustained_osc | 0.23 | 0.27 (0.20–0.36) | 33 / 26 | 0.435 | run |
| 312 | COMP_I1_negfeedback | n6 | Kx2m | x2.0 | is_memory | 50.52 | 56.32 (55.55–57.08) | 1228 / 279 | 5.5e-142 | run |
| 313 | COMP_I1_negfeedback | n6 | Kx2m | x2.0 | is_pulse | 30.76 | 24.25 (23.59–24.91) | 330 / 1397 | 6.05e-156 | run |
| 314 | COMP_I1_negfeedback | n6 | Kx2m | x2.0 | is_sustained_osc | 0.23 | 0.16 (0.10–0.23) | 15 / 26 | 0.117 | run |
| 322 | COMP_I1_negfeedback | n6 | nx12 | x0.5 | is_memory | 50.52 | 51.23 (50.46–51.99) | 637 / 522 | 0.000805 | run |
| 323 | COMP_I1_negfeedback | n6 | nx12 | x0.5 | is_pulse | 30.76 | 29.14 (28.44–29.84) | 471 / 737 | 1.93e-14 | run |
| 324 | COMP_I1_negfeedback | n6 | nx12 | x0.5 | is_sustained_osc | 0.23 | 0.01 (0.00–0.04) | 2 / 37 | 2.84e-09 | run |
| 332 | COMP_I1_negfeedback | n6 | nx12 | x2.0 | is_memory | 50.52 | 49.76 (48.99–50.52) | 498 / 624 | 0.000187 | run |
| 333 | COMP_I1_negfeedback | n6 | nx12 | x2.0 | is_pulse | 30.76 | 31.77 (31.06–32.49) | 697 / 531 | 2.41e-06 | run |
| 334 | COMP_I1_negfeedback | n6 | nx12 | x2.0 | is_sustained_osc | 0.23 | 1.12 (0.96–1.29) | 175 / 29 | 1.24e-26 | run |
| 342 | COMP_I1_negfeedback | n6 | Kmx | x0.5 | is_memory | 50.52 | 47.14 (46.37–47.91) | 442 / 997 | 1.6e-49 | run |
| 343 | COMP_I1_negfeedback | n6 | Kmx | x0.5 | is_pulse | 30.76 | 33.11 (32.39–33.84) | 996 / 611 | 6.26e-22 | run |
| 344 | COMP_I1_negfeedback | n6 | Kmx | x0.5 | is_sustained_osc | 0.23 | 0.31 (0.23–0.41) | 38 / 24 | 0.098 | run |
| 352 | COMP_I1_negfeedback | n6 | Kmx | x2.0 | is_memory | 50.52 | 53.11 (52.34–53.87) | 802 / 379 | 2.14e-35 | run |
| 353 | COMP_I1_negfeedback | n6 | Kmx | x2.0 | is_pulse | 30.76 | 27.67 (26.98–28.36) | 520 / 1027 | 1.52e-38 | run |
| 354 | COMP_I1_negfeedback | n6 | Kmx | x2.0 | is_sustained_osc | 0.23 | 0.12 (0.07–0.18) | 11 / 29 | 0.00643 | run |
| 362 | COMP_I1_negfeedback | n6 | nmx | x0.5 | is_memory | 50.52 | 48.43 (47.66–49.20) | 489 / 832 | 2.99e-21 | run |
| 363 | COMP_I1_negfeedback | n6 | nmx | x0.5 | is_pulse | 30.76 | 31.08 (30.37–31.79) | 621 / 569 | 0.139 | run |
| 364 | COMP_I1_negfeedback | n6 | nmx | x0.5 | is_sustained_osc | 0.23 | 0.01 (0.00–0.04) | 2 / 37 | 2.84e-09 | run |
| 372 | COMP_I1_negfeedback | n6 | nmx | x2.0 | is_memory | 50.52 | 51.58 (50.81–52.35) | 503 / 330 | 2.23e-09 | run |
| 373 | COMP_I1_negfeedback | n6 | nmx | x2.0 | is_pulse | 30.76 | 29.71 (29.01–30.42) | 385 / 557 | 2.32e-08 | run |
| 374 | COMP_I1_negfeedback | n6 | nmx | x2.0 | is_sustained_osc | 0.23 | 0.71 (0.59–0.85) | 106 / 27 | 2.91e-12 | run |
| 382 | COMP_I1_negfeedback | n6 | kappa | x0.5 | is_memory | 50.52 | 43.12 (42.36–43.88) | 68 / 1282 | 4.54e-291 | run |
| 383 | COMP_I1_negfeedback | n6 | kappa | x0.5 | is_pulse | 30.76 | 31.54 (30.83–32.25) | 155 / 28 | 1.65e-22 | run |
| 384 | COMP_I1_negfeedback | n6 | kappa | x0.5 | is_sustained_osc | 0.23 | 0.23 (0.16–0.31) | 0 / 0 | 1 | run |
| 392 | COMP_I1_negfeedback | n6 | kappa | x2.0 | is_memory | 50.52 | 57.03 (56.27–57.79) | 1143 / 77 | 3.15e-244 | run |
| 393 | COMP_I1_negfeedback | n6 | kappa | x2.0 | is_pulse | 30.76 | 29.82 (29.12–30.52) | 29 / 184 | 9.21e-29 | run |
| 394 | COMP_I1_negfeedback | n6 | kappa | x2.0 | is_sustained_osc | 0.23 | 0.23 (0.16–0.31) | 0 / 0 | 1 | run |
| 402 | I1_miRNA_FFL | n6 | gp | x0.5 | is_memory | 59.34 | 74.41 (73.73–75.07) | 2594 / 125 | 0 | run |
| 403 | I1_miRNA_FFL | n6 | gp | x0.5 | is_pulse | 16.33 | 12.56 (12.06–13.08) | 44 / 661 | 2.54e-142 | run |
| 404 | I1_miRNA_FFL | n6 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 412 | I1_miRNA_FFL | n6 | gp | x2.0 | is_memory | 59.34 | 42.30 (41.54–43.06) | 78 / 2870 | 0 | run |
| 413 | I1_miRNA_FFL | n6 | gp | x2.0 | is_pulse | 16.33 | 20.10 (19.49–20.73) | 654 / 35 | 7.2e-149 | run |
| 414 | I1_miRNA_FFL | n6 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 422 | I1_miRNA_FFL | n6 | KS | x0.5 | is_memory | 59.34 | 60.42 (59.67–61.17) | 993 / 816 | 3.46e-05 | run |
| 423 | I1_miRNA_FFL | n6 | KS | x0.5 | is_pulse | 16.33 | 11.63 (11.14–12.13) | 287 / 1057 | 9.71e-104 | run |
| 424 | I1_miRNA_FFL | n6 | KS | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 432 | I1_miRNA_FFL | n6 | KS | x2.0 | is_memory | 59.34 | 56.31 (55.55–57.07) | 673 / 1169 | 3.92e-31 | run |
| 433 | I1_miRNA_FFL | n6 | KS | x2.0 | is_pulse | 16.33 | 20.45 (19.83–21.07) | 1119 / 444 | 2.74e-67 | run |
| 434 | I1_miRNA_FFL | n6 | KS | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 442 | I1_miRNA_FFL | n6 | Kxy | x0.5 | is_memory | 59.34 | 66.09 (65.36–66.82) | 1407 / 300 | 3.87e-171 | run |
| 443 | I1_miRNA_FFL | n6 | Kxy | x0.5 | is_pulse | 16.33 | 14.42 (13.88–14.96) | 231 / 544 | 6.1e-30 | run |
| 444 | I1_miRNA_FFL | n6 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 452 | I1_miRNA_FFL | n6 | Kxy | x2.0 | is_memory | 59.34 | 50.49 (49.72–51.26) | 276 / 1726 | 6.84e-256 | run |
| 453 | I1_miRNA_FFL | n6 | Kxy | x2.0 | is_pulse | 16.33 | 21.26 (20.64–21.90) | 1063 / 254 | 7.4e-118 | run |
| 454 | I1_miRNA_FFL | n6 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 462 | I1_miRNA_FFL | n6 | Kx2y | x0.5 | is_memory | 59.34 | 64.18 (63.44–64.92) | 1012 / 218 | 1.48e-122 | run |
| 463 | I1_miRNA_FFL | n6 | Kx2y | x0.5 | is_pulse | 16.33 | 14.98 (14.43–15.53) | 267 / 488 | 7.37e-16 | run |
| 464 | I1_miRNA_FFL | n6 | Kx2y | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 472 | I1_miRNA_FFL | n6 | Kx2y | x2.0 | is_memory | 59.34 | 52.82 (52.05–53.59) | 262 / 1330 | 5.79e-172 | run |
| 473 | I1_miRNA_FFL | n6 | Kx2y | x2.0 | is_pulse | 16.33 | 20.21 (19.60–20.84) | 934 / 297 | 5.22e-77 | run |
| 474 | I1_miRNA_FFL | n6 | Kx2y | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 482 | I1_miRNA_FFL | n6 | Kmp | x0.5 | is_memory | 59.34 | 58.50 (57.74–59.26) | 463 / 600 | 2.96e-05 | run |
| 483 | I1_miRNA_FFL | n6 | Kmp | x0.5 | is_pulse | 16.33 | 19.32 (18.72–19.93) | 731 / 241 | 5.84e-58 | run |
| 484 | I1_miRNA_FFL | n6 | Kmp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 492 | I1_miRNA_FFL | n6 | Kmp | x2.0 | is_memory | 59.34 | 60.58 (59.83–61.33) | 531 / 327 | 3.34e-12 | run |
| 493 | I1_miRNA_FFL | n6 | Kmp | x2.0 | is_pulse | 16.33 | 12.90 (12.39–13.42) | 152 / 714 | 9.42e-88 | run |
| 494 | I1_miRNA_FFL | n6 | Kmp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 502 | I1_miRNA_FFL | n6 | lam_y | x0.5 | is_memory | 59.34 | 61.69 (60.94–62.44) | 592 / 206 | 5.01e-44 | run |
| 503 | I1_miRNA_FFL | n6 | lam_y | x0.5 | is_pulse | 16.33 | 13.07 (12.56–13.60) | 101 / 634 | 3.23e-95 | run |
| 504 | I1_miRNA_FFL | n6 | lam_y | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 512 | I1_miRNA_FFL | n6 | lam_y | x2.0 | is_memory | 59.34 | 57.06 (56.29–57.82) | 263 / 637 | 1.59e-36 | run |
| 513 | I1_miRNA_FFL | n6 | lam_y | x2.0 | is_pulse | 16.33 | 20.82 (20.20–21.45) | 859 / 123 | 1.72e-136 | run |
| 514 | I1_miRNA_FFL | n6 | lam_y | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 522 | I1_miRNA_FFL | n6 | theta | x0.5 | is_memory | 59.34 | 58.30 (57.54–59.06) | 347 / 517 | 8.04e-09 | run |
| 523 | I1_miRNA_FFL | n6 | theta | x0.5 | is_pulse | 16.33 | 20.47 (19.86–21.10) | 820 / 141 | 4.82e-117 | run |
| 524 | I1_miRNA_FFL | n6 | theta | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 532 | I1_miRNA_FFL | n6 | theta | x2.0 | is_memory | 59.34 | 60.44 (59.69–61.19) | 515 / 334 | 5.63e-10 | run |
| 533 | I1_miRNA_FFL | n6 | theta | x2.0 | is_pulse | 16.33 | 13.50 (12.98–14.03) | 139 / 602 | 1.84e-69 | run |
| 534 | I1_miRNA_FFL | n6 | theta | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 542 | I1_miRNA_FFL | n6 | Kc | x0.5 | is_memory | 59.34 | 60.38 (59.62–61.13) | 355 / 185 | 2.25e-13 | run |
| 543 | I1_miRNA_FFL | n6 | Kc | x0.5 | is_pulse | 16.33 | 13.62 (13.10–14.16) | 83 / 526 | 1.07e-79 | run |
| 544 | I1_miRNA_FFL | n6 | Kc | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 552 | I1_miRNA_FFL | n6 | Kc | x2.0 | is_memory | 59.34 | 58.84 (58.08–59.59) | 175 / 257 | 9.31e-05 | run |
| 553 | I1_miRNA_FFL | n6 | Kc | x2.0 | is_pulse | 16.33 | 18.36 (17.77–18.96) | 423 / 90 | 1.28e-52 | run |
| 554 | I1_miRNA_FFL | n6 | Kc | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 562 | I1_miRNA_FFL | n6 | kappa | x0.5 | is_memory | 59.34 | 51.49 (50.72–52.26) | 65 / 1351 | 2.08e-313 | run |
| 563 | I1_miRNA_FFL | n6 | kappa | x0.5 | is_pulse | 16.33 | 17.15 (16.58–17.74) | 168 / 33 | 5.58e-23 | run |
| 564 | I1_miRNA_FFL | n6 | kappa | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 572 | I1_miRNA_FFL | n6 | kappa | x2.0 | is_memory | 59.34 | 66.11 (65.38–66.84) | 1218 / 108 | 2.2e-238 | run |
| 573 | I1_miRNA_FFL | n6 | kappa | x2.0 | is_pulse | 16.33 | 15.42 (14.87–15.98) | 37 / 186 | 4.37e-25 | run |
| 574 | I1_miRNA_FFL | n6 | kappa | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 582 | I1_miRNA_FFL | n6 | Kg12 | x0.5 | is_memory | 59.34 | 63.61 (62.87–64.35) | 739 / 39 | 1.38e-168 | run |
| 583 | I1_miRNA_FFL | n6 | Kg12 | x0.5 | is_pulse | 16.33 | 15.56 (15.01–16.13) | 23 / 148 | 1.48e-23 | run |
| 584 | I1_miRNA_FFL | n6 | Kg12 | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 592 | I1_miRNA_FFL | n6 | Kg12 | x2.0 | is_memory | 59.34 | 54.71 (53.95–55.48) | 56 / 814 | 2.58e-173 | run |
| 593 | I1_miRNA_FFL | n6 | Kg12 | x2.0 | is_pulse | 16.33 | 17.12 (16.55–17.71) | 151 / 21 | 1.87e-25 | run |
| 594 | I1_miRNA_FFL | n6 | Kg12 | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 602 | COMP_C2_toggle | n5 | gp | x0.5 | is_memory | 48.55 | 65.78 (65.05–66.50) | 2945 / 122 | 0 | run |
| 603 | COMP_C2_toggle | n5 | gp | x0.5 | is_pulse | 0.00 | 0.01 (0.00–0.03) | 1 / 0 | 1 | run |
| 604 | COMP_C2_toggle | n5 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 612 | COMP_C2_toggle | n5 | gp | x2.0 | is_memory | 48.55 | 31.85 (31.13–32.57) | 83 / 2819 | 0 | run |
| 613 | COMP_C2_toggle | n5 | gp | x2.0 | is_pulse | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 614 | COMP_C2_toggle | n5 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 622 | COMP_C2_toggle | n5 | Kxm | x0.5 | is_memory | 48.55 | 57.14 (56.38–57.90) | 1535 / 127 | 2.64e-307 | run |
| 623 | COMP_C2_toggle | n5 | Kxm | x0.5 | is_pulse | 0.00 | 0.01 (0.00–0.04) | 2 / 0 | 0.5 | run |
| 624 | COMP_C2_toggle | n5 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 632 | COMP_C2_toggle | n5 | Kxm | x2.0 | is_memory | 48.55 | 40.88 (40.13–41.64) | 102 / 1358 | 1.1e-280 | run |
| 633 | COMP_C2_toggle | n5 | Kxm | x2.0 | is_pulse | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 634 | COMP_C2_toggle | n5 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 642 | COMP_C2_toggle | n5 | Kxy | x0.5 | is_memory | 48.55 | 53.38 (52.61–54.14) | 942 / 151 | 3.35e-140 | run |
| 643 | COMP_C2_toggle | n5 | Kxy | x0.5 | is_pulse | 0.00 | 0.01 (0.00–0.03) | 1 / 0 | 1 | run |
| 644 | COMP_C2_toggle | n5 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 652 | COMP_C2_toggle | n5 | Kxy | x2.0 | is_memory | 48.55 | 42.50 (41.74–43.26) | 107 / 1098 | 1.19e-207 | run |
| 653 | COMP_C2_toggle | n5 | Kxy | x2.0 | is_pulse | 0.00 | 0.01 (0.00–0.03) | 1 / 0 | 1 | run |
| 654 | COMP_C2_toggle | n5 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 662 | COMP_C2_toggle | n5 | Kmx | x0.5 | is_memory | 48.55 | 38.60 (37.85–39.35) | 48 / 1678 | 0 | run |
| 663 | COMP_C2_toggle | n5 | Kmx | x0.5 | is_pulse | 0.00 | 0.02 (0.00–0.05) | 3 / 0 | 0.25 | run |
| 664 | COMP_C2_toggle | n5 | Kmx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 672 | COMP_C2_toggle | n5 | Kmx | x2.0 | is_memory | 48.55 | 57.10 (56.34–57.86) | 1451 / 50 | 0 | run |
| 673 | COMP_C2_toggle | n5 | Kmx | x2.0 | is_pulse | 0.00 | 0.01 (0.00–0.04) | 2 / 0 | 0.5 | run |
| 674 | COMP_C2_toggle | n5 | Kmx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 682 | COMP_C2_toggle | n5 | kappa | x0.5 | is_memory | 48.55 | 40.64 (39.88–41.39) | 85 / 1381 | 4.09e-302 | run |
| 683 | COMP_C2_toggle | n5 | kappa | x0.5 | is_pulse | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 684 | COMP_C2_toggle | n5 | kappa | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 692 | COMP_C2_toggle | n5 | kappa | x2.0 | is_memory | 48.55 | 54.36 (53.60–55.13) | 1085 / 132 | 1.02e-186 | run |
| 693 | COMP_C2_toggle | n5 | kappa | x2.0 | is_pulse | 0.00 | 0.01 (0.00–0.03) | 1 / 0 | 1 | run |
| 694 | COMP_C2_toggle | n5 | kappa | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 702 | COMP_I1_negfeedback | n5 | gp | x0.5 | is_memory | 53.96 | 68.71 (68.00–69.42) | 2562 / 144 | 0 | run |
| 703 | COMP_I1_negfeedback | n5 | gp | x0.5 | is_pulse | 25.66 | 19.95 (19.34–20.57) | 49 / 985 | 3.05e-227 | run |
| 704 | COMP_I1_negfeedback | n5 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 712 | COMP_I1_negfeedback | n5 | gp | x2.0 | is_memory | 53.96 | 39.03 (38.28–39.78) | 92 / 2537 | 0 | run |
| 713 | COMP_I1_negfeedback | n5 | gp | x2.0 | is_pulse | 25.66 | 30.67 (29.96–31.38) | 864 / 43 | 1.76e-199 | run |
| 714 | COMP_I1_negfeedback | n5 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 722 | COMP_I1_negfeedback | n5 | KS | x0.5 | is_memory | 53.96 | 53.74 (52.97–54.51) | 944 / 979 | 0.438 | run |
| 723 | COMP_I1_negfeedback | n5 | KS | x0.5 | is_pulse | 25.66 | 22.16 (21.53–22.81) | 367 / 940 | 2.87e-58 | run |
| 724 | COMP_I1_negfeedback | n5 | KS | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 732 | COMP_I1_negfeedback | n5 | KS | x2.0 | is_memory | 53.96 | 53.53 (52.76–54.29) | 968 / 1038 | 0.123 | run |
| 733 | COMP_I1_negfeedback | n5 | KS | x2.0 | is_pulse | 25.66 | 29.00 (28.30–29.70) | 980 / 433 | 5.13e-49 | run |
| 734 | COMP_I1_negfeedback | n5 | KS | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 742 | COMP_I1_negfeedback | n5 | Kxm | x0.5 | is_memory | 53.96 | 48.49 (47.72–49.25) | 726 / 1622 | 4.44e-78 | run |
| 743 | COMP_I1_negfeedback | n5 | Kxm | x0.5 | is_pulse | 25.66 | 31.88 (31.17–32.60) | 1560 / 540 | 3.27e-114 | run |
| 744 | COMP_I1_negfeedback | n5 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 752 | COMP_I1_negfeedback | n5 | Kxm | x2.0 | is_memory | 53.96 | 60.57 (59.82–61.32) | 1770 / 686 | 1.61e-109 | run |
| 753 | COMP_I1_negfeedback | n5 | Kxm | x2.0 | is_pulse | 25.66 | 18.66 (18.07–19.27) | 409 / 1555 | 8.92e-157 | run |
| 754 | COMP_I1_negfeedback | n5 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 762 | COMP_I1_negfeedback | n5 | Kxy | x0.5 | is_memory | 53.96 | 62.46 (61.71–63.20) | 1912 / 519 | 2.49e-186 | run |
| 763 | COMP_I1_negfeedback | n5 | Kxy | x0.5 | is_pulse | 25.66 | 16.74 (16.17–17.32) | 262 / 1723 | 3.9e-263 | run |
| 764 | COMP_I1_negfeedback | n5 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 772 | COMP_I1_negfeedback | n5 | Kxy | x2.0 | is_memory | 53.96 | 45.23 (44.46–45.99) | 486 / 1916 | 1.05e-199 | run |
| 773 | COMP_I1_negfeedback | n5 | Kxy | x2.0 | is_pulse | 25.66 | 34.97 (34.24–35.71) | 1896 / 370 | 4.16e-246 | run |
| 774 | COMP_I1_negfeedback | n5 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 782 | COMP_I1_negfeedback | n5 | nxy | x0.5 | is_memory | 53.96 | 62.18 (61.43–62.93) | 2211 / 863 | 7.81e-135 | run |
| 783 | COMP_I1_negfeedback | n5 | nxy | x0.5 | is_pulse | 25.66 | 13.48 (12.96–14.01) | 438 / 2434 | 0 | run |
| 784 | COMP_I1_negfeedback | n5 | nxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 792 | COMP_I1_negfeedback | n5 | nxy | x2.0 | is_memory | 53.96 | 48.59 (47.82–49.36) | 710 / 1589 | 8.78e-77 | run |
| 793 | COMP_I1_negfeedback | n5 | nxy | x2.0 | is_pulse | 25.66 | 31.78 (31.07–32.50) | 1987 / 984 | 6.25e-77 | run |
| 794 | COMP_I1_negfeedback | n5 | nxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 802 | COMP_I1_negfeedback | n5 | rho_g2 | x0.5 | is_memory | 53.96 | 49.81 (49.04–50.58) | 238 / 917 | 2.33e-94 | run |
| 803 | COMP_I1_negfeedback | n5 | rho_g2 | x0.5 | is_pulse | 25.66 | 28.04 (27.35–28.73) | 472 / 82 | 1.45e-67 | run |
| 804 | COMP_I1_negfeedback | n5 | rho_g2 | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 812 | COMP_I1_negfeedback | n5 | rho_g2 | x2.0 | is_memory | 53.96 | 58.30 (57.54–59.06) | 963 / 251 | 1.33e-98 | run |
| 813 | COMP_I1_negfeedback | n5 | rho_g2 | x2.0 | is_pulse | 25.66 | 22.69 (22.05–23.34) | 92 / 579 | 3.21e-87 | run |
| 814 | COMP_I1_negfeedback | n5 | rho_g2 | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 822 | COMP_I1_negfeedback | n5 | kappa | x0.5 | is_memory | 53.96 | 46.24 (45.47–47.01) | 73 / 1337 | 1.96e-301 | run |
| 823 | COMP_I1_negfeedback | n5 | kappa | x0.5 | is_pulse | 25.66 | 26.37 (25.70–27.06) | 153 / 36 | 2.26e-18 | run |
| 824 | COMP_I1_negfeedback | n5 | kappa | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 832 | COMP_I1_negfeedback | n5 | kappa | x2.0 | is_memory | 53.96 | 60.24 (59.48–60.99) | 1144 / 115 | 1.13e-213 | run |
| 833 | COMP_I1_negfeedback | n5 | kappa | x2.0 | is_pulse | 25.66 | 24.64 (23.98–25.31) | 33 / 200 | 2.45e-30 | run |
| 834 | COMP_I1_negfeedback | n5 | kappa | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 842 | I1_miRNA_FFL | n5 | gx | x0.5 | is_memory | 62.32 | 69.34 (68.62–70.04) | 1245 / 95 | 3.42e-256 | run |
| 843 | I1_miRNA_FFL | n5 | gx | x0.5 | is_pulse | 11.25 | 11.45 (10.97–11.95) | 188 / 156 | 0.0945 | run |
| 844 | I1_miRNA_FFL | n5 | gx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 852 | I1_miRNA_FFL | n5 | gx | x2.0 | is_memory | 62.32 | 57.75 (56.98–58.50) | 101 / 850 | 3.21e-148 | run |
| 853 | I1_miRNA_FFL | n5 | gx | x2.0 | is_pulse | 11.25 | 11.16 (10.68–11.66) | 131 / 146 | 0.4 | run |
| 854 | I1_miRNA_FFL | n5 | gx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 862 | I1_miRNA_FFL | n5 | gp | x0.5 | is_memory | 62.32 | 76.30 (75.64–76.95) | 2435 / 144 | 0 | run |
| 863 | I1_miRNA_FFL | n5 | gp | x0.5 | is_pulse | 11.25 | 8.34 (7.92–8.77) | 54 / 532 | 8.96e-100 | run |
| 864 | I1_miRNA_FFL | n5 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 872 | I1_miRNA_FFL | n5 | gp | x2.0 | is_memory | 62.32 | 45.82 (45.05–46.59) | 88 / 2791 | 0 | run |
| 873 | I1_miRNA_FFL | n5 | gp | x2.0 | is_pulse | 11.25 | 14.18 (13.65–14.72) | 519 / 40 | 2.65e-107 | run |
| 874 | I1_miRNA_FFL | n5 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 882 | I1_miRNA_FFL | n5 | KS | x0.5 | is_memory | 62.32 | 62.16 (61.42–62.91) | 959 / 984 | 0.586 | run |
| 883 | I1_miRNA_FFL | n5 | KS | x0.5 | is_pulse | 11.25 | 7.39 (6.99–7.80) | 272 / 906 | 5.38e-80 | run |
| 884 | I1_miRNA_FFL | n5 | KS | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 892 | I1_miRNA_FFL | n5 | KS | x2.0 | is_memory | 62.32 | 61.32 (60.57–62.07) | 932 / 1095 | 0.000318 | run |
| 893 | I1_miRNA_FFL | n5 | KS | x2.0 | is_pulse | 11.25 | 15.78 (15.23–16.35) | 1095 / 353 | 1.72e-88 | run |
| 894 | I1_miRNA_FFL | n5 | KS | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 902 | I1_miRNA_FFL | n5 | Kxm | x0.5 | is_memory | 62.32 | 58.13 (57.37–58.89) | 458 / 1144 | 9.32e-68 | run |
| 903 | I1_miRNA_FFL | n5 | Kxm | x0.5 | is_pulse | 11.25 | 9.89 (9.43–10.35) | 429 / 653 | 1.02e-11 | run |
| 904 | I1_miRNA_FFL | n5 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 912 | I1_miRNA_FFL | n5 | Kxm | x2.0 | is_memory | 62.32 | 67.35 (66.62–68.06) | 1377 / 553 | 8.41e-81 | run |
| 913 | I1_miRNA_FFL | n5 | Kxm | x2.0 | is_pulse | 11.25 | 10.69 (10.22–11.17) | 603 / 696 | 0.0107 | run |
| 914 | I1_miRNA_FFL | n5 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 922 | I1_miRNA_FFL | n5 | Kxy | x0.5 | is_memory | 62.32 | 69.67 (68.96–70.38) | 1725 / 520 | 4.39e-150 | run |
| 923 | I1_miRNA_FFL | n5 | Kxy | x0.5 | is_pulse | 11.25 | 9.59 (9.14–10.05) | 370 / 643 | 8.11e-18 | run |
| 924 | I1_miRNA_FFL | n5 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 932 | I1_miRNA_FFL | n5 | Kxy | x2.0 | is_memory | 62.32 | 52.84 (52.07–53.61) | 442 / 1995 | 1.68e-234 | run |
| 933 | I1_miRNA_FFL | n5 | Kxy | x2.0 | is_pulse | 11.25 | 15.31 (14.76–15.87) | 1013 / 349 | 3.13e-75 | run |
| 934 | I1_miRNA_FFL | n5 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 942 | I1_miRNA_FFL | n5 | lam_y | x0.5 | is_memory | 62.32 | 64.44 (63.70–65.17) | 650 / 302 | 5.49e-30 | run |
| 943 | I1_miRNA_FFL | n5 | lam_y | x0.5 | is_pulse | 11.25 | 8.58 (8.15–9.01) | 124 / 563 | 1.17e-67 | run |
| 944 | I1_miRNA_FFL | n5 | lam_y | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 952 | I1_miRNA_FFL | n5 | lam_y | x2.0 | is_memory | 62.32 | 59.85 (59.10–60.60) | 364 / 768 | 9.39e-34 | run |
| 953 | I1_miRNA_FFL | n5 | lam_y | x2.0 | is_pulse | 11.25 | 14.94 (14.39–15.49) | 750 / 147 | 4.9e-98 | run |
| 954 | I1_miRNA_FFL | n5 | lam_y | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 962 | I1_miRNA_FFL | n5 | theta | x0.5 | is_memory | 62.32 | 61.06 (60.31–61.81) | 429 / 635 | 2.9e-10 | run |
| 963 | I1_miRNA_FFL | n5 | theta | x0.5 | is_pulse | 11.25 | 14.38 (13.85–14.93) | 657 / 145 | 1.43e-78 | run |
| 964 | I1_miRNA_FFL | n5 | theta | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 972 | I1_miRNA_FFL | n5 | theta | x2.0 | is_memory | 62.32 | 63.75 (63.01–64.49) | 626 / 391 | 1.71e-13 | run |
| 973 | I1_miRNA_FFL | n5 | theta | x2.0 | is_pulse | 11.25 | 8.58 (8.15–9.01) | 126 / 565 | 3.29e-67 | run |
| 974 | I1_miRNA_FFL | n5 | theta | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 982 | I1_miRNA_FFL | n5 | Kc | x0.5 | is_memory | 62.32 | 63.62 (62.87–64.35) | 440 / 227 | 1.24e-16 | run |
| 983 | I1_miRNA_FFL | n5 | Kc | x0.5 | is_pulse | 11.25 | 9.44 (8.99–9.89) | 113 / 411 | 9.32e-41 | run |
| 984 | I1_miRNA_FFL | n5 | Kc | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 992 | I1_miRNA_FFL | n5 | Kc | x2.0 | is_memory | 62.32 | 61.58 (60.83–62.32) | 190 / 311 | 7.19e-08 | run |
| 993 | I1_miRNA_FFL | n5 | Kc | x2.0 | is_pulse | 11.25 | 12.69 (12.18–13.21) | 362 / 127 | 3.42e-27 | run |
| 994 | I1_miRNA_FFL | n5 | Kc | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1002 | I1_miRNA_FFL | n5 | kappa | x0.5 | is_memory | 62.32 | 53.77 (53.00–54.53) | 77 / 1478 | 0 | run |
| 1003 | I1_miRNA_FFL | n5 | kappa | x0.5 | is_pulse | 11.25 | 11.93 (11.43–12.43) | 145 / 35 | 4.13e-17 | run |
| 1004 | I1_miRNA_FFL | n5 | kappa | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1012 | I1_miRNA_FFL | n5 | kappa | x2.0 | is_memory | 62.32 | 68.97 (68.25–69.68) | 1223 / 133 | 4.62e-221 | run |
| 1013 | I1_miRNA_FFL | n5 | kappa | x2.0 | is_pulse | 11.25 | 10.68 (10.21–11.16) | 56 / 151 | 2.89e-11 | run |
| 1014 | I1_miRNA_FFL | n5 | kappa | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1022 | COMP_C2_toggle | n3 | gx | x0.5 | is_memory | 24.79 | 33.46 (32.74–34.19) | 1575 / 154 | 8.8e-297 | run |
| 1023 | COMP_C2_toggle | n3 | gx | x0.5 | is_pulse | 0.02 | 0.01 (0.00–0.04) | 1 / 2 | 1 | run |
| 1024 | COMP_C2_toggle | n3 | gx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1032 | COMP_C2_toggle | n3 | gx | x2.0 | is_memory | 24.79 | 20.16 (19.55–20.78) | 69 / 827 | 8.41e-166 | run |
| 1033 | COMP_C2_toggle | n3 | gx | x2.0 | is_pulse | 0.02 | 0.02 (0.00–0.05) | 2 / 2 | 1 | run |
| 1034 | COMP_C2_toggle | n3 | gx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1042 | COMP_C2_toggle | n3 | gm | x0.5 | is_memory | 24.79 | 28.14 (27.46–28.84) | 689 / 139 | 2.7e-88 | run |
| 1043 | COMP_C2_toggle | n3 | gm | x0.5 | is_pulse | 0.02 | 0.57 (0.46–0.69) | 93 / 3 | 3.72e-24 | run |
| 1044 | COMP_C2_toggle | n3 | gm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1052 | COMP_C2_toggle | n3 | gm | x2.0 | is_memory | 24.79 | 21.54 (20.91–22.18) | 103 / 635 | 2.43e-94 | run |
| 1053 | COMP_C2_toggle | n3 | gm | x2.0 | is_pulse | 0.02 | 0.01 (0.00–0.03) | 1 / 3 | 0.625 | run |
| 1054 | COMP_C2_toggle | n3 | gm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1062 | COMP_C2_toggle | n3 | gp | x0.5 | is_memory | 24.79 | 47.51 (46.74–48.28) | 3763 / 40 | 0 | run |
| 1063 | COMP_C2_toggle | n3 | gp | x0.5 | is_pulse | 0.02 | 0.04 (0.02–0.09) | 5 / 1 | 0.219 | run |
| 1064 | COMP_C2_toggle | n3 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1072 | COMP_C2_toggle | n3 | gp | x2.0 | is_memory | 24.79 | 6.79 (6.41–7.19) | 8 / 2956 | 0 | run |
| 1073 | COMP_C2_toggle | n3 | gp | x2.0 | is_pulse | 0.02 | 0.01 (0.00–0.04) | 0 / 1 | 1 | run |
| 1074 | COMP_C2_toggle | n3 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1082 | COMP_C2_toggle | n3 | Kxm | x0.5 | is_memory | 24.79 | 27.26 (26.58–27.95) | 592 / 186 | 4.23e-50 | run |
| 1083 | COMP_C2_toggle | n3 | Kxm | x0.5 | is_pulse | 0.02 | 0.01 (0.00–0.04) | 2 / 3 | 1 | run |
| 1084 | COMP_C2_toggle | n3 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1092 | COMP_C2_toggle | n3 | Kxm | x2.0 | is_memory | 24.79 | 22.47 (21.83–23.11) | 158 / 538 | 2.51e-49 | run |
| 1093 | COMP_C2_toggle | n3 | Kxm | x2.0 | is_pulse | 0.02 | 0.02 (0.00–0.05) | 2 / 2 | 1 | run |
| 1094 | COMP_C2_toggle | n3 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1102 | COMP_C2_toggle | n3 | Kxy | x0.5 | is_memory | 24.79 | 28.11 (27.43–28.81) | 656 / 111 | 6.63e-95 | run |
| 1103 | COMP_C2_toggle | n3 | Kxy | x0.5 | is_pulse | 0.02 | 0.01 (0.00–0.03) | 0 / 2 | 0.5 | run |
| 1104 | COMP_C2_toggle | n3 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1112 | COMP_C2_toggle | n3 | Kxy | x2.0 | is_memory | 24.79 | 21.72 (21.09–22.36) | 80 / 582 | 5.44e-95 | run |
| 1113 | COMP_C2_toggle | n3 | Kxy | x2.0 | is_pulse | 0.02 | 0.04 (0.01–0.08) | 5 / 2 | 0.453 | run |
| 1114 | COMP_C2_toggle | n3 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1122 | COMP_I1_negfeedback | n3 | gx | x0.5 | is_memory | 29.50 | 37.92 (37.18–38.67) | 1433 / 53 | 0 | run |
| 1123 | COMP_I1_negfeedback | n3 | gx | x0.5 | is_pulse | 27.50 | 28.45 (27.76–29.15) | 288 / 131 | 1.24e-14 | run |
| 1124 | COMP_I1_negfeedback | n3 | gx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1132 | COMP_I1_negfeedback | n3 | gx | x2.0 | is_memory | 29.50 | 24.44 (23.78–25.10) | 40 / 869 | 5.47e-204 | run |
| 1133 | COMP_I1_negfeedback | n3 | gx | x2.0 | is_pulse | 27.50 | 27.14 (26.46–27.82) | 104 / 163 | 0.000367 | run |
| 1134 | COMP_I1_negfeedback | n3 | gx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1142 | COMP_I1_negfeedback | n3 | gp | x0.5 | is_memory | 29.50 | 49.58 (48.82–50.35) | 3347 / 56 | 0 | run |
| 1143 | COMP_I1_negfeedback | n3 | gp | x0.5 | is_pulse | 27.50 | 21.25 (20.62–21.88) | 29 / 1053 | 0 | run |
| 1144 | COMP_I1_negfeedback | n3 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1152 | COMP_I1_negfeedback | n3 | gp | x2.0 | is_memory | 29.50 | 13.04 (12.53–13.57) | 20 / 2716 | 0 | run |
| 1153 | COMP_I1_negfeedback | n3 | gp | x2.0 | is_pulse | 27.50 | 33.26 (32.54–33.99) | 948 / 3 | 1.51e-278 | run |
| 1154 | COMP_I1_negfeedback | n3 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1162 | COMP_I1_negfeedback | n3 | KS | x0.5 | is_memory | 29.50 | 28.43 (27.74–29.13) | 820 / 995 | 4.37e-05 | run |
| 1163 | COMP_I1_negfeedback | n3 | KS | x0.5 | is_pulse | 27.50 | 24.19 (23.54–24.86) | 510 / 1051 | 2.24e-43 | run |
| 1164 | COMP_I1_negfeedback | n3 | KS | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1172 | COMP_I1_negfeedback | n3 | KS | x2.0 | is_memory | 29.50 | 29.88 (29.18–30.59) | 965 / 902 | 0.151 | run |
| 1173 | COMP_I1_negfeedback | n3 | KS | x2.0 | is_pulse | 27.50 | 31.05 (30.34–31.76) | 1128 / 546 | 1.09e-46 | run |
| 1174 | COMP_I1_negfeedback | n3 | KS | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1182 | COMP_I1_negfeedback | n3 | Kxm | x0.5 | is_memory | 29.50 | 27.92 (27.23–28.61) | 712 / 971 | 2.96e-10 | run |
| 1183 | COMP_I1_negfeedback | n3 | Kxm | x0.5 | is_pulse | 27.50 | 33.75 (33.03–34.48) | 1623 / 598 | 1.16e-108 | run |
| 1184 | COMP_I1_negfeedback | n3 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1192 | COMP_I1_negfeedback | n3 | Kxm | x2.0 | is_memory | 29.50 | 31.10 (30.39–31.81) | 975 / 713 | 1.95e-10 | run |
| 1193 | COMP_I1_negfeedback | n3 | Kxm | x2.0 | is_pulse | 27.50 | 21.08 (20.45–21.71) | 537 / 1589 | 4.04e-120 | run |
| 1194 | COMP_I1_negfeedback | n3 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1202 | COMP_I1_negfeedback | n3 | Kxy | x0.5 | is_memory | 29.50 | 33.05 (32.33–33.78) | 1217 / 635 | 3.27e-42 | run |
| 1203 | COMP_I1_negfeedback | n3 | Kxy | x0.5 | is_pulse | 27.50 | 19.08 (18.48–19.69) | 393 / 1772 | 2.63e-208 | run |
| 1204 | COMP_I1_negfeedback | n3 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1212 | COMP_I1_negfeedback | n3 | Kxy | x2.0 | is_memory | 29.50 | 26.14 (25.46–26.82) | 533 / 1084 | 2.08e-43 | run |
| 1213 | COMP_I1_negfeedback | n3 | Kxy | x2.0 | is_pulse | 27.50 | 37.96 (37.21–38.71) | 2092 / 378 | 1.76e-286 | run |
| 1214 | COMP_I1_negfeedback | n3 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1222 | COMP_I1_negfeedback | n3 | Kmx | x0.5 | is_memory | 29.50 | 26.93 (26.26–27.62) | 442 / 862 | 1.19e-31 | run |
| 1223 | COMP_I1_negfeedback | n3 | Kmx | x0.5 | is_pulse | 27.50 | 30.48 (29.78–31.19) | 1177 / 688 | 6.2e-30 | run |
| 1224 | COMP_I1_negfeedback | n3 | Kmx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1232 | COMP_I1_negfeedback | n3 | Kmx | x2.0 | is_memory | 29.50 | 31.87 (31.15–32.59) | 766 / 378 | 6.97e-31 | run |
| 1233 | COMP_I1_negfeedback | n3 | Kmx | x2.0 | is_pulse | 27.50 | 24.04 (23.38–24.70) | 557 / 1124 | 3.61e-44 | run |
| 1234 | COMP_I1_negfeedback | n3 | Kmx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1242 | COMP_I1_negfeedback | n3 | Kmp | x0.5 | is_memory | 29.50 | 32.18 (31.47–32.91) | 822 / 382 | 1.68e-37 | run |
| 1243 | COMP_I1_negfeedback | n3 | Kmp | x0.5 | is_pulse | 27.50 | 27.05 (26.37–27.74) | 486 / 559 | 0.0259 | run |
| 1244 | COMP_I1_negfeedback | n3 | Kmp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1252 | COMP_I1_negfeedback | n3 | Kmp | x2.0 | is_memory | 29.50 | 27.54 (26.86–28.23) | 260 / 581 | 5.47e-29 | run |
| 1253 | COMP_I1_negfeedback | n3 | Kmp | x2.0 | is_pulse | 27.50 | 27.52 (26.84–28.21) | 402 / 398 | 0.916 | run |
| 1254 | COMP_I1_negfeedback | n3 | Kmp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1262 | COMP_I1_negfeedback | n3 | theta | x0.5 | is_memory | 29.50 | 29.23 (28.53–29.93) | 463 / 507 | 0.167 | run |
| 1263 | COMP_I1_negfeedback | n3 | theta | x0.5 | is_pulse | 27.50 | 32.23 (31.52–32.95) | 982 / 206 | 1.72e-121 | run |
| 1264 | COMP_I1_negfeedback | n3 | theta | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1272 | COMP_I1_negfeedback | n3 | theta | x2.0 | is_memory | 29.50 | 30.15 (29.44–30.85) | 523 / 417 | 0.000608 | run |
| 1273 | COMP_I1_negfeedback | n3 | theta | x2.0 | is_pulse | 27.50 | 23.55 (22.90–24.20) | 222 / 869 | 7.08e-91 | run |
| 1274 | COMP_I1_negfeedback | n3 | theta | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1282 | I1_miRNA_FFL | n3 | gx | x0.5 | is_memory | 36.44 | 49.87 (49.10–50.64) | 2250 / 49 | 0 | run |
| 1283 | I1_miRNA_FFL | n3 | gx | x0.5 | is_pulse | 13.91 | 14.58 (14.04–15.13) | 303 / 193 | 8.95e-07 | run |
| 1284 | I1_miRNA_FFL | n3 | gx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1292 | I1_miRNA_FFL | n3 | gx | x2.0 | is_memory | 36.44 | 28.03 (27.34–28.72) | 50 / 1428 | 0 | run |
| 1293 | I1_miRNA_FFL | n3 | gx | x2.0 | is_pulse | 13.91 | 13.80 (13.28–14.34) | 173 / 191 | 0.373 | run |
| 1294 | I1_miRNA_FFL | n3 | gx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1302 | I1_miRNA_FFL | n3 | gp | x0.5 | is_memory | 36.44 | 56.62 (55.85–57.38) | 3344 / 38 | 0 | run |
| 1303 | I1_miRNA_FFL | n3 | gp | x0.5 | is_pulse | 13.91 | 10.05 (9.59–10.52) | 20 / 653 | 5.91e-165 | run |
| 1304 | I1_miRNA_FFL | n3 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1312 | I1_miRNA_FFL | n3 | gp | x2.0 | is_memory | 36.44 | 17.80 (17.21–18.39) | 9 / 3063 | 0 | run |
| 1313 | I1_miRNA_FFL | n3 | gp | x2.0 | is_pulse | 13.91 | 17.90 (17.31–18.49) | 657 / 4 | 1.66e-189 | run |
| 1314 | I1_miRNA_FFL | n3 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1322 | I1_miRNA_FFL | n3 | KS | x0.5 | is_memory | 36.44 | 34.88 (34.15–35.62) | 823 / 1078 | 5.41e-09 | run |
| 1323 | I1_miRNA_FFL | n3 | KS | x0.5 | is_pulse | 13.91 | 9.66 (9.21–10.12) | 392 / 1088 | 8.09e-76 | run |
| 1324 | I1_miRNA_FFL | n3 | KS | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1332 | I1_miRNA_FFL | n3 | KS | x2.0 | is_memory | 36.44 | 35.82 (35.08–36.56) | 943 / 1045 | 0.0235 | run |
| 1333 | I1_miRNA_FFL | n3 | KS | x2.0 | is_pulse | 13.91 | 18.48 (17.88–19.08) | 1255 / 507 | 4.55e-73 | run |
| 1334 | I1_miRNA_FFL | n3 | KS | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1342 | I1_miRNA_FFL | n3 | Kxm | x0.5 | is_memory | 36.44 | 34.60 (33.87–35.33) | 623 / 924 | 1.99e-14 | run |
| 1343 | I1_miRNA_FFL | n3 | Kxm | x0.5 | is_pulse | 13.91 | 12.24 (11.75–12.76) | 511 / 784 | 3.29e-14 | run |
| 1344 | I1_miRNA_FFL | n3 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1352 | I1_miRNA_FFL | n3 | Kxm | x2.0 | is_memory | 36.44 | 37.38 (36.64–38.13) | 947 / 792 | 0.00022 | run |
| 1353 | I1_miRNA_FFL | n3 | Kxm | x2.0 | is_pulse | 13.91 | 13.45 (12.93–13.98) | 806 / 881 | 0.0716 | run |
| 1354 | I1_miRNA_FFL | n3 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1362 | I1_miRNA_FFL | n3 | Kxy | x0.5 | is_memory | 36.44 | 41.53 (40.78–42.29) | 1438 / 603 | 2.64e-78 | run |
| 1363 | I1_miRNA_FFL | n3 | Kxy | x0.5 | is_pulse | 13.91 | 12.64 (12.14–13.16) | 584 / 792 | 2.26e-08 | run |
| 1364 | I1_miRNA_FFL | n3 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1372 | I1_miRNA_FFL | n3 | Kxy | x2.0 | is_memory | 36.44 | 30.51 (29.80–31.22) | 460 / 1432 | 1.01e-115 | run |
| 1373 | I1_miRNA_FFL | n3 | Kxy | x2.0 | is_pulse | 13.91 | 17.91 (17.32–18.50) | 1147 / 492 | 2.61e-60 | run |
| 1374 | I1_miRNA_FFL | n3 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1382 | I1_miRNA_FFL | n3 | Kmy | x0.5 | is_memory | 36.44 | 35.88 (35.15–36.62) | 431 / 522 | 0.00353 | run |
| 1383 | I1_miRNA_FFL | n3 | Kmy | x0.5 | is_pulse | 13.91 | 16.41 (15.85–16.99) | 762 / 352 | 2.68e-35 | run |
| 1384 | I1_miRNA_FFL | n3 | Kmy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1392 | I1_miRNA_FFL | n3 | Kmy | x2.0 | is_memory | 36.44 | 36.04 (35.31–36.78) | 340 / 405 | 0.019 | run |
| 1393 | I1_miRNA_FFL | n3 | Kmy | x2.0 | is_pulse | 13.91 | 11.08 (10.61–11.57) | 269 / 732 | 4.27e-50 | run |
| 1394 | I1_miRNA_FFL | n3 | Kmy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1402 | I1_miRNA_FFL | n3 | Kmp | x0.5 | is_memory | 36.44 | 37.95 (37.20–38.69) | 632 / 385 | 9.07e-15 | run |
| 1403 | I1_miRNA_FFL | n3 | Kmp | x0.5 | is_pulse | 13.91 | 15.66 (15.11–16.23) | 653 / 366 | 1.89e-19 | run |
| 1404 | I1_miRNA_FFL | n3 | Kmp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1412 | I1_miRNA_FFL | n3 | Kmp | x2.0 | is_memory | 36.44 | 34.92 (34.19–35.65) | 246 / 495 | 3.64e-20 | run |
| 1413 | I1_miRNA_FFL | n3 | Kmp | x2.0 | is_pulse | 13.91 | 12.03 (11.54–12.54) | 274 / 582 | 2.68e-26 | run |
| 1414 | I1_miRNA_FFL | n3 | Kmp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1422 | I1_miRNA_FFL | n3 | lam_y | x0.5 | is_memory | 36.44 | 35.55 (34.81–36.29) | 251 / 397 | 1.07e-08 | run |
| 1423 | I1_miRNA_FFL | n3 | lam_y | x0.5 | is_pulse | 13.91 | 10.66 (10.19–11.15) | 177 / 709 | 4.62e-76 | run |
| 1424 | I1_miRNA_FFL | n3 | lam_y | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1432 | I1_miRNA_FFL | n3 | lam_y | x2.0 | is_memory | 36.44 | 37.18 (36.44–37.92) | 474 / 353 | 2.92e-05 | run |
| 1433 | I1_miRNA_FFL | n3 | lam_y | x2.0 | is_pulse | 13.91 | 17.90 (17.31–18.49) | 855 / 202 | 4.55e-96 | run |
| 1434 | I1_miRNA_FFL | n3 | lam_y | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1442 | I1_miRNA_FFL | n3 | theta | x0.5 | is_memory | 36.44 | 37.68 (36.94–38.43) | 566 / 362 | 2.23e-11 | run |
| 1443 | I1_miRNA_FFL | n3 | theta | x0.5 | is_pulse | 13.91 | 18.68 (18.09–19.29) | 969 / 187 | 1.27e-127 | run |
| 1444 | I1_miRNA_FFL | n3 | theta | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1452 | I1_miRNA_FFL | n3 | theta | x2.0 | is_memory | 36.44 | 35.05 (34.32–35.78) | 314 / 542 | 5.97e-15 | run |
| 1453 | I1_miRNA_FFL | n3 | theta | x2.0 | is_pulse | 13.91 | 10.06 (9.61–10.54) | 156 / 786 | 9.92e-102 | run |
| 1454 | I1_miRNA_FFL | n3 | theta | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1462 | COMP_C2_toggle | n4 | gm | x0.5 | is_memory | 59.11 | 66.65 (65.92–67.37) | 1357 / 121 | 7.05e-265 | run |
| 1463 | COMP_C2_toggle | n4 | gm | x0.5 | is_pulse | 0.00 | 0.43 (0.34–0.55) | 71 / 0 | 8.47e-22 | run |
| 1464 | COMP_C2_toggle | n4 | gm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1472 | COMP_C2_toggle | n4 | gm | x2.0 | is_memory | 59.11 | 51.14 (50.37–51.90) | 113 / 1419 | 8.04e-288 | run |
| 1473 | COMP_C2_toggle | n4 | gm | x2.0 | is_pulse | 0.00 | 0.01 (0.00–0.04) | 2 / 0 | 0.5 | run |
| 1474 | COMP_C2_toggle | n4 | gm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1482 | COMP_C2_toggle | n4 | gp | x0.5 | is_memory | 59.11 | 74.50 (73.82–75.17) | 2685 / 163 | 0 | run |
| 1483 | COMP_C2_toggle | n4 | gp | x0.5 | is_pulse | 0.00 | 0.01 (0.00–0.04) | 2 / 0 | 0.5 | run |
| 1484 | COMP_C2_toggle | n4 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1492 | COMP_C2_toggle | n4 | gp | x2.0 | is_memory | 59.11 | 41.41 (40.66–42.17) | 124 / 3023 | 0 | run |
| 1493 | COMP_C2_toggle | n4 | gp | x2.0 | is_pulse | 0.00 | 0.01 (0.00–0.04) | 2 / 0 | 0.5 | run |
| 1494 | COMP_C2_toggle | n4 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1502 | COMP_C2_toggle | n4 | Kxm | x0.5 | is_memory | 59.11 | 66.03 (65.30–66.75) | 1278 / 144 | 2.02e-227 | run |
| 1503 | COMP_C2_toggle | n4 | Kxm | x0.5 | is_pulse | 0.00 | 0.01 (0.00–0.04) | 2 / 0 | 0.5 | run |
| 1504 | COMP_C2_toggle | n4 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1512 | COMP_C2_toggle | n4 | Kxm | x2.0 | is_memory | 59.11 | 52.23 (51.47–53.00) | 116 / 1242 | 1.72e-238 | run |
| 1513 | COMP_C2_toggle | n4 | Kxm | x2.0 | is_pulse | 0.00 | 0.01 (0.00–0.03) | 1 / 0 | 1 | run |
| 1514 | COMP_C2_toggle | n4 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1522 | COMP_C2_toggle | n4 | Kmx | x0.5 | is_memory | 59.11 | 49.84 (49.07–50.61) | 62 / 1580 | 0 | run |
| 1523 | COMP_C2_toggle | n4 | Kmx | x0.5 | is_pulse | 0.00 | 0.01 (0.00–0.03) | 1 / 0 | 1 | run |
| 1524 | COMP_C2_toggle | n4 | Kmx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1532 | COMP_C2_toggle | n4 | Kmx | x2.0 | is_memory | 59.11 | 65.89 (65.16–66.61) | 1179 / 68 | 1.81e-262 | run |
| 1533 | COMP_C2_toggle | n4 | Kmx | x2.0 | is_pulse | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1534 | COMP_C2_toggle | n4 | Kmx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1542 | COMP_C2_toggle | n4 | kappa | x0.5 | is_memory | 59.11 | 49.26 (48.49–50.03) | 128 / 1741 | 0 | run |
| 1543 | COMP_C2_toggle | n4 | kappa | x0.5 | is_pulse | 0.00 | 0.01 (0.00–0.04) | 2 / 0 | 0.5 | run |
| 1544 | COMP_C2_toggle | n4 | kappa | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1552 | COMP_C2_toggle | n4 | kappa | x2.0 | is_memory | 59.11 | 66.36 (65.63–67.08) | 1338 / 150 | 1.5e-238 | run |
| 1553 | COMP_C2_toggle | n4 | kappa | x2.0 | is_pulse | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1554 | COMP_C2_toggle | n4 | kappa | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1562 | COMP_I1_negfeedback | n4 | gm | x0.5 | is_memory | 60.27 | 66.38 (65.65–67.11) | 1361 / 359 | 4.85e-137 | run |
| 1563 | COMP_I1_negfeedback | n4 | gm | x0.5 | is_pulse | 19.84 | 18.49 (17.90–19.10) | 630 / 851 | 1.02e-08 | run |
| 1564 | COMP_I1_negfeedback | n4 | gm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1572 | COMP_I1_negfeedback | n4 | gm | x2.0 | is_memory | 60.27 | 54.49 (53.73–55.26) | 401 / 1347 | 2.84e-119 | run |
| 1573 | COMP_I1_negfeedback | n4 | gm | x2.0 | is_pulse | 19.84 | 18.09 (17.50–18.69) | 681 / 968 | 1.66e-12 | run |
| 1574 | COMP_I1_negfeedback | n4 | gm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1582 | COMP_I1_negfeedback | n4 | gp | x0.5 | is_memory | 60.27 | 74.32 (73.64–74.98) | 2497 / 195 | 0 | run |
| 1583 | COMP_I1_negfeedback | n4 | gp | x0.5 | is_pulse | 19.84 | 15.48 (14.93–16.05) | 45 / 759 | 2.59e-168 | run |
| 1584 | COMP_I1_negfeedback | n4 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1592 | COMP_I1_negfeedback | n4 | gp | x2.0 | is_memory | 60.27 | 43.57 (42.81–44.34) | 105 / 2840 | 0 | run |
| 1593 | COMP_I1_negfeedback | n4 | gp | x2.0 | is_pulse | 19.84 | 23.52 (22.87–24.17) | 632 / 30 | 8.97e-148 | run |
| 1594 | COMP_I1_negfeedback | n4 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1602 | COMP_I1_negfeedback | n4 | Kxm | x0.5 | is_memory | 60.27 | 54.82 (54.05–55.58) | 467 / 1360 | 6.98e-101 | run |
| 1603 | COMP_I1_negfeedback | n4 | Kxm | x0.5 | is_pulse | 19.84 | 25.81 (25.14–26.48) | 1360 / 383 | 5.19e-128 | run |
| 1604 | COMP_I1_negfeedback | n4 | Kxm | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1612 | COMP_I1_negfeedback | n4 | Kxm | x2.0 | is_memory | 60.27 | 66.06 (65.33–66.78) | 1419 / 470 | 2.57e-110 | run |
| 1613 | COMP_I1_negfeedback | n4 | Kxm | x2.0 | is_pulse | 19.84 | 13.89 (13.36–14.42) | 291 / 1267 | 4.29e-145 | run |
| 1614 | COMP_I1_negfeedback | n4 | Kxm | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1622 | COMP_I1_negfeedback | n4 | Kxy | x0.5 | is_memory | 60.27 | 68.99 (68.27–69.70) | 1882 / 453 | 5.47e-206 | run |
| 1623 | COMP_I1_negfeedback | n4 | Kxy | x0.5 | is_pulse | 19.84 | 12.13 (11.63–12.64) | 212 / 1476 | 5.58e-233 | run |
| 1624 | COMP_I1_negfeedback | n4 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1632 | COMP_I1_negfeedback | n4 | Kxy | x2.0 | is_memory | 60.27 | 50.04 (49.27–50.81) | 400 / 2076 | 8.32e-272 | run |
| 1633 | COMP_I1_negfeedback | n4 | Kxy | x2.0 | is_pulse | 19.84 | 29.98 (29.28–30.69) | 1931 / 270 | 1.01e-308 | run |
| 1634 | COMP_I1_negfeedback | n4 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1642 | COMP_I1_negfeedback | n4 | nxy | x0.5 | is_memory | 60.27 | 66.49 (65.76–67.21) | 1724 / 704 | 6.04e-98 | run |
| 1643 | COMP_I1_negfeedback | n4 | nxy | x0.5 | is_pulse | 19.84 | 10.68 (10.21–11.16) | 318 / 1820 | 4.18e-255 | run |
| 1644 | COMP_I1_negfeedback | n4 | nxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1652 | COMP_I1_negfeedback | n4 | nxy | x2.0 | is_memory | 60.27 | 55.76 (55.00–56.52) | 540 / 1278 | 7.58e-69 | run |
| 1653 | COMP_I1_negfeedback | n4 | nxy | x2.0 | is_pulse | 19.84 | 26.10 (25.43–26.78) | 1720 / 694 | 1.02e-99 | run |
| 1654 | COMP_I1_negfeedback | n4 | nxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1662 | COMP_I1_negfeedback | n4 | Kmx | x0.5 | is_memory | 60.27 | 56.45 (55.68–57.21) | 400 / 1026 | 1.35e-63 | run |
| 1663 | COMP_I1_negfeedback | n4 | Kmx | x0.5 | is_pulse | 19.84 | 23.38 (22.74–24.04) | 1044 / 464 | 1.55e-51 | run |
| 1664 | COMP_I1_negfeedback | n4 | Kmx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1672 | COMP_I1_negfeedback | n4 | Kmx | x2.0 | is_memory | 60.27 | 63.27 (62.53–64.01) | 838 / 346 | 1.47e-47 | run |
| 1673 | COMP_I1_negfeedback | n4 | Kmx | x2.0 | is_pulse | 19.84 | 16.06 (15.51–16.64) | 353 / 972 | 3.85e-67 | run |
| 1674 | COMP_I1_negfeedback | n4 | Kmx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1682 | COMP_I1_negfeedback | n4 | kappa | x0.5 | is_memory | 60.27 | 51.05 (50.28–51.82) | 98 / 1608 | 0 | run |
| 1683 | COMP_I1_negfeedback | n4 | kappa | x0.5 | is_pulse | 19.84 | 20.51 (19.90–21.14) | 138 / 28 | 1.14e-18 | run |
| 1684 | COMP_I1_negfeedback | n4 | kappa | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1692 | COMP_I1_negfeedback | n4 | kappa | x2.0 | is_memory | 60.27 | 67.22 (66.50–67.94) | 1284 / 144 | 5.98e-229 | run |
| 1693 | COMP_I1_negfeedback | n4 | kappa | x2.0 | is_pulse | 19.84 | 19.16 (18.56–19.77) | 41 / 153 | 2.14e-16 | run |
| 1694 | COMP_I1_negfeedback | n4 | kappa | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1702 | I1_miRNA_FFL | n4 | gx | x0.5 | is_memory | 68.48 | 74.83 (74.16–75.49) | 1125 / 84 | 3.24e-233 | run |
| 1703 | I1_miRNA_FFL | n4 | gx | x0.5 | is_pulse | 6.77 | 7.15 (6.76–7.55) | 156 / 95 | 0.000142 | run |
| 1704 | I1_miRNA_FFL | n4 | gx | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1712 | I1_miRNA_FFL | n4 | gx | x2.0 | is_memory | 68.48 | 64.09 (63.35–64.82) | 84 / 803 | 4.76e-148 | run |
| 1713 | I1_miRNA_FFL | n4 | gx | x2.0 | is_pulse | 6.77 | 6.89 (6.51–7.29) | 113 / 94 | 0.211 | run |
| 1714 | I1_miRNA_FFL | n4 | gx | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1722 | I1_miRNA_FFL | n4 | gp | x0.5 | is_memory | 68.48 | 81.43 (80.83–82.03) | 2303 / 180 | 0 | run |
| 1723 | I1_miRNA_FFL | n4 | gp | x0.5 | is_pulse | 6.77 | 5.02 (4.69–5.36) | 49 / 337 | 5.59e-54 | run |
| 1724 | I1_miRNA_FFL | n4 | gp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1732 | I1_miRNA_FFL | n4 | gp | x2.0 | is_memory | 68.48 | 50.84 (50.07–51.61) | 111 / 3000 | 0 | run |
| 1733 | I1_miRNA_FFL | n4 | gp | x2.0 | is_pulse | 6.77 | 9.01 (8.58–9.46) | 400 / 33 | 3.27e-81 | run |
| 1734 | I1_miRNA_FFL | n4 | gp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1742 | I1_miRNA_FFL | n4 | KS | x0.5 | is_memory | 68.48 | 68.35 (67.63–69.07) | 780 / 800 | 0.633 | run |
| 1743 | I1_miRNA_FFL | n4 | KS | x0.5 | is_pulse | 6.77 | 4.57 (4.25–4.90) | 204 / 566 | 3.57e-40 | run |
| 1744 | I1_miRNA_FFL | n4 | KS | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1752 | I1_miRNA_FFL | n4 | KS | x2.0 | is_memory | 68.48 | 66.51 (65.78–67.23) | 649 / 971 | 1.23e-15 | run |
| 1753 | I1_miRNA_FFL | n4 | KS | x2.0 | is_pulse | 6.77 | 9.81 (9.36–10.28) | 750 / 252 | 5.36e-58 | run |
| 1754 | I1_miRNA_FFL | n4 | KS | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1762 | I1_miRNA_FFL | n4 | Kxy | x0.5 | is_memory | 68.48 | 75.00 (74.33–75.66) | 1532 / 463 | 4.3e-133 | run |
| 1763 | I1_miRNA_FFL | n4 | Kxy | x0.5 | is_pulse | 6.77 | 6.05 (5.69–6.43) | 300 / 418 | 1.21e-05 | run |
| 1764 | I1_miRNA_FFL | n4 | Kxy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1772 | I1_miRNA_FFL | n4 | Kxy | x2.0 | is_memory | 68.48 | 58.50 (57.74–59.26) | 363 / 1997 | 1.93e-272 | run |
| 1773 | I1_miRNA_FFL | n4 | Kxy | x2.0 | is_pulse | 6.77 | 10.93 (10.45–11.41) | 909 / 229 | 3e-96 | run |
| 1774 | I1_miRNA_FFL | n4 | Kxy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1782 | I1_miRNA_FFL | n4 | Kmy | x0.5 | is_memory | 68.48 | 65.46 (64.73–66.19) | 247 / 741 | 6.51e-58 | run |
| 1783 | I1_miRNA_FFL | n4 | Kmy | x0.5 | is_pulse | 6.77 | 8.57 (8.15–9.01) | 480 / 186 | 7.44e-31 | run |
| 1784 | I1_miRNA_FFL | n4 | Kmy | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1792 | I1_miRNA_FFL | n4 | Kmy | x2.0 | is_memory | 68.48 | 70.88 (70.18–71.58) | 606 / 212 | 1.04e-44 | run |
| 1793 | I1_miRNA_FFL | n4 | Kmy | x2.0 | is_pulse | 6.77 | 5.02 (4.69–5.36) | 144 / 432 | 2.17e-34 | run |
| 1794 | I1_miRNA_FFL | n4 | Kmy | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1802 | I1_miRNA_FFL | n4 | Kmp | x0.5 | is_memory | 68.48 | 66.42 (65.70–67.15) | 346 / 682 | 5.43e-26 | run |
| 1803 | I1_miRNA_FFL | n4 | Kmp | x0.5 | is_pulse | 6.77 | 8.58 (8.15–9.01) | 478 / 183 | 2.67e-31 | run |
| 1804 | I1_miRNA_FFL | n4 | Kmp | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1812 | I1_miRNA_FFL | n4 | Kmp | x2.0 | is_memory | 68.48 | 69.94 (69.23–70.64) | 500 / 260 | 2.28e-18 | run |
| 1813 | I1_miRNA_FFL | n4 | Kmp | x2.0 | is_pulse | 6.77 | 5.41 (5.07–5.77) | 132 / 356 | 8.38e-25 | run |
| 1814 | I1_miRNA_FFL | n4 | Kmp | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1822 | I1_miRNA_FFL | n4 | lam_y | x0.5 | is_memory | 68.48 | 70.27 (69.56–70.97) | 537 / 243 | 2.41e-26 | run |
| 1823 | I1_miRNA_FFL | n4 | lam_y | x0.5 | is_pulse | 6.77 | 4.65 (4.33–4.98) | 103 / 451 | 6.76e-53 | run |
| 1824 | I1_miRNA_FFL | n4 | lam_y | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1832 | I1_miRNA_FFL | n4 | lam_y | x2.0 | is_memory | 68.48 | 65.80 (65.07–66.53) | 269 / 707 | 5.23e-46 | run |
| 1833 | I1_miRNA_FFL | n4 | lam_y | x2.0 | is_pulse | 6.77 | 10.04 (9.58–10.51) | 664 / 129 | 1.58e-87 | run |
| 1834 | I1_miRNA_FFL | n4 | lam_y | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1842 | I1_miRNA_FFL | n4 | theta | x0.5 | is_memory | 68.48 | 66.60 (65.87–67.32) | 344 / 652 | 1.05e-22 | run |
| 1843 | I1_miRNA_FFL | n4 | theta | x0.5 | is_pulse | 6.77 | 10.04 (9.58–10.51) | 651 / 116 | 4.16e-91 | run |
| 1844 | I1_miRNA_FFL | n4 | theta | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1852 | I1_miRNA_FFL | n4 | theta | x2.0 | is_memory | 68.48 | 69.90 (69.19–70.60) | 525 / 292 | 3.03e-16 | run |
| 1853 | I1_miRNA_FFL | n4 | theta | x2.0 | is_pulse | 6.77 | 4.74 (4.42–5.08) | 110 / 443 | 2.51e-48 | run |
| 1854 | I1_miRNA_FFL | n4 | theta | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1862 | I1_miRNA_FFL | n4 | kappa | x0.5 | is_memory | 68.48 | 58.12 (57.36–58.88) | 111 / 1807 | 0 | run |
| 1863 | I1_miRNA_FFL | n4 | kappa | x0.5 | is_pulse | 6.77 | 7.45 (7.05–7.87) | 144 / 33 | 9.79e-18 | run |
| 1864 | I1_miRNA_FFL | n4 | kappa | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1872 | I1_miRNA_FFL | n4 | kappa | x2.0 | is_memory | 68.48 | 75.97 (75.31–76.62) | 1400 / 172 | 2.4e-239 | run |
| 1873 | I1_miRNA_FFL | n4 | kappa | x2.0 | is_pulse | 6.77 | 6.26 (5.90–6.64) | 49 / 133 | 3.73e-10 | run |
| 1874 | I1_miRNA_FFL | n4 | kappa | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1882 | I1_miRNA_FFL | n4 | Kg12 | x0.5 | is_memory | 68.48 | 71.96 (71.27–72.65) | 631 / 60 | 4.32e-121 | run |
| 1883 | I1_miRNA_FFL | n4 | Kg12 | x0.5 | is_pulse | 6.77 | 6.32 (5.96–6.71) | 22 / 96 | 3.26e-12 | run |
| 1884 | I1_miRNA_FFL | n4 | Kg12 | x0.5 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |
| 1892 | I1_miRNA_FFL | n4 | Kg12 | x2.0 | is_memory | 68.48 | 63.40 (62.65–64.14) | 78 / 910 | 1.27e-180 | run |
| 1893 | I1_miRNA_FFL | n4 | Kg12 | x2.0 | is_pulse | 6.77 | 7.37 (6.97–7.78) | 115 / 18 | 1.72e-18 | run |
| 1894 | I1_miRNA_FFL | n4 | Kg12 | x2.0 | is_sustained_osc | 0.00 | 0.00 (0.00–0.02) | 0 / 0 | 1 | run |

### Range over the excursions, per module and behaviour (`excursion_range_by_module.csv`, descriptive)

Columns baseline_pct, min_pct (min_at), max_pct (max_at), median_abs_change_pp, max_abs_change_pp, appears (0 at
baseline, above 0 in an excursion), vanishes (above 0 at baseline, 0 in an excursion).

| line | family | module | behaviour | excursions | baseline % | lowest % (excursion) | highest % (excursion) | median abs change (pp) | max abs change (pp) | appears | vanishes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | COMP_C2_toggle | n6 | is_memory | 20 | 58.55 | 40.63 (gp x2) | 73.51 (gp x0.5) | 4.72 | 17.92 | False | False |
| 3 | COMP_C2_toggle | n6 | is_pulse | 20 | 25.12 | 16.23 (Kxm x0.5) | 33.59 (Kxm x2) | 3.48 | 8.89 | False | False |
| 4 | COMP_C2_toggle | n6 | is_sustained_osc | 20 | 0.66 | 0.08 (nmx x0.5) | 1.81 (nmx x2) | 0.23 | 1.15 | False | False |
| 5 | COMP_C2_toggle | n6 | is_damped_osc | 20 | 5.26 | 3.05 (Kmx x2) | 7.99 (gx x0.5) | 1.40 | 2.73 | False | False |
| 6 | COMP_C2_toggle | n6 | is_bistable | 20 | 12.49 | 9.58 (Kx2m x2) | 15.12 (Kx2m x0.5) | 1.15 | 2.91 | False | False |
| 7 | COMP_C2_toggle | n6 | bistability_unresolved | 20 | 8.66 | 4.71 (gp x2) | 21.23 (gp x0.5) | 1.55 | 12.57 | False | False |
| 8 | COMP_C2_toggle | n6 | is_ultrasensitive | 20 | 66.15 | 55.26 (gp x0.5) | 74.50 (gp x2) | 1.71 | 10.89 | False | False |
| 9 | COMP_C2_toggle | n6 | is_very_ultrasensitive | 20 | 17.42 | 13.54 (gp x0.5) | 21.59 (gp x2) | 1.35 | 4.17 | False | False |
| 10 | COMP_C2_toggle | n6 | is_noise_rejecting | 20 | 80.22 | 73.27 (Kmx x0.5) | 86.15 (Kmx x2) | 1.94 | 6.95 | False | False |
| 11 | COMP_C2_toggle | n6 | nonmonotone_dose | 20 | 16.22 | 11.83 (Kxm x0.5) | 20.22 (Kxm x2) | 1.19 | 4.39 | False | False |
| 12 | COMP_I1_negfeedback | n6 | is_memory | 20 | 50.52 | 34.41 (gp x2) | 66.31 (gp x0.5) | 4.97 | 16.12 | False | False |
| 13 | COMP_I1_negfeedback | n6 | is_pulse | 20 | 30.76 | 22.99 (Kxm x2) | 39.39 (Kxy x2) | 3.40 | 8.62 | False | False |
| 14 | COMP_I1_negfeedback | n6 | is_sustained_osc | 20 | 0.23 | 0.01 (nx12 x0.5) | 1.12 (nx12 x2) | 0.06 | 0.89 | False | False |
| 15 | COMP_I1_negfeedback | n6 | is_damped_osc | 20 | 3.88 | 1.40 (nx12 x0.5) | 6.85 (gx x0.5) | 0.91 | 2.98 | False | False |
| 16 | COMP_I1_negfeedback | n6 | is_bistable | 20 | 2.60 | 1.50 (kappa x0.5) | 3.31 (kappa x2) | 0.26 | 1.10 | False | False |
| 17 | COMP_I1_negfeedback | n6 | bistability_unresolved | 20 | 5.03 | 1.21 (gp x2) | 16.53 (gp x0.5) | 0.59 | 11.50 | False | False |
| 18 | COMP_I1_negfeedback | n6 | is_ultrasensitive | 20 | 76.46 | 67.41 (gp x0.5) | 83.34 (gp x2) | 0.98 | 9.06 | False | False |
| 19 | COMP_I1_negfeedback | n6 | is_very_ultrasensitive | 20 | 32.20 | 26.39 (gp x0.5) | 37.20 (gp x2) | 2.91 | 5.82 | False | False |
| 20 | COMP_I1_negfeedback | n6 | is_noise_rejecting | 20 | 73.10 | 69.05 (Kxm x0.5) | 78.58 (Kxm x2) | 1.86 | 5.48 | False | False |
| 21 | COMP_I1_negfeedback | n6 | nonmonotone_dose | 20 | 31.80 | 26.36 (Kmx x0.5) | 34.90 (Kxy x0.5) | 2.18 | 5.44 | False | False |
| 22 | I1_miRNA_FFL | n6 | is_memory | 20 | 59.34 | 42.30 (gp x2) | 74.41 (gp x0.5) | 3.65 | 17.04 | False | False |
| 23 | I1_miRNA_FFL | n6 | is_pulse | 20 | 16.33 | 11.63 (KS x0.5) | 21.26 (Kxy x2) | 3.12 | 4.94 | False | False |
| 24 | I1_miRNA_FFL | n6 | is_sustained_osc | 20 | 0.00 | 0.00 (gp x0.5) | 0.00 (gp x0.5) | 0.00 | 0.00 | False | False |
| 25 | I1_miRNA_FFL | n6 | is_damped_osc | 20 | 0.00 | 0.00 (gp x0.5) | 0.00 (gp x0.5) | 0.00 | 0.00 | False | False |
| 26 | I1_miRNA_FFL | n6 | is_bistable | 20 | 2.25 | 0.87 (kappa x0.5) | 3.41 (kappa x2) | 0.20 | 1.38 | False | False |
| 27 | I1_miRNA_FFL | n6 | bistability_unresolved | 20 | 5.44 | 0.60 (gp x2) | 17.66 (gp x0.5) | 0.54 | 12.22 | False | False |
| 28 | I1_miRNA_FFL | n6 | is_ultrasensitive | 20 | 77.36 | 66.69 (gp x0.5) | 85.23 (gp x2) | 1.46 | 10.68 | False | False |
| 29 | I1_miRNA_FFL | n6 | is_very_ultrasensitive | 20 | 32.27 | 26.17 (gp x0.5) | 38.82 (Kxy x0.5) | 2.90 | 6.55 | False | False |
| 30 | I1_miRNA_FFL | n6 | is_noise_rejecting | 20 | 82.62 | 78.38 (Kmp x0.5) | 86.02 (Kxy x2) | 1.80 | 4.24 | False | False |
| 31 | I1_miRNA_FFL | n6 | nonmonotone_dose | 20 | 33.01 | 27.07 (Kmp x2) | 40.34 (Kmp x0.5) | 3.28 | 7.33 | False | False |
| 32 | COMP_C2_toggle | n5 | is_memory | 10 | 48.55 | 31.85 (gp x2) | 65.78 (gp x0.5) | 8.23 | 17.23 | False | False |
| 33 | COMP_C2_toggle | n5 | is_pulse | 10 | 0.00 | 0.00 (gp x2) | 0.02 (Kmx x0.5) | 0.01 | 0.02 | True | False |
| 34 | COMP_C2_toggle | n5 | is_sustained_osc | 10 | 0.00 | 0.00 (gp x0.5) | 0.00 (gp x0.5) | 0.00 | 0.00 | False | False |
| 35 | COMP_C2_toggle | n5 | is_damped_osc | 10 | 0.00 | 0.00 (gp x0.5) | 0.00 (gp x0.5) | 0.00 | 0.00 | False | False |
| 36 | COMP_C2_toggle | n5 | is_bistable | 10 | 26.07 | 19.64 (Kxm x2) | 30.63 (Kxm x0.5) | 2.37 | 6.44 | False | False |
| 37 | COMP_C2_toggle | n5 | bistability_unresolved | 10 | 3.13 | 0.20 (gp x2) | 15.06 (gp x0.5) | 0.78 | 11.93 | False | False |
| 38 | COMP_C2_toggle | n5 | is_ultrasensitive | 10 | 42.61 | 32.31 (Kmx x0.5) | 52.31 (Kmx x2) | 4.46 | 10.30 | False | False |
| 39 | COMP_C2_toggle | n5 | is_very_ultrasensitive | 10 | 7.45 | 4.50 (gp x0.5) | 10.72 (gp x2) | 2.34 | 3.27 | False | False |
| 40 | COMP_C2_toggle | n5 | is_noise_rejecting | 10 | 72.65 | 61.59 (Kmx x0.5) | 81.77 (Kmx x2) | 1.87 | 11.06 | False | False |
| 41 | COMP_C2_toggle | n5 | nonmonotone_dose | 10 | 6.30 | 3.29 (gp x2) | 8.87 (gp x0.5) | 0.20 | 3.02 | False | False |
| 42 | COMP_I1_negfeedback | n5 | is_memory | 14 | 53.96 | 39.03 (gp x2) | 68.71 (gp x0.5) | 6.45 | 14.92 | False | False |
| 43 | COMP_I1_negfeedback | n5 | is_pulse | 14 | 25.66 | 13.48 (nxy x0.5) | 34.97 (Kxy x2) | 5.36 | 12.18 | False | False |
| 44 | COMP_I1_negfeedback | n5 | is_sustained_osc | 14 | 0.00 | 0.00 (gp x0.5) | 0.00 (gp x0.5) | 0.00 | 0.00 | False | False |
| 45 | COMP_I1_negfeedback | n5 | is_damped_osc | 14 | 0.25 | 0.04 (gp x0.5) | 0.77 (gp x2) | 0.05 | 0.52 | False | False |
| 46 | COMP_I1_negfeedback | n5 | is_bistable | 14 | 3.07 | 1.97 (kappa x0.5) | 4.85 (nxy x2) | 0.55 | 1.78 | False | False |
| 47 | COMP_I1_negfeedback | n5 | bistability_unresolved | 14 | 4.77 | 0.52 (gp x2) | 16.76 (gp x0.5) | 0.98 | 11.99 | False | False |
| 48 | COMP_I1_negfeedback | n5 | is_ultrasensitive | 14 | 72.74 | 63.73 (gp x0.5) | 79.25 (gp x2) | 2.27 | 9.01 | False | False |
| 49 | COMP_I1_negfeedback | n5 | is_very_ultrasensitive | 14 | 28.39 | 22.52 (gp x0.5) | 33.27 (gp x2) | 2.22 | 5.87 | False | False |
| 50 | COMP_I1_negfeedback | n5 | is_noise_rejecting | 14 | 63.62 | 58.28 (Kxm x0.5) | 70.07 (Kxm x2) | 2.56 | 6.46 | False | False |
| 51 | COMP_I1_negfeedback | n5 | nonmonotone_dose | 14 | 38.78 | 31.66 (nxy x2) | 49.04 (nxy x0.5) | 4.01 | 10.26 | False | False |
| 52 | I1_miRNA_FFL | n5 | is_memory | 18 | 62.32 | 45.82 (gp x2) | 76.30 (gp x0.5) | 4.38 | 16.50 | False | False |
| 53 | I1_miRNA_FFL | n5 | is_pulse | 18 | 11.25 | 7.39 (KS x0.5) | 15.78 (KS x2) | 2.25 | 4.53 | False | False |
| 54 | I1_miRNA_FFL | n5 | is_sustained_osc | 18 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 55 | I1_miRNA_FFL | n5 | is_damped_osc | 18 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 56 | I1_miRNA_FFL | n5 | is_bistable | 18 | 2.40 | 0.88 (kappa x0.5) | 3.60 (kappa x2) | 0.25 | 1.52 | False | False |
| 57 | I1_miRNA_FFL | n5 | bistability_unresolved | 18 | 5.53 | 0.62 (gp x2) | 18.49 (gp x0.5) | 0.26 | 12.96 | False | False |
| 58 | I1_miRNA_FFL | n5 | is_ultrasensitive | 18 | 75.58 | 65.11 (gp x0.5) | 83.22 (gp x2) | 1.66 | 10.47 | False | False |
| 59 | I1_miRNA_FFL | n5 | is_very_ultrasensitive | 18 | 31.23 | 24.40 (Kxy x2) | 38.06 (Kxy x0.5) | 2.55 | 6.83 | False | False |
| 60 | I1_miRNA_FFL | n5 | is_noise_rejecting | 18 | 76.72 | 71.03 (Kxy x0.5) | 81.22 (Kxy x2) | 1.58 | 5.69 | False | False |
| 61 | I1_miRNA_FFL | n5 | nonmonotone_dose | 18 | 40.23 | 31.73 (Kxm x0.5) | 49.13 (Kxy x0.5) | 3.58 | 8.90 | False | False |
| 62 | COMP_C2_toggle | n3 | is_memory | 10 | 24.79 | 6.79 (gp x2) | 47.51 (gp x0.5) | 3.34 | 22.72 | False | False |
| 63 | COMP_C2_toggle | n3 | is_pulse | 10 | 0.02 | 0.01 (gm x2) | 0.57 (gm x0.5) | 0.01 | 0.55 | False | False |
| 64 | COMP_C2_toggle | n3 | is_sustained_osc | 10 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 65 | COMP_C2_toggle | n3 | is_damped_osc | 10 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 66 | COMP_C2_toggle | n3 | is_bistable | 10 | 25.21 | 19.01 (Kxm x2) | 30.01 (Kxm x0.5) | 1.61 | 6.20 | False | False |
| 67 | COMP_C2_toggle | n3 | bistability_unresolved | 10 | 0.05 | 0.01 (gm x2) | 0.19 (gm x0.5) | 0.01 | 0.14 | False | False |
| 68 | COMP_C2_toggle | n3 | is_ultrasensitive | 10 | 60.91 | 55.16 (Kxy x2) | 66.43 (Kxy x0.5) | 1.45 | 5.75 | False | False |
| 69 | COMP_C2_toggle | n3 | is_very_ultrasensitive | 10 | 19.57 | 14.52 (Kxy x2) | 24.62 (Kxy x0.5) | 3.64 | 5.05 | False | False |
| 70 | COMP_C2_toggle | n3 | is_noise_rejecting | 10 | 85.22 | 81.45 (gp x0.5) | 87.72 (Kxy x0.5) | 0.80 | 3.77 | False | False |
| 71 | COMP_C2_toggle | n3 | nonmonotone_dose | 10 | 5.54 | 2.81 (gp x2) | 7.81 (gp x0.5) | 0.56 | 2.72 | False | False |
| 72 | COMP_I1_negfeedback | n3 | is_memory | 16 | 29.50 | 13.04 (gp x2) | 49.58 (gp x0.5) | 2.47 | 20.09 | False | False |
| 73 | COMP_I1_negfeedback | n3 | is_pulse | 16 | 27.50 | 19.08 (Kxy x0.5) | 37.96 (Kxy x2) | 3.75 | 10.46 | False | False |
| 74 | COMP_I1_negfeedback | n3 | is_sustained_osc | 16 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 75 | COMP_I1_negfeedback | n3 | is_damped_osc | 16 | 2.55 | 0.35 (gp x0.5) | 6.51 (gx x0.5) | 0.27 | 3.96 | False | False |
| 76 | COMP_I1_negfeedback | n3 | is_bistable | 16 | 1.50 | 0.95 (Kxm x2) | 2.60 (Kxm x0.5) | 0.05 | 1.10 | False | False |
| 77 | COMP_I1_negfeedback | n3 | bistability_unresolved | 16 | 0.01 | 0.00 (gx x2) | 0.02 (Kmx x0.5) | 0.00 | 0.01 | False | True |
| 78 | COMP_I1_negfeedback | n3 | is_ultrasensitive | 16 | 82.95 | 78.70 (Kxm x0.5) | 86.24 (Kxm x2) | 1.15 | 4.25 | False | False |
| 79 | COMP_I1_negfeedback | n3 | is_very_ultrasensitive | 16 | 33.72 | 27.64 (Kxy x2) | 40.58 (Kxy x0.5) | 1.39 | 6.85 | False | False |
| 80 | COMP_I1_negfeedback | n3 | is_noise_rejecting | 16 | 69.32 | 60.03 (gp x2) | 75.68 (gp x0.5) | 3.19 | 9.29 | False | False |
| 81 | COMP_I1_negfeedback | n3 | nonmonotone_dose | 16 | 42.40 | 35.01 (Kxy x2) | 52.70 (Kxy x0.5) | 2.06 | 10.30 | False | False |
| 82 | I1_miRNA_FFL | n3 | is_memory | 18 | 36.44 | 17.80 (gp x2) | 56.62 (gp x0.5) | 1.51 | 20.18 | False | False |
| 83 | I1_miRNA_FFL | n3 | is_pulse | 18 | 13.91 | 9.66 (KS x0.5) | 18.68 (theta x0.5) | 3.04 | 4.77 | False | False |
| 84 | I1_miRNA_FFL | n3 | is_sustained_osc | 18 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 85 | I1_miRNA_FFL | n3 | is_damped_osc | 18 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 86 | I1_miRNA_FFL | n3 | is_bistable | 18 | 0.48 | 0.10 (lam_y x0.5) | 1.08 (lam_y x2) | 0.03 | 0.60 | False | False |
| 87 | I1_miRNA_FFL | n3 | bistability_unresolved | 18 | 0.00 | 0.00 (gx x0.5) | 0.01 (KS x0.5) | 0.00 | 0.01 | True | False |
| 88 | I1_miRNA_FFL | n3 | is_ultrasensitive | 18 | 86.89 | 83.03 (Kxy x2) | 89.70 (Kxy x0.5) | 0.44 | 3.86 | False | False |
| 89 | I1_miRNA_FFL | n3 | is_very_ultrasensitive | 18 | 35.94 | 27.20 (Kxy x2) | 45.98 (Kxy x0.5) | 1.49 | 10.05 | False | False |
| 90 | I1_miRNA_FFL | n3 | is_noise_rejecting | 18 | 79.39 | 72.34 (gp x2) | 84.47 (KS x0.5) | 2.96 | 7.06 | False | False |
| 91 | I1_miRNA_FFL | n3 | nonmonotone_dose | 18 | 43.32 | 33.40 (Kxm x0.5) | 55.40 (Kxy x0.5) | 3.16 | 12.08 | False | False |
| 92 | COMP_C2_toggle | n4 | is_memory | 10 | 59.11 | 41.41 (gp x2) | 74.50 (gp x0.5) | 7.76 | 17.69 | False | False |
| 93 | COMP_C2_toggle | n4 | is_pulse | 10 | 0.00 | 0.00 (Kmx x2) | 0.43 (gm x0.5) | 0.01 | 0.43 | True | False |
| 94 | COMP_C2_toggle | n4 | is_sustained_osc | 10 | 0.00 | 0.00 (gm x0.5) | 0.00 (gm x0.5) | 0.00 | 0.00 | False | False |
| 95 | COMP_C2_toggle | n4 | is_damped_osc | 10 | 0.00 | 0.00 (gm x0.5) | 0.00 (gm x0.5) | 0.00 | 0.00 | False | False |
| 96 | COMP_C2_toggle | n4 | is_bistable | 10 | 22.83 | 15.89 (Kmx x2) | 28.84 (Kmx x0.5) | 2.22 | 6.95 | False | False |
| 97 | COMP_C2_toggle | n4 | bistability_unresolved | 10 | 3.89 | 0.19 (gp x2) | 17.13 (gp x0.5) | 0.80 | 13.24 | False | False |
| 98 | COMP_C2_toggle | n4 | is_ultrasensitive | 10 | 51.64 | 41.55 (Kmx x0.5) | 59.68 (Kmx x2) | 5.84 | 10.08 | False | False |
| 99 | COMP_C2_toggle | n4 | is_very_ultrasensitive | 10 | 9.93 | 5.98 (gp x0.5) | 14.42 (gp x2) | 2.21 | 4.49 | False | False |
| 100 | COMP_C2_toggle | n4 | is_noise_rejecting | 10 | 80.46 | 70.82 (Kmx x0.5) | 87.39 (Kmx x2) | 2.15 | 9.64 | False | False |
| 101 | COMP_C2_toggle | n4 | nonmonotone_dose | 10 | 5.40 | 2.77 (gp x2) | 7.59 (gp x0.5) | 1.19 | 2.63 | False | False |
| 102 | COMP_I1_negfeedback | n4 | is_memory | 14 | 60.27 | 43.57 (gp x2) | 74.32 (gp x0.5) | 6.17 | 16.69 | False | False |
| 103 | COMP_I1_negfeedback | n4 | is_pulse | 14 | 19.84 | 10.68 (nxy x0.5) | 29.98 (Kxy x2) | 4.07 | 10.14 | False | False |
| 104 | COMP_I1_negfeedback | n4 | is_sustained_osc | 14 | 0.00 | 0.00 (gm x0.5) | 0.00 (gm x0.5) | 0.00 | 0.00 | False | False |
| 105 | COMP_I1_negfeedback | n4 | is_damped_osc | 14 | 1.00 | 0.11 (gp x0.5) | 2.92 (gp x2) | 0.19 | 1.92 | False | False |
| 106 | COMP_I1_negfeedback | n4 | is_bistable | 14 | 3.45 | 2.42 (kappa x0.5) | 5.55 (nxy x2) | 0.71 | 2.11 | False | False |
| 107 | COMP_I1_negfeedback | n4 | bistability_unresolved | 14 | 5.07 | 0.44 (gp x2) | 18.27 (gp x0.5) | 0.68 | 13.20 | False | False |
| 108 | COMP_I1_negfeedback | n4 | is_ultrasensitive | 14 | 71.10 | 59.19 (gp x0.5) | 79.99 (gp x2) | 1.80 | 11.91 | False | False |
| 109 | COMP_I1_negfeedback | n4 | is_very_ultrasensitive | 14 | 23.91 | 17.77 (gp x0.5) | 29.46 (gp x2) | 3.48 | 6.14 | False | False |
| 110 | COMP_I1_negfeedback | n4 | is_noise_rejecting | 14 | 72.19 | 67.92 (Kxm x0.5) | 77.18 (gm x0.5) | 3.31 | 4.99 | False | False |
| 111 | COMP_I1_negfeedback | n4 | nonmonotone_dose | 14 | 28.55 | 20.51 (gm x0.5) | 36.21 (nxy x0.5) | 3.36 | 8.04 | False | False |
| 112 | I1_miRNA_FFL | n4 | is_memory | 20 | 68.48 | 50.84 (gp x2) | 81.43 (gp x0.5) | 3.25 | 17.63 | False | False |
| 113 | I1_miRNA_FFL | n4 | is_pulse | 20 | 6.77 | 4.57 (KS x0.5) | 10.93 (Kxy x2) | 1.78 | 4.15 | False | False |
| 114 | I1_miRNA_FFL | n4 | is_sustained_osc | 20 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 115 | I1_miRNA_FFL | n4 | is_damped_osc | 20 | 0.00 | 0.00 (gx x0.5) | 0.00 (gx x0.5) | 0.00 | 0.00 | False | False |
| 116 | I1_miRNA_FFL | n4 | is_bistable | 20 | 2.39 | 1.02 (kappa x0.5) | 3.80 (kappa x2) | 0.26 | 1.40 | False | False |
| 117 | I1_miRNA_FFL | n4 | bistability_unresolved | 20 | 6.02 | 0.64 (gp x2) | 20.24 (gp x0.5) | 0.25 | 14.22 | False | False |
| 118 | I1_miRNA_FFL | n4 | is_ultrasensitive | 20 | 71.96 | 58.60 (gp x0.5) | 82.71 (gp x2) | 1.42 | 13.36 | False | False |
| 119 | I1_miRNA_FFL | n4 | is_very_ultrasensitive | 20 | 24.00 | 17.52 (gp x0.5) | 30.23 (gp x2) | 2.47 | 6.48 | False | False |
| 120 | I1_miRNA_FFL | n4 | is_noise_rejecting | 20 | 81.37 | 74.99 (Kxy x0.5) | 86.69 (Kxy x2) | 2.05 | 6.37 | False | False |
| 121 | I1_miRNA_FFL | n4 | nonmonotone_dose | 20 | 27.82 | 22.17 (Kxy x2) | 35.91 (Kxy x0.5) | 3.67 | 8.09 | False | False |

## (d) Parameter table (`parameter_table.csv`)

| line | parameter | meaning | low | high | distribution | acts in |
|---|---|---|---|---|---|---|
| 2 | gx | TF protein turnover relative to target mRNA | 0.5 | 5.0 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 3 | gm | mature miRNA turnover relative to target mRNA | 0.1 | 1.5 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 4 | gp | target protein turnover relative to target mRNA | 0.1 | 2.0 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 5 | gx2 | second TF turnover (n=6) | 0.5 | 5.0 | log-uniform | I1_miRNA_FFL:n6;COMP_C2_toggle:n6;COMP_I1_negfeedback:n6 |
| 6 | gm2 | second miRNA turnover (n=5) | 0.1 | 1.5 | log-uniform | I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6 |
| 7 | KS | signal level at half-maximal TF activation | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 8 | nS | cooperativity of signal -> TF | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 9 | Kxm | TF level at half-maximal miRNA promoter activity | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 10 | nxm | cooperativity of TF -> miRNA | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 11 | Kxy | TF level at half-maximal target promoter activity | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 12 | nxy | cooperativity of TF -> target | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 13 | Kx2y | TF2 level at half-maximal target promoter activity (n=6) | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n6;COMP_C2_toggle:n6;COMP_I1_negfeedback:n6 |
| 14 | nx2y | cooperativity of TF2 -> target (n=6) | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n6;COMP_C2_toggle:n6;COMP_I1_negfeedback:n6 |
| 15 | Kx2m | TF2 level at half-maximal miRNA promoter activity (n=6) | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n6;COMP_C2_toggle:n6;COMP_I1_negfeedback:n6 |
| 16 | nx2m | cooperativity of TF2 -> miRNA (n=6) | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n6;COMP_C2_toggle:n6;COMP_I1_negfeedback:n6 |
| 17 | Kx12 | TF1 level at half-maximal TF2 promoter activity (n=6) | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n6;COMP_C2_toggle:n6;COMP_I1_negfeedback:n6 |
| 18 | nx12 | cooperativity of TF1 -> TF2 (n=6) | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n6;COMP_C2_toggle:n6;COMP_I1_negfeedback:n6 |
| 19 | Kmx | miRNA level at half-maximal repression of the TF (composite) | 0.05 | 1.5 | log-uniform | COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 20 | nmx | cooperativity of miRNA -| TF | 1.0 | 4.0 | uniform | COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 21 | Kmy | miRNA level at half-maximal acceleration of target-mRNA decay | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 22 | nmy | cooperativity of RISC action on mRNA stability | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 23 | Kmp | miRNA level at half-maximal translational repression | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 24 | nmp | cooperativity of translational repression | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 25 | lam_y | max fold-acceleration of target-mRNA decay by miRNA | 0.5 | 20.0 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 26 | lam_x | max fold-acceleration of TF decay by miRNA (composite) | 0.5 | 20.0 | log-uniform | COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 27 | theta | stoichiometric titration rate (0 = purely catalytic miRNA) | 0.0 | 10.0 | uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 28 | Kc | shared RISC/AGO loading capacity | 0.2 | 50.0 | log-uniform | I1_miRNA_FFL:n3;I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n3;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n3;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 29 | b | TF-independent basal transcription of the target (OR gate) | 0.02 | 0.5 | log-uniform | none of the modules (OR-gate parameter) |
| 30 | rho | miRNA2 : miRNA1 production from the shared promoter (n=5) | 0.2 | 3.0 | log-uniform | I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6 |
| 31 | phi | miRNA2 : miRNA1 targeting efficacy (n=5) | 0.2 | 3.0 | log-uniform | I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6 |
| 32 | rho_g2 | G2 : G1 transcription rate (n>=4) | 0.2 | 3.0 | log-uniform | I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6 |
| 33 | kappa | max fold protection of G1/G2 protein by mutual complex formation (n>=4 reciprocal); saturable, so protection <= 1+kappa | 0.0 | 10.0 | uniform | I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6 |
| 34 | w2 | G2 mRNA weight in miRNA titration | 0.2 | 3.0 | log-uniform | I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 35 | Kg12 | partner-protein level at half-maximal effect of the gene-gene edge (half-maximal G2 transcription in the directed variant; half-maximal mutual stabilisation in the reciprocal variant) | 0.05 | 1.5 | log-uniform | I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |
| 36 | ng12 | cooperativity of the gene-gene edge (directed transcription or reciprocal stabilisation) | 1.0 | 4.0 | uniform | I1_miRNA_FFL:n4;I1_miRNA_FFL:n5;I1_miRNA_FFL:n6;I1_miRNA_FFL:n4d;COMP_C2_toggle:n4;COMP_C2_toggle:n5;COMP_C2_toggle:n6;COMP_C2_toggle:n4d;COMP_I1_negfeedback:n4;COMP_I1_negfeedback:n5;COMP_I1_negfeedback:n6;COMP_I1_negfeedback:n4d |

Sampler (every row): scrambled Sobol (scipy.stats.qmc.Sobol, d=35, scramble=True, seed=20260908), first 16,384 points (random_base2(14)); u in [0,1) mapped as lo*(hi/lo)^u (log) or lo+u*(hi-lo) (uniform).

