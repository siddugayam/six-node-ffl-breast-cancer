# VALUES: tier- and sign-aware null models for the three-node FFL test

Lines count the header as line 1. Tables are written by `summarise_tier_sign_nulls.py` from the replicate files in
`runs/` (written by `motif_null_tier_sign.c`, run by `run_tier_sign_nulls.sh`). 1,000 randomised networks per run.

## 1. Controls (`control_reproduction.csv`, column identical)

- line 2: analysed, published, the degree-preserving null: identical to `analyses/census_and_motif_nulls/nulls/null_A_dep_noexemplar.tsv` = True (1002 lines)
- line 3: analysed, published, the reciprocity-preserving null: identical to `analyses/census_and_motif_nulls/nulls/null_B_dep_noexemplar.tsv` = True (1002 lines)
- line 4: analysed, published, the bipartite-layer null: identical to `analyses/census_and_motif_nulls/nulls/null_C_dep_noexemplar.tsv` = True (1002 lines)
- line 5: census, published, the degree-preserving null: identical to `analyses/census_and_motif_nulls/nulls/null_A_pub_noexemplar.tsv` = True (1002 lines)
- line 6: census, published, the reciprocity-preserving null: identical to `analyses/census_and_motif_nulls/nulls/null_B_pub_noexemplar.tsv` = True (1002 lines)
- line 7: census, published, the bipartite-layer null: identical to `analyses/census_and_motif_nulls/nulls/null_C_pub_noexemplar.tsv` = True (1002 lines)
- line 8: analysed, retyped, the degree-preserving null: identical to `analyses/six_node_pattern/E/E5/null_A_dep_noexemplar_retyped22.tsv` = True (1002 lines)
- line 9: analysed, retyped, the reciprocity-preserving null: identical to `analyses/six_node_pattern/E/E5/null_B_dep_noexemplar_retyped22.tsv` = True (1002 lines)
- line 10: analysed, retyped, the bipartite-layer null: identical to `analyses/six_node_pattern/E/E5/null_C_dep_noexemplar_retyped22.tsv` = True (1002 lines)

## 2. Arcs per tier and sign (`arc_attribute_counts.csv`, column arcs)

| line | graph | typing | attribute | value | arc class | arcs |
|---|---|---|---|---|---|---|
| 2 | analysed | published | sign | + | TF->Gene | 267 |
| 3 | analysed | published | sign | + | TF->TF | 34 |
| 4 | analysed | published | sign | - | TF->Gene | 113 |
| 5 | analysed | published | sign | - | TF->TF | 23 |
| 6 | analysed | published | sign | unannotated | TF->Gene | 301 |
| 7 | analysed | published | sign | unannotated | TF->TF | 32 |
| 8 | analysed | published | tier | strong | miRNA->Gene | 404 |
| 9 | analysed | published | tier | strong | miRNA->TF | 303 |
| 10 | analysed | published | tier | weak | miRNA->Gene | 1677 |
| 11 | analysed | published | tier | weak | miRNA->TF | 632 |
| 12 | analysed | published | tier | predicted_only | miRNA->Gene | 1185 |
| 13 | analysed | published | tier | predicted_only | miRNA->TF | 618 |
| 14 | analysed | retyped | sign | + | TF->Gene | 247 |
| 15 | analysed | retyped | sign | + | TF->TF | 30 |
| 16 | analysed | retyped | sign | - | TF->Gene | 88 |
| 17 | analysed | retyped | sign | - | TF->TF | 13 |
| 18 | analysed | retyped | sign | unannotated | TF->Gene | 285 |
| 19 | analysed | retyped | sign | unannotated | TF->TF | 30 |
| 20 | analysed | retyped | tier | strong | miRNA->Gene | 428 |
| 21 | analysed | retyped | tier | strong | miRNA->TF | 279 |
| 22 | analysed | retyped | tier | weak | miRNA->Gene | 1728 |
| 23 | analysed | retyped | tier | weak | miRNA->TF | 581 |
| 24 | analysed | retyped | tier | predicted_only | miRNA->Gene | 1224 |
| 25 | analysed | retyped | tier | predicted_only | miRNA->TF | 579 |
| 26 | census | published | sign | + | TF->Gene | 267 |
| 27 | census | published | sign | + | TF->TF | 194 |
| 28 | census | published | sign | - | TF->Gene | 113 |
| 29 | census | published | sign | - | TF->TF | 163 |
| 30 | census | published | sign | unannotated | TF->Gene | 350 |
| 31 | census | published | sign | unannotated | TF->TF | 837 |
| 32 | census | published | tier | strong | miRNA->Gene | 404 |
| 33 | census | published | tier | strong | miRNA->TF | 303 |
| 34 | census | published | tier | weak | miRNA->Gene | 1677 |
| 35 | census | published | tier | weak | miRNA->TF | 632 |
| 36 | census | published | tier | predicted_only | miRNA->Gene | 1185 |
| 37 | census | published | tier | predicted_only | miRNA->TF | 618 |
| 38 | census | retyped | sign | + | TF->Gene | 264 |
| 39 | census | retyped | sign | + | TF->TF | 152 |
| 40 | census | retyped | sign | - | TF->Gene | 101 |
| 41 | census | retyped | sign | - | TF->TF | 100 |
| 42 | census | retyped | sign | unannotated | TF->Gene | 412 |
| 43 | census | retyped | sign | unannotated | TF->TF | 492 |
| 44 | census | retyped | tier | strong | miRNA->Gene | 428 |
| 45 | census | retyped | tier | strong | miRNA->TF | 279 |
| 46 | census | retyped | tier | weak | miRNA->Gene | 1728 |
| 47 | census | retyped | tier | weak | miRNA->TF | 581 |
| 48 | census | retyped | tier | predicted_only | miRNA->Gene | 1224 |
| 49 | census | retyped | tier | predicted_only | miRNA->TF | 579 |

Tiers are counted over all miRNA -> {Gene, TF} arcs. The reciprocal arcs keep their tier, but they are frozen and
never swapped.

## 3. Every run (`null_model_summary.csv`)

Columns observed, expected (mean of 1,000), sd, fold, excess_pct, p_upper, p_lower.

| line | graph | typing | null | metric | observed | expected (SD) | fold | excess % | p upper | p lower |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | analysed | published | Degree-preserving null | all three-node FFLs | 1774 | 2638.1 (80.6) | 0.672 | -32.8 | 1.0000 | 0.0010 |
| 3 | analysed | published | Degree-preserving null | composite FFLs | 1296 | 398.9 (34.1) | 3.249 | +224.9 | 0.0010 | 1.0000 |
| 4 | analysed | published | Degree-preserving null | TF-FFLs | 9 | 934.9 (42.6) | 0.010 | -99.0 | 1.0000 | 0.0010 |
| 5 | analysed | published | Degree-preserving null | miRNA-FFLs | 202 | 1141.9 (47.9) | 0.177 | -82.3 | 1.0000 | 0.0010 |
| 6 | analysed | published | Reciprocity-preserving null | all three-node FFLs | 1774 | 1644.3 (53.5) | 1.079 | +7.9 | 0.0130 | 0.9880 |
| 7 | analysed | published | Reciprocity-preserving null | composite FFLs | 1296 | 1217.3 (44.3) | 1.065 | +6.5 | 0.0430 | 0.9590 |
| 8 | analysed | published | Reciprocity-preserving null | TF-FFLs | 9 | 16.6 (4.8) | 0.542 | -45.8 | 0.9650 | 0.0569 |
| 9 | analysed | published | Reciprocity-preserving null | miRNA-FFLs | 202 | 232.4 (19.7) | 0.869 | -13.1 | 0.9451 | 0.0629 |
| 10 | analysed | published | Bipartite-layer null | all three-node FFLs | 1774 | 1595.1 (56.3) | 1.112 | +11.2 | 0.0010 | 1.0000 |
| 11 | analysed | published | Bipartite-layer null | composite FFLs | 1296 | 1201.6 (45.9) | 1.079 | +7.9 | 0.0270 | 0.9750 |
| 12 | analysed | published | Bipartite-layer null | TF-FFLs | 9 | 15.8 (4.8) | 0.570 | -43.0 | 0.9550 | 0.0809 |
| 13 | analysed | published | Bipartite-layer null | miRNA-FFLs | 202 | 188.6 (19.1) | 1.071 | +7.1 | 0.2607 | 0.7602 |
| 14 | analysed | published | (i) the reciprocity-preserving null + tier | all three-node FFLs | 1774 | 1554.7 (51.4) | 1.141 | +14.1 | 0.0010 | 1.0000 |
| 15 | analysed | published | (i) the reciprocity-preserving null + tier | composite FFLs | 1296 | 1183.3 (45.0) | 1.095 | +9.5 | 0.0060 | 0.9960 |
| 16 | analysed | published | (i) the reciprocity-preserving null + tier | TF-FFLs | 9 | 14.2 (4.5) | 0.635 | -36.5 | 0.9051 | 0.1439 |
| 17 | analysed | published | (i) the reciprocity-preserving null + tier | miRNA-FFLs | 202 | 177.7 (18.3) | 1.137 | +13.7 | 0.1059 | 0.9081 |
| 18 | analysed | published | (ii) the reciprocity-preserving null + sign | all three-node FFLs | 1774 | 1644.0 (53.4) | 1.079 | +7.9 | 0.0140 | 0.9880 |
| 19 | analysed | published | (ii) the reciprocity-preserving null + sign | composite FFLs | 1296 | 1208.1 (44.0) | 1.073 | +7.3 | 0.0260 | 0.9750 |
| 20 | analysed | published | (ii) the reciprocity-preserving null + sign | TF-FFLs | 9 | 16.7 (4.5) | 0.540 | -46.0 | 0.9750 | 0.0480 |
| 21 | analysed | published | (ii) the reciprocity-preserving null + sign | miRNA-FFLs | 202 | 231.7 (19.7) | 0.872 | -12.8 | 0.9431 | 0.0659 |
| 22 | analysed | published | (iii) the reciprocity-preserving null + tier + sign | all three-node FFLs | 1774 | 1555.3 (54.6) | 1.141 | +14.1 | 0.0010 | 1.0000 |
| 23 | analysed | published | (iii) the reciprocity-preserving null + tier + sign | composite FFLs | 1296 | 1174.8 (45.4) | 1.103 | +10.3 | 0.0050 | 0.9970 |
| 24 | analysed | published | (iii) the reciprocity-preserving null + tier + sign | TF-FFLs | 9 | 14.0 (4.3) | 0.644 | -35.6 | 0.8981 | 0.1638 |
| 25 | analysed | published | (iii) the reciprocity-preserving null + tier + sign | miRNA-FFLs | 202 | 177.7 (18.9) | 1.136 | +13.6 | 0.1149 | 0.8931 |
| 26 | analysed | published | Bipartite-layer null + tier | all three-node FFLs | 1774 | 1563.3 (56.2) | 1.135 | +13.5 | 0.0020 | 0.9990 |
| 27 | analysed | published | Bipartite-layer null + tier | composite FFLs | 1296 | 1184.1 (44.3) | 1.095 | +9.5 | 0.0050 | 0.9960 |
| 28 | analysed | published | Bipartite-layer null + tier | TF-FFLs | 9 | 14.4 (4.4) | 0.626 | -37.4 | 0.9131 | 0.1239 |
| 29 | analysed | published | Bipartite-layer null + tier | miRNA-FFLs | 202 | 175.9 (18.6) | 1.149 | +14.9 | 0.0949 | 0.9181 |
| 30 | analysed | published | Bipartite-layer null + tier + sign | all three-node FFLs | 1774 | 1477.2 (55.2) | 1.201 | +20.1 | 0.0010 | 1.0000 |
| 31 | analysed | published | Bipartite-layer null + tier + sign | composite FFLs | 1296 | 1112.7 (44.7) | 1.165 | +16.5 | 0.0010 | 1.0000 |
| 32 | analysed | published | Bipartite-layer null + tier + sign | TF-FFLs | 9 | 13.4 (4.3) | 0.673 | -32.7 | 0.8821 | 0.1798 |
| 33 | analysed | published | Bipartite-layer null + tier + sign | miRNA-FFLs | 202 | 168.3 (18.0) | 1.200 | +20.0 | 0.0360 | 0.9680 |
| 34 | analysed | retyped | Degree-preserving null | all three-node FFLs | 1687 | 2650.0 (83.2) | 0.637 | -36.3 | 1.0000 | 0.0010 |
| 35 | analysed | retyped | Degree-preserving null | composite FFLs | 1240 | 390.7 (35.0) | 3.174 | +217.4 | 0.0010 | 1.0000 |
| 36 | analysed | retyped | Degree-preserving null | TF-FFLs | 9 | 882.2 (42.0) | 0.010 | -99.0 | 1.0000 | 0.0010 |
| 37 | analysed | retyped | Degree-preserving null | miRNA-FFLs | 199 | 1085.4 (46.4) | 0.183 | -81.7 | 1.0000 | 0.0010 |
| 38 | analysed | retyped | Reciprocity-preserving null | all three-node FFLs | 1687 | 1712.9 (55.8) | 0.985 | -1.5 | 0.6783 | 0.3337 |
| 39 | analysed | retyped | Reciprocity-preserving null | composite FFLs | 1240 | 1179.8 (44.4) | 1.051 | +5.1 | 0.0909 | 0.9131 |
| 40 | analysed | retyped | Reciprocity-preserving null | TF-FFLs | 9 | 15.0 (4.6) | 0.599 | -40.1 | 0.9361 | 0.1049 |
| 41 | analysed | retyped | Reciprocity-preserving null | miRNA-FFLs | 199 | 223.6 (18.5) | 0.890 | -11.0 | 0.9121 | 0.0969 |
| 42 | analysed | retyped | Bipartite-layer null | all three-node FFLs | 1687 | 1663.6 (55.1) | 1.014 | +1.4 | 0.3407 | 0.6633 |
| 43 | analysed | retyped | Bipartite-layer null | composite FFLs | 1240 | 1164.3 (42.8) | 1.065 | +6.5 | 0.0380 | 0.9650 |
| 44 | analysed | retyped | Bipartite-layer null | TF-FFLs | 9 | 14.1 (4.4) | 0.639 | -36.1 | 0.9161 | 0.1499 |
| 45 | analysed | retyped | Bipartite-layer null | miRNA-FFLs | 199 | 180.6 (19.2) | 1.102 | +10.2 | 0.1788 | 0.8302 |
| 46 | analysed | retyped | (i) the reciprocity-preserving null + tier | all three-node FFLs | 1687 | 1626.6 (54.5) | 1.037 | +3.7 | 0.1389 | 0.8661 |
| 47 | analysed | retyped | (i) the reciprocity-preserving null + tier | composite FFLs | 1240 | 1150.6 (42.2) | 1.078 | +7.8 | 0.0200 | 0.9820 |
| 48 | analysed | retyped | (i) the reciprocity-preserving null + tier | TF-FFLs | 9 | 13.1 (4.3) | 0.688 | -31.2 | 0.8621 | 0.2058 |
| 49 | analysed | retyped | (i) the reciprocity-preserving null + tier | miRNA-FFLs | 199 | 170.9 (18.3) | 1.164 | +16.4 | 0.0709 | 0.9381 |
| 50 | analysed | retyped | (ii) the reciprocity-preserving null + sign | all three-node FFLs | 1687 | 1719.9 (54.8) | 0.981 | -1.9 | 0.7303 | 0.2797 |
| 51 | analysed | retyped | (ii) the reciprocity-preserving null + sign | composite FFLs | 1240 | 1176.9 (41.3) | 1.054 | +5.4 | 0.0649 | 0.9381 |
| 52 | analysed | retyped | (ii) the reciprocity-preserving null + sign | TF-FFLs | 9 | 14.9 (4.3) | 0.604 | -39.6 | 0.9361 | 0.1019 |
| 53 | analysed | retyped | (ii) the reciprocity-preserving null + sign | miRNA-FFLs | 199 | 224.7 (19.7) | 0.886 | -11.4 | 0.9161 | 0.0919 |
| 54 | analysed | retyped | (iii) the reciprocity-preserving null + tier + sign | all three-node FFLs | 1687 | 1628.5 (54.1) | 1.036 | +3.6 | 0.1419 | 0.8621 |
| 55 | analysed | retyped | (iii) the reciprocity-preserving null + tier + sign | composite FFLs | 1240 | 1146.6 (44.0) | 1.081 | +8.1 | 0.0200 | 0.9810 |
| 56 | analysed | retyped | (iii) the reciprocity-preserving null + tier + sign | TF-FFLs | 9 | 12.5 (4.0) | 0.721 | -27.9 | 0.8511 | 0.2348 |
| 57 | analysed | retyped | (iii) the reciprocity-preserving null + tier + sign | miRNA-FFLs | 199 | 169.7 (17.4) | 1.172 | +17.2 | 0.0470 | 0.9560 |
| 58 | analysed | retyped | Bipartite-layer null + tier | all three-node FFLs | 1687 | 1632.8 (58.3) | 1.033 | +3.3 | 0.1828 | 0.8262 |
| 59 | analysed | retyped | Bipartite-layer null + tier | composite FFLs | 1240 | 1149.0 (45.2) | 1.079 | +7.9 | 0.0230 | 0.9800 |
| 60 | analysed | retyped | Bipartite-layer null + tier | TF-FFLs | 9 | 12.9 (4.2) | 0.699 | -30.1 | 0.8571 | 0.2328 |
| 61 | analysed | retyped | Bipartite-layer null + tier | miRNA-FFLs | 199 | 169.3 (18.2) | 1.175 | +17.5 | 0.0509 | 0.9560 |
| 62 | analysed | retyped | Bipartite-layer null + tier + sign | all three-node FFLs | 1687 | 1549.5 (58.9) | 1.089 | +8.9 | 0.0120 | 0.9900 |
| 63 | analysed | retyped | Bipartite-layer null + tier + sign | composite FFLs | 1240 | 1081.2 (45.6) | 1.147 | +14.7 | 0.0010 | 1.0000 |
| 64 | analysed | retyped | Bipartite-layer null + tier + sign | TF-FFLs | 9 | 12.1 (4.2) | 0.744 | -25.6 | 0.8012 | 0.2857 |
| 65 | analysed | retyped | Bipartite-layer null + tier + sign | miRNA-FFLs | 199 | 162.5 (17.3) | 1.225 | +22.5 | 0.0200 | 0.9810 |
| 66 | census | published | Degree-preserving null | all three-node FFLs | 5833 | 9702.6 (257.6) | 0.601 | -39.9 | 1.0000 | 0.0010 |
| 67 | census | published | Degree-preserving null | composite FFLs | 1375 | 702.4 (48.2) | 1.958 | +95.8 | 0.0010 | 1.0000 |
| 68 | census | published | Degree-preserving null | TF-FFLs | 9 | 1793.4 (70.2) | 0.005 | -99.5 | 1.0000 | 0.0010 |
| 69 | census | published | Degree-preserving null | miRNA-FFLs | 289 | 2148.4 (73.6) | 0.135 | -86.5 | 1.0000 | 0.0010 |
| 70 | census | published | Reciprocity-preserving null | all three-node FFLs | 5833 | 6778.1 (188.5) | 0.861 | -13.9 | 1.0000 | 0.0010 |
| 71 | census | published | Reciprocity-preserving null | composite FFLs | 1375 | 1538.0 (53.0) | 0.894 | -10.6 | 1.0000 | 0.0010 |
| 72 | census | published | Reciprocity-preserving null | TF-FFLs | 9 | 20.2 (5.4) | 0.446 | -55.4 | 0.9910 | 0.0170 |
| 73 | census | published | Reciprocity-preserving null | miRNA-FFLs | 289 | 320.7 (22.6) | 0.901 | -9.9 | 0.9331 | 0.0739 |
| 74 | census | published | Bipartite-layer null | all three-node FFLs | 5833 | 6506.0 (193.0) | 0.897 | -10.3 | 1.0000 | 0.0010 |
| 75 | census | published | Bipartite-layer null | composite FFLs | 1375 | 1418.9 (56.9) | 0.969 | -3.1 | 0.7772 | 0.2318 |
| 76 | census | published | Bipartite-layer null | TF-FFLs | 9 | 20.3 (5.5) | 0.443 | -55.7 | 0.9900 | 0.0180 |
| 77 | census | published | Bipartite-layer null | miRNA-FFLs | 289 | 232.4 (23.4) | 1.243 | +24.3 | 0.0090 | 0.9930 |
| 78 | census | published | (i) the reciprocity-preserving null + tier | all three-node FFLs | 5833 | 6526.2 (190.3) | 0.894 | -10.6 | 1.0000 | 0.0010 |
| 79 | census | published | (i) the reciprocity-preserving null + tier | composite FFLs | 1375 | 1493.5 (52.5) | 0.921 | -7.9 | 0.9850 | 0.0160 |
| 80 | census | published | (i) the reciprocity-preserving null + tier | TF-FFLs | 9 | 19.1 (5.2) | 0.470 | -53.0 | 0.9900 | 0.0210 |
| 81 | census | published | (i) the reciprocity-preserving null + tier | miRNA-FFLs | 289 | 228.2 (21.8) | 1.266 | +26.6 | 0.0070 | 0.9940 |
| 82 | census | published | (ii) the reciprocity-preserving null + sign | all three-node FFLs | 5833 | 6793.1 (188.5) | 0.859 | -14.1 | 1.0000 | 0.0010 |
| 83 | census | published | (ii) the reciprocity-preserving null + sign | composite FFLs | 1375 | 1535.5 (51.6) | 0.895 | -10.5 | 1.0000 | 0.0010 |
| 84 | census | published | (ii) the reciprocity-preserving null + sign | TF-FFLs | 9 | 20.4 (5.5) | 0.442 | -55.8 | 0.9940 | 0.0110 |
| 85 | census | published | (ii) the reciprocity-preserving null + sign | miRNA-FFLs | 289 | 321.6 (24.0) | 0.899 | -10.1 | 0.9061 | 0.1019 |
| 86 | census | published | (iii) the reciprocity-preserving null + tier + sign | all three-node FFLs | 5833 | 6545.9 (197.8) | 0.891 | -10.9 | 1.0000 | 0.0010 |
| 87 | census | published | (iii) the reciprocity-preserving null + tier + sign | composite FFLs | 1375 | 1489.6 (50.7) | 0.923 | -7.7 | 0.9930 | 0.0080 |
| 88 | census | published | (iii) the reciprocity-preserving null + tier + sign | TF-FFLs | 9 | 19.3 (5.1) | 0.467 | -53.3 | 0.9860 | 0.0250 |
| 89 | census | published | (iii) the reciprocity-preserving null + tier + sign | miRNA-FFLs | 289 | 228.0 (20.6) | 1.268 | +26.8 | 0.0020 | 0.9990 |
| 90 | census | published | Bipartite-layer null + tier | all three-node FFLs | 5833 | 6380.3 (182.5) | 0.914 | -8.6 | 0.9990 | 0.0020 |
| 91 | census | published | Bipartite-layer null + tier | composite FFLs | 1375 | 1403.4 (57.2) | 0.980 | -2.0 | 0.6893 | 0.3207 |
| 92 | census | published | Bipartite-layer null + tier | TF-FFLs | 9 | 18.9 (5.3) | 0.477 | -52.3 | 0.9830 | 0.0310 |
| 93 | census | published | Bipartite-layer null + tier | miRNA-FFLs | 289 | 216.1 (21.8) | 1.337 | +33.7 | 0.0050 | 0.9980 |
| 94 | census | published | Bipartite-layer null + tier + sign | all three-node FFLs | 5833 | 6306.8 (208.6) | 0.925 | -7.5 | 0.9930 | 0.0080 |
| 95 | census | published | Bipartite-layer null + tier + sign | composite FFLs | 1375 | 1359.3 (55.1) | 1.012 | +1.2 | 0.3866 | 0.6204 |
| 96 | census | published | Bipartite-layer null + tier + sign | TF-FFLs | 9 | 17.8 (5.2) | 0.504 | -49.6 | 0.9750 | 0.0370 |
| 97 | census | published | Bipartite-layer null + tier + sign | miRNA-FFLs | 289 | 213.7 (21.4) | 1.352 | +35.2 | 0.0020 | 1.0000 |
| 98 | census | retyped | Degree-preserving null | all three-node FFLs | 5352 | 9563.1 (255.6) | 0.560 | -44.0 | 1.0000 | 0.0010 |
| 99 | census | retyped | Degree-preserving null | composite FFLs | 1297 | 665.9 (46.6) | 1.948 | +94.8 | 0.0010 | 1.0000 |
| 100 | census | retyped | Degree-preserving null | TF-FFLs | 9 | 1591.9 (67.3) | 0.006 | -99.4 | 1.0000 | 0.0010 |
| 101 | census | retyped | Degree-preserving null | miRNA-FFLs | 282 | 1948.2 (72.0) | 0.145 | -85.5 | 1.0000 | 0.0010 |
| 102 | census | retyped | Reciprocity-preserving null | all three-node FFLs | 5352 | 6835.0 (191.2) | 0.783 | -21.7 | 1.0000 | 0.0010 |
| 103 | census | retyped | Reciprocity-preserving null | composite FFLs | 1297 | 1537.8 (51.3) | 0.843 | -15.7 | 1.0000 | 0.0010 |
| 104 | census | retyped | Reciprocity-preserving null | TF-FFLs | 9 | 18.2 (4.8) | 0.495 | -50.5 | 0.9930 | 0.0240 |
| 105 | census | retyped | Reciprocity-preserving null | miRNA-FFLs | 282 | 319.7 (22.5) | 0.882 | -11.8 | 0.9650 | 0.0390 |
| 106 | census | retyped | Bipartite-layer null | all three-node FFLs | 5352 | 6727.9 (189.0) | 0.795 | -20.5 | 1.0000 | 0.0010 |
| 107 | census | retyped | Bipartite-layer null | composite FFLs | 1297 | 1474.8 (56.0) | 0.879 | -12.1 | 1.0000 | 0.0010 |
| 108 | census | retyped | Bipartite-layer null | TF-FFLs | 9 | 18.9 (5.2) | 0.475 | -52.5 | 0.9870 | 0.0230 |
| 109 | census | retyped | Bipartite-layer null | miRNA-FFLs | 282 | 253.1 (23.4) | 1.114 | +11.4 | 0.1149 | 0.8941 |
| 110 | census | retyped | (i) the reciprocity-preserving null + tier | all three-node FFLs | 5352 | 6587.8 (197.2) | 0.812 | -18.8 | 1.0000 | 0.0010 |
| 111 | census | retyped | (i) the reciprocity-preserving null + tier | composite FFLs | 1297 | 1487.2 (52.7) | 0.872 | -12.8 | 1.0000 | 0.0010 |
| 112 | census | retyped | (i) the reciprocity-preserving null + tier | TF-FFLs | 9 | 17.8 (4.9) | 0.506 | -49.4 | 0.9810 | 0.0380 |
| 113 | census | retyped | (i) the reciprocity-preserving null + tier | miRNA-FFLs | 282 | 227.0 (21.6) | 1.242 | +24.2 | 0.0110 | 0.9920 |
| 114 | census | retyped | (ii) the reciprocity-preserving null + sign | all three-node FFLs | 5352 | 6844.5 (194.4) | 0.782 | -21.8 | 1.0000 | 0.0010 |
| 115 | census | retyped | (ii) the reciprocity-preserving null + sign | composite FFLs | 1297 | 1537.1 (52.8) | 0.844 | -15.6 | 1.0000 | 0.0010 |
| 116 | census | retyped | (ii) the reciprocity-preserving null + sign | TF-FFLs | 9 | 18.4 (5.0) | 0.489 | -51.1 | 0.9850 | 0.0330 |
| 117 | census | retyped | (ii) the reciprocity-preserving null + sign | miRNA-FFLs | 282 | 319.4 (22.5) | 0.883 | -11.7 | 0.9600 | 0.0440 |
| 118 | census | retyped | (iii) the reciprocity-preserving null + tier + sign | all three-node FFLs | 5352 | 6598.0 (186.3) | 0.811 | -18.9 | 1.0000 | 0.0010 |
| 119 | census | retyped | (iii) the reciprocity-preserving null + tier + sign | composite FFLs | 1297 | 1490.4 (51.0) | 0.870 | -13.0 | 1.0000 | 0.0010 |
| 120 | census | retyped | (iii) the reciprocity-preserving null + tier + sign | TF-FFLs | 9 | 17.4 (5.0) | 0.518 | -48.2 | 0.9790 | 0.0420 |
| 121 | census | retyped | (iii) the reciprocity-preserving null + tier + sign | miRNA-FFLs | 282 | 229.3 (21.2) | 1.230 | +23.0 | 0.0120 | 0.9900 |
| 122 | census | retyped | Bipartite-layer null + tier | all three-node FFLs | 5352 | 6580.5 (202.2) | 0.813 | -18.7 | 1.0000 | 0.0010 |
| 123 | census | retyped | Bipartite-layer null + tier | composite FFLs | 1297 | 1455.5 (54.4) | 0.891 | -10.9 | 0.9990 | 0.0020 |
| 124 | census | retyped | Bipartite-layer null + tier | TF-FFLs | 9 | 17.8 (5.2) | 0.504 | -49.6 | 0.9780 | 0.0380 |
| 125 | census | retyped | Bipartite-layer null + tier | miRNA-FFLs | 282 | 234.4 (22.6) | 1.203 | +20.3 | 0.0220 | 0.9810 |
| 126 | census | retyped | Bipartite-layer null + tier + sign | all three-node FFLs | 5352 | 6518.2 (203.5) | 0.821 | -17.9 | 1.0000 | 0.0010 |
| 127 | census | retyped | Bipartite-layer null + tier + sign | composite FFLs | 1297 | 1404.0 (53.7) | 0.924 | -7.6 | 0.9820 | 0.0210 |
| 128 | census | retyped | Bipartite-layer null + tier + sign | TF-FFLs | 9 | 17.2 (4.8) | 0.523 | -47.7 | 0.9760 | 0.0400 |
| 129 | census | retyped | Bipartite-layer null + tier + sign | miRNA-FFLs | 282 | 228.9 (22.9) | 1.232 | +23.2 | 0.0120 | 0.9900 |

## 4. Waterfall of the residual excess (`waterfall.csv`, columns excess_pct and change_pp)

| line | graph | typing | metric | step | excess % | change (points) |
|---|---|---|---|---|---|---|
| 2 | analysed | published | all three-node FFLs | Degree-preserving null | -32.8 |  |
| 3 | analysed | published | all three-node FFLs | + reciprocity (the reciprocity-preserving null) | +7.9 | +40.6 |
| 4 | analysed | published | all three-node FFLs | + bipartite curveball (the bipartite-layer null) | +11.2 | +3.3 |
| 5 | analysed | published | all three-node FFLs | + tier (the bipartite-layer null + tier) | +13.5 | +2.3 |
| 6 | analysed | published | all three-node FFLs | + sign (the bipartite-layer null + tier + sign) | +20.1 | +6.6 |
| 7 | analysed | published | composite FFLs | Degree-preserving null | +224.9 |  |
| 8 | analysed | published | composite FFLs | + reciprocity (the reciprocity-preserving null) | +6.5 | -218.4 |
| 9 | analysed | published | composite FFLs | + bipartite curveball (the bipartite-layer null) | +7.9 | +1.4 |
| 10 | analysed | published | composite FFLs | + tier (the bipartite-layer null + tier) | +9.5 | +1.6 |
| 11 | analysed | published | composite FFLs | + sign (the bipartite-layer null + tier + sign) | +16.5 | +7.0 |
| 12 | analysed | published | TF-FFLs | Degree-preserving null | -99.0 |  |
| 13 | analysed | published | TF-FFLs | + reciprocity (the reciprocity-preserving null) | -45.8 | +53.3 |
| 14 | analysed | published | TF-FFLs | + bipartite curveball (the bipartite-layer null) | -43.0 | +2.8 |
| 15 | analysed | published | TF-FFLs | + tier (the bipartite-layer null + tier) | -37.4 | +5.6 |
| 16 | analysed | published | TF-FFLs | + sign (the bipartite-layer null + tier + sign) | -32.7 | +4.6 |
| 17 | analysed | published | miRNA-FFLs | Degree-preserving null | -82.3 |  |
| 18 | analysed | published | miRNA-FFLs | + reciprocity (the reciprocity-preserving null) | -13.1 | +69.2 |
| 19 | analysed | published | miRNA-FFLs | + bipartite curveball (the bipartite-layer null) | +7.1 | +20.2 |
| 20 | analysed | published | miRNA-FFLs | + tier (the bipartite-layer null + tier) | +14.9 | +7.7 |
| 21 | analysed | published | miRNA-FFLs | + sign (the bipartite-layer null + tier + sign) | +20.0 | +5.1 |
| 22 | analysed | retyped | all three-node FFLs | Degree-preserving null | -36.3 |  |
| 23 | analysed | retyped | all three-node FFLs | + reciprocity (the reciprocity-preserving null) | -1.5 | +34.8 |
| 24 | analysed | retyped | all three-node FFLs | + bipartite curveball (the bipartite-layer null) | +1.4 | +2.9 |
| 25 | analysed | retyped | all three-node FFLs | + tier (the bipartite-layer null + tier) | +3.3 | +1.9 |
| 26 | analysed | retyped | all three-node FFLs | + sign (the bipartite-layer null + tier + sign) | +8.9 | +5.6 |
| 27 | analysed | retyped | composite FFLs | Degree-preserving null | +217.4 |  |
| 28 | analysed | retyped | composite FFLs | + reciprocity (the reciprocity-preserving null) | +5.1 | -212.3 |
| 29 | analysed | retyped | composite FFLs | + bipartite curveball (the bipartite-layer null) | +6.5 | +1.4 |
| 30 | analysed | retyped | composite FFLs | + tier (the bipartite-layer null + tier) | +7.9 | +1.4 |
| 31 | analysed | retyped | composite FFLs | + sign (the bipartite-layer null + tier + sign) | +14.7 | +6.8 |
| 32 | analysed | retyped | TF-FFLs | Degree-preserving null | -99.0 |  |
| 33 | analysed | retyped | TF-FFLs | + reciprocity (the reciprocity-preserving null) | -40.1 | +58.8 |
| 34 | analysed | retyped | TF-FFLs | + bipartite curveball (the bipartite-layer null) | -36.1 | +4.1 |
| 35 | analysed | retyped | TF-FFLs | + tier (the bipartite-layer null + tier) | -30.1 | +5.9 |
| 36 | analysed | retyped | TF-FFLs | + sign (the bipartite-layer null + tier + sign) | -25.6 | +4.6 |
| 37 | analysed | retyped | miRNA-FFLs | Degree-preserving null | -81.7 |  |
| 38 | analysed | retyped | miRNA-FFLs | + reciprocity (the reciprocity-preserving null) | -11.0 | +70.6 |
| 39 | analysed | retyped | miRNA-FFLs | + bipartite curveball (the bipartite-layer null) | +10.2 | +21.2 |
| 40 | analysed | retyped | miRNA-FFLs | + tier (the bipartite-layer null + tier) | +17.5 | +7.3 |
| 41 | analysed | retyped | miRNA-FFLs | + sign (the bipartite-layer null + tier + sign) | +22.5 | +5.0 |
| 42 | census | published | all three-node FFLs | Degree-preserving null | -39.9 |  |
| 43 | census | published | all three-node FFLs | + reciprocity (the reciprocity-preserving null) | -13.9 | +25.9 |
| 44 | census | published | all three-node FFLs | + bipartite curveball (the bipartite-layer null) | -10.3 | +3.6 |
| 45 | census | published | all three-node FFLs | + tier (the bipartite-layer null + tier) | -8.6 | +1.8 |
| 46 | census | published | all three-node FFLs | + sign (the bipartite-layer null + tier + sign) | -7.5 | +1.1 |
| 47 | census | published | composite FFLs | Degree-preserving null | +95.8 |  |
| 48 | census | published | composite FFLs | + reciprocity (the reciprocity-preserving null) | -10.6 | -106.4 |
| 49 | census | published | composite FFLs | + bipartite curveball (the bipartite-layer null) | -3.1 | +7.5 |
| 50 | census | published | composite FFLs | + tier (the bipartite-layer null + tier) | -2.0 | +1.1 |
| 51 | census | published | composite FFLs | + sign (the bipartite-layer null + tier + sign) | +1.2 | +3.2 |
| 52 | census | published | TF-FFLs | Degree-preserving null | -99.5 |  |
| 53 | census | published | TF-FFLs | + reciprocity (the reciprocity-preserving null) | -55.4 | +44.1 |
| 54 | census | published | TF-FFLs | + bipartite curveball (the bipartite-layer null) | -55.7 | -0.3 |
| 55 | census | published | TF-FFLs | + tier (the bipartite-layer null + tier) | -52.3 | +3.4 |
| 56 | census | published | TF-FFLs | + sign (the bipartite-layer null + tier + sign) | -49.6 | +2.7 |
| 57 | census | published | miRNA-FFLs | Degree-preserving null | -86.5 |  |
| 58 | census | published | miRNA-FFLs | + reciprocity (the reciprocity-preserving null) | -9.9 | +76.7 |
| 59 | census | published | miRNA-FFLs | + bipartite curveball (the bipartite-layer null) | +24.3 | +34.2 |
| 60 | census | published | miRNA-FFLs | + tier (the bipartite-layer null + tier) | +33.7 | +9.4 |
| 61 | census | published | miRNA-FFLs | + sign (the bipartite-layer null + tier + sign) | +35.2 | +1.5 |
| 62 | census | retyped | all three-node FFLs | Degree-preserving null | -44.0 |  |
| 63 | census | retyped | all three-node FFLs | + reciprocity (the reciprocity-preserving null) | -21.7 | +22.3 |
| 64 | census | retyped | all three-node FFLs | + bipartite curveball (the bipartite-layer null) | -20.5 | +1.2 |
| 65 | census | retyped | all three-node FFLs | + tier (the bipartite-layer null + tier) | -18.7 | +1.8 |
| 66 | census | retyped | all three-node FFLs | + sign (the bipartite-layer null + tier + sign) | -17.9 | +0.8 |
| 67 | census | retyped | composite FFLs | Degree-preserving null | +94.8 |  |
| 68 | census | retyped | composite FFLs | + reciprocity (the reciprocity-preserving null) | -15.7 | -110.4 |
| 69 | census | retyped | composite FFLs | + bipartite curveball (the bipartite-layer null) | -12.1 | +3.6 |
| 70 | census | retyped | composite FFLs | + tier (the bipartite-layer null + tier) | -10.9 | +1.2 |
| 71 | census | retyped | composite FFLs | + sign (the bipartite-layer null + tier + sign) | -7.6 | +3.3 |
| 72 | census | retyped | TF-FFLs | Degree-preserving null | -99.4 |  |
| 73 | census | retyped | TF-FFLs | + reciprocity (the reciprocity-preserving null) | -50.5 | +48.9 |
| 74 | census | retyped | TF-FFLs | + bipartite curveball (the bipartite-layer null) | -52.5 | -2.0 |
| 75 | census | retyped | TF-FFLs | + tier (the bipartite-layer null + tier) | -49.6 | +2.9 |
| 76 | census | retyped | TF-FFLs | + sign (the bipartite-layer null + tier + sign) | -47.7 | +1.8 |
| 77 | census | retyped | miRNA-FFLs | Degree-preserving null | -85.5 |  |
| 78 | census | retyped | miRNA-FFLs | + reciprocity (the reciprocity-preserving null) | -11.8 | +73.7 |
| 79 | census | retyped | miRNA-FFLs | + bipartite curveball (the bipartite-layer null) | +11.4 | +23.2 |
| 80 | census | retyped | miRNA-FFLs | + tier (the bipartite-layer null + tier) | +20.3 | +8.9 |
| 81 | census | retyped | miRNA-FFLs | + sign (the bipartite-layer null + tier + sign) | +23.2 | +2.9 |
