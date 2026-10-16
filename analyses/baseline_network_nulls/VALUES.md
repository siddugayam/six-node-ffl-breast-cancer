# VALUES: a baseline network without FFL pre-selection

Lines count the header as line 1. Networks: `build_baseline_network.R` (log `runs/build_baseline_network.log`); tests:
`summarise_baseline_nulls.py` from the replicate files in `runs/`.

The curated lists from which the network was assembled (1,057 genes, 421 miRNAs) are not available. The substitute keeps the 587 nodes
of the analysed network and takes every edge among them from the primary tables.

## 1. Networks (`network_arc_counts.csv`)

| line | network | typing | arcs | reciprocal_TF_miRNA_pairs | Gene->Gene | Gene->miRNA | Gene->TF | miRNA->Gene | miRNA->TF | TF->Gene | TF->miRNA | TF->TF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | validated | published | 30604 | 809 | 15 | 16 | 21 | 15216 | 12641 | 681 | 1357 | 657 |
| 3 | validated | retyped | 30604 | 734 | 97 | 131 | 122 | 17422 | 10435 | 685 | 1242 | 470 |
| 4 | validated_plus_predicted | published | 32868 | 840 | 15 | 16 | 21 | 16347 | 13774 | 681 | 1357 | 657 |
| 5 | validated_plus_predicted | retyped | 32868 | 759 | 97 | 131 | 122 | 18656 | 11465 | 685 | 1242 | 470 |

Checks (`runs/build_baseline_network.log`):
- line 1: TRRUST arcs among the nodes: 1374 ; equal to the project's TRRUST layer: TRUE
- line 2: TransmiR arcs among the nodes: 1373 ; equal to the project's TransmiR layer: TRUE

## 2. Overlap with the analysed network (`analysed_network_overlap.csv`)

| line | network | edge type | analysed arcs | in the primary tables | baseline arcs |
|---|---|---|---|---|---|
| 2 | validated | gene_gene | 1 | 0 |  |
| 3 | validated | miRNA_target | 4819 | 3016 |  |
| 4 | validated | TF_miRNA | 1239 | 1236 |  |
| 5 | validated | TF_target | 770 | 770 |  |
| 6 | validated_plus_predicted | gene_gene | 1 | 0 |  |
| 7 | validated_plus_predicted | miRNA_target | 4819 | 4257 |  |
| 8 | validated_plus_predicted | TF_miRNA | 1239 | 1236 |  |
| 9 | validated_plus_predicted | TF_target | 770 | 770 |  |
| 10 | validated | all | 6829 | 5022 | 30604 |
| 11 | validated_plus_predicted | all | 6829 | 6263 | 32868 |

## 3. Three-node FFL motif test (`null_model_summary.csv`)

| line | network | typing | null model | metric | observed | expected (SD) | fold | excess % | p upper | p lower |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | validated | published typing (157 TFs) | degree-preserving | all three-node FFLs | 63890 | 66642.9 (671.7) | 0.959 | -4.1 | 1.0000 | 0.0010 |
| 3 | validated | published typing (157 TFs) | degree-preserving | composite FFLs | 8196 | 7593.2 (225.1) | 1.079 | +7.9 | 0.0040 | 0.9970 |
| 4 | validated | published typing (157 TFs) | degree-preserving | TF-FFLs | 2281 | 2797.2 (151.5) | 0.815 | -18.5 | 1.0000 | 0.0010 |
| 5 | validated | published typing (157 TFs) | degree-preserving | miRNA-FFLs | 50164 | 52742.6 (545.9) | 0.951 | -4.9 | 1.0000 | 0.0010 |
| 6 | validated | published typing (157 TFs) | reciprocity-preserving | all three-node FFLs | 63890 | 64053.5 (644.7) | 0.997 | -0.3 | 0.5994 | 0.4026 |
| 7 | validated | published typing (157 TFs) | reciprocity-preserving | composite FFLs | 8196 | 9314.9 (182.7) | 0.880 | -12.0 | 1.0000 | 0.0010 |
| 8 | validated | published typing (157 TFs) | reciprocity-preserving | TF-FFLs | 2281 | 1234.5 (81.8) | 1.848 | +84.8 | 0.0010 | 1.0000 |
| 9 | validated | published typing (157 TFs) | reciprocity-preserving | miRNA-FFLs | 50164 | 50013.7 (491.6) | 1.003 | +0.3 | 0.3746 | 0.6264 |
| 10 | validated | published typing (157 TFs) | bipartite-layer | all three-node FFLs | 63890 | 58351.5 (591.5) | 1.095 | +9.5 | 0.0010 | 1.0000 |
| 11 | validated | published typing (157 TFs) | bipartite-layer | composite FFLs | 8196 | 8728.9 (166.8) | 0.939 | -6.1 | 0.9980 | 0.0030 |
| 12 | validated | published typing (157 TFs) | bipartite-layer | TF-FFLs | 2281 | 1319.8 (83.7) | 1.728 | +72.8 | 0.0010 | 1.0000 |
| 13 | validated | published typing (157 TFs) | bipartite-layer | miRNA-FFLs | 50164 | 44839.6 (463.3) | 1.119 | +11.9 | 0.0010 | 1.0000 |
| 14 | validated_plus_predicted | published typing (157 TFs) | degree-preserving | all three-node FFLs | 68049 | 70643.6 (711.2) | 0.963 | -3.7 | 1.0000 | 0.0010 |
| 15 | validated_plus_predicted | published typing (157 TFs) | degree-preserving | composite FFLs | 8486 | 7811.8 (229.5) | 1.086 | +8.6 | 0.0040 | 0.9970 |
| 16 | validated_plus_predicted | published typing (157 TFs) | degree-preserving | TF-FFLs | 2353 | 2985.9 (160.4) | 0.788 | -21.2 | 1.0000 | 0.0010 |
| 17 | validated_plus_predicted | published typing (157 TFs) | degree-preserving | miRNA-FFLs | 53889 | 56256.7 (561.9) | 0.958 | -4.2 | 1.0000 | 0.0010 |
| 18 | validated_plus_predicted | published typing (157 TFs) | reciprocity-preserving | all three-node FFLs | 68049 | 67937.3 (673.7) | 1.002 | +0.2 | 0.4356 | 0.5664 |
| 19 | validated_plus_predicted | published typing (157 TFs) | reciprocity-preserving | composite FFLs | 8486 | 9666.5 (170.3) | 0.878 | -12.2 | 1.0000 | 0.0010 |
| 20 | validated_plus_predicted | published typing (157 TFs) | reciprocity-preserving | TF-FFLs | 2353 | 1285.9 (84.3) | 1.830 | +83.0 | 0.0010 | 1.0000 |
| 21 | validated_plus_predicted | published typing (157 TFs) | reciprocity-preserving | miRNA-FFLs | 53889 | 53424.1 (516.7) | 1.009 | +0.9 | 0.1938 | 0.8072 |
| 22 | validated_plus_predicted | published typing (157 TFs) | bipartite-layer | all three-node FFLs | 68049 | 62176.2 (614.2) | 1.094 | +9.4 | 0.0010 | 1.0000 |
| 23 | validated_plus_predicted | published typing (157 TFs) | bipartite-layer | composite FFLs | 8486 | 9056.6 (174.0) | 0.937 | -6.3 | 0.9990 | 0.0020 |
| 24 | validated_plus_predicted | published typing (157 TFs) | bipartite-layer | TF-FFLs | 2353 | 1379.7 (86.2) | 1.706 | +70.6 | 0.0010 | 1.0000 |
| 25 | validated_plus_predicted | published typing (157 TFs) | bipartite-layer | miRNA-FFLs | 53889 | 48210.2 (479.7) | 1.118 | +11.8 | 0.0010 | 1.0000 |
| 26 | validated | 22 co-regulators typed as genes | degree-preserving | all three-node FFLs | 63414 | 66354.0 (696.7) | 0.956 | -4.4 | 1.0000 | 0.0010 |
| 27 | validated | 22 co-regulators typed as genes | degree-preserving | composite FFLs | 7834 | 7290.7 (219.5) | 1.075 | +7.5 | 0.0050 | 0.9960 |
| 28 | validated | 22 co-regulators typed as genes | degree-preserving | TF-FFLs | 2158 | 2626.5 (148.3) | 0.822 | -17.8 | 0.9990 | 0.0020 |
| 29 | validated | 22 co-regulators typed as genes | degree-preserving | miRNA-FFLs | 42096 | 44087.4 (468.5) | 0.955 | -4.5 | 1.0000 | 0.0010 |
| 30 | validated | 22 co-regulators typed as genes | reciprocity-preserving | all three-node FFLs | 63414 | 63960.9 (621.2) | 0.991 | -0.9 | 0.7982 | 0.2028 |
| 31 | validated | 22 co-regulators typed as genes | reciprocity-preserving | composite FFLs | 7834 | 8952.6 (171.5) | 0.875 | -12.5 | 1.0000 | 0.0010 |
| 32 | validated | 22 co-regulators typed as genes | reciprocity-preserving | TF-FFLs | 2158 | 1158.3 (81.9) | 1.863 | +86.3 | 0.0010 | 1.0000 |
| 33 | validated | 22 co-regulators typed as genes | reciprocity-preserving | miRNA-FFLs | 42096 | 41781.8 (397.2) | 1.008 | +0.8 | 0.2338 | 0.7682 |
| 34 | validated | 22 co-regulators typed as genes | bipartite-layer | all three-node FFLs | 63414 | 58475.2 (613.0) | 1.084 | +8.4 | 0.0010 | 1.0000 |
| 35 | validated | 22 co-regulators typed as genes | bipartite-layer | composite FFLs | 7834 | 8411.3 (168.4) | 0.931 | -6.9 | 1.0000 | 0.0010 |
| 36 | validated | 22 co-regulators typed as genes | bipartite-layer | TF-FFLs | 2158 | 1256.2 (81.7) | 1.718 | +71.8 | 0.0010 | 1.0000 |
| 37 | validated | 22 co-regulators typed as genes | bipartite-layer | miRNA-FFLs | 42096 | 37310.2 (411.1) | 1.128 | +12.8 | 0.0010 | 1.0000 |
| 38 | validated_plus_predicted | 22 co-regulators typed as genes | degree-preserving | all three-node FFLs | 67546 | 70317.9 (708.8) | 0.961 | -3.9 | 1.0000 | 0.0010 |
| 39 | validated_plus_predicted | 22 co-regulators typed as genes | degree-preserving | composite FFLs | 8106 | 7488.7 (219.6) | 1.082 | +8.2 | 0.0020 | 0.9990 |
| 40 | validated_plus_predicted | 22 co-regulators typed as genes | degree-preserving | TF-FFLs | 2237 | 2826.1 (164.5) | 0.792 | -20.8 | 1.0000 | 0.0010 |
| 41 | validated_plus_predicted | 22 co-regulators typed as genes | degree-preserving | miRNA-FFLs | 45334 | 47197.4 (481.9) | 0.961 | -3.9 | 1.0000 | 0.0010 |
| 42 | validated_plus_predicted | 22 co-regulators typed as genes | reciprocity-preserving | all three-node FFLs | 67546 | 67782.8 (642.7) | 0.997 | -0.3 | 0.6454 | 0.3556 |
| 43 | validated_plus_predicted | 22 co-regulators typed as genes | reciprocity-preserving | composite FFLs | 8106 | 9279.8 (170.3) | 0.874 | -12.6 | 1.0000 | 0.0010 |
| 44 | validated_plus_predicted | 22 co-regulators typed as genes | reciprocity-preserving | TF-FFLs | 2237 | 1210.0 (84.7) | 1.849 | +84.9 | 0.0010 | 1.0000 |
| 45 | validated_plus_predicted | 22 co-regulators typed as genes | reciprocity-preserving | miRNA-FFLs | 45334 | 44797.3 (422.2) | 1.012 | +1.2 | 0.0939 | 0.9081 |
| 46 | validated_plus_predicted | 22 co-regulators typed as genes | bipartite-layer | all three-node FFLs | 67546 | 62241.4 (626.8) | 1.085 | +8.5 | 0.0010 | 1.0000 |
| 47 | validated_plus_predicted | 22 co-regulators typed as genes | bipartite-layer | composite FFLs | 8106 | 8701.7 (168.8) | 0.932 | -6.8 | 1.0000 | 0.0010 |
| 48 | validated_plus_predicted | 22 co-regulators typed as genes | bipartite-layer | TF-FFLs | 2237 | 1306.3 (83.9) | 1.712 | +71.2 | 0.0010 | 1.0000 |
| 49 | validated_plus_predicted | 22 co-regulators typed as genes | bipartite-layer | miRNA-FFLs | 45334 | 40340.5 (424.7) | 1.124 | +12.4 | 0.0010 | 1.0000 |
| 50 | analysed network | published typing (157 TFs) | degree-preserving | all three-node FFLs | 1774 | 2638.1 (80.6) | 0.672 | -32.8 | 1.0000 | 0.0010 |
| 51 | analysed network | published typing (157 TFs) | degree-preserving | composite FFLs | 1296 | 398.9 (34.1) | 3.249 | +224.9 | 0.0010 | 1.0000 |
| 52 | analysed network | published typing (157 TFs) | degree-preserving | TF-FFLs | 9 | 934.9 (42.6) | 0.010 | -99.0 | 1.0000 | 0.0010 |
| 53 | analysed network | published typing (157 TFs) | degree-preserving | miRNA-FFLs | 202 | 1141.9 (47.9) | 0.177 | -82.3 | 1.0000 | 0.0010 |
| 54 | analysed network | published typing (157 TFs) | reciprocity-preserving | all three-node FFLs | 1774 | 1644.3 (53.5) | 1.079 | +7.9 | 0.0130 | 0.9880 |
| 55 | analysed network | published typing (157 TFs) | reciprocity-preserving | composite FFLs | 1296 | 1217.3 (44.3) | 1.065 | +6.5 | 0.0430 | 0.9590 |
| 56 | analysed network | published typing (157 TFs) | reciprocity-preserving | TF-FFLs | 9 | 16.6 (4.8) | 0.542 | -45.8 | 0.9650 | 0.0569 |
| 57 | analysed network | published typing (157 TFs) | reciprocity-preserving | miRNA-FFLs | 202 | 232.4 (19.7) | 0.869 | -13.1 | 0.9451 | 0.0629 |
| 58 | analysed network | published typing (157 TFs) | bipartite-layer | all three-node FFLs | 1774 | 1595.1 (56.3) | 1.112 | +11.2 | 0.0010 | 1.0000 |
| 59 | analysed network | published typing (157 TFs) | bipartite-layer | composite FFLs | 1296 | 1201.6 (45.9) | 1.079 | +7.9 | 0.0270 | 0.9750 |
| 60 | analysed network | published typing (157 TFs) | bipartite-layer | TF-FFLs | 9 | 15.8 (4.8) | 0.570 | -43.0 | 0.9550 | 0.0809 |
| 61 | analysed network | published typing (157 TFs) | bipartite-layer | miRNA-FFLs | 202 | 188.6 (19.1) | 1.071 | +7.1 | 0.2607 | 0.7602 |
| 62 | analysed network | 22 co-regulators typed as genes | degree-preserving | all three-node FFLs | 1687 | 2650.0 (83.2) | 0.637 | -36.3 | 1.0000 | 0.0010 |
| 63 | analysed network | 22 co-regulators typed as genes | degree-preserving | composite FFLs | 1240 | 390.7 (35.0) | 3.174 | +217.4 | 0.0010 | 1.0000 |
| 64 | analysed network | 22 co-regulators typed as genes | degree-preserving | TF-FFLs | 9 | 882.2 (42.0) | 0.010 | -99.0 | 1.0000 | 0.0010 |
| 65 | analysed network | 22 co-regulators typed as genes | degree-preserving | miRNA-FFLs | 199 | 1085.4 (46.4) | 0.183 | -81.7 | 1.0000 | 0.0010 |
| 66 | analysed network | 22 co-regulators typed as genes | reciprocity-preserving | all three-node FFLs | 1687 | 1712.9 (55.8) | 0.985 | -1.5 | 0.6783 | 0.3337 |
| 67 | analysed network | 22 co-regulators typed as genes | reciprocity-preserving | composite FFLs | 1240 | 1179.8 (44.4) | 1.051 | +5.1 | 0.0909 | 0.9131 |
| 68 | analysed network | 22 co-regulators typed as genes | reciprocity-preserving | TF-FFLs | 9 | 15.0 (4.6) | 0.599 | -40.1 | 0.9361 | 0.1049 |
| 69 | analysed network | 22 co-regulators typed as genes | reciprocity-preserving | miRNA-FFLs | 199 | 223.6 (18.5) | 0.890 | -11.0 | 0.9121 | 0.0969 |
| 70 | analysed network | 22 co-regulators typed as genes | bipartite-layer | all three-node FFLs | 1687 | 1663.6 (55.1) | 1.014 | +1.4 | 0.3407 | 0.6633 |
| 71 | analysed network | 22 co-regulators typed as genes | bipartite-layer | composite FFLs | 1240 | 1164.3 (42.8) | 1.065 | +6.5 | 0.0380 | 0.9650 |
| 72 | analysed network | 22 co-regulators typed as genes | bipartite-layer | TF-FFLs | 9 | 14.1 (4.4) | 0.639 | -36.1 | 0.9161 | 0.1499 |
| 73 | analysed network | 22 co-regulators typed as genes | bipartite-layer | miRNA-FFLs | 199 | 180.6 (19.2) | 1.102 | +10.2 | 0.1788 | 0.8302 |

## 4. Decision (`decision.csv`; persists = p_upper < 0.05 under both the reciprocity-preserving and the
bipartite-layer null; reverses = p_lower < 0.05 under both; vanishes otherwise)

| line | network | typing | metric | reciprocity-preserving: excess % (p upper) | bipartite-layer: excess % (p upper) | decision |
|---|---|---|---|---|---|---|
| 2 | validated | published typing (157 TFs) | all three-node FFLs | -0.3 (0.5994) | +9.5 (0.0010) | vanishes |
| 3 | validated | published typing (157 TFs) | composite FFLs | -12.0 (1.0000) | -6.1 (0.9980) | reverses |
| 4 | validated | published typing (157 TFs) | TF-FFLs | +84.8 (0.0010) | +72.8 (0.0010) | persists |
| 5 | validated | published typing (157 TFs) | miRNA-FFLs | +0.3 (0.3746) | +11.9 (0.0010) | vanishes |
| 6 | validated_plus_predicted | published typing (157 TFs) | all three-node FFLs | +0.2 (0.4356) | +9.4 (0.0010) | vanishes |
| 7 | validated_plus_predicted | published typing (157 TFs) | composite FFLs | -12.2 (1.0000) | -6.3 (0.9990) | reverses |
| 8 | validated_plus_predicted | published typing (157 TFs) | TF-FFLs | +83.0 (0.0010) | +70.6 (0.0010) | persists |
| 9 | validated_plus_predicted | published typing (157 TFs) | miRNA-FFLs | +0.9 (0.1938) | +11.8 (0.0010) | vanishes |
| 10 | validated | 22 co-regulators typed as genes | all three-node FFLs | -0.9 (0.7982) | +8.4 (0.0010) | vanishes |
| 11 | validated | 22 co-regulators typed as genes | composite FFLs | -12.5 (1.0000) | -6.9 (1.0000) | reverses |
| 12 | validated | 22 co-regulators typed as genes | TF-FFLs | +86.3 (0.0010) | +71.8 (0.0010) | persists |
| 13 | validated | 22 co-regulators typed as genes | miRNA-FFLs | +0.8 (0.2338) | +12.8 (0.0010) | vanishes |
| 14 | validated_plus_predicted | 22 co-regulators typed as genes | all three-node FFLs | -0.3 (0.6454) | +8.5 (0.0010) | vanishes |
| 15 | validated_plus_predicted | 22 co-regulators typed as genes | composite FFLs | -12.6 (1.0000) | -6.8 (1.0000) | reverses |
| 16 | validated_plus_predicted | 22 co-regulators typed as genes | TF-FFLs | +84.9 (0.0010) | +71.2 (0.0010) | persists |
| 17 | validated_plus_predicted | 22 co-regulators typed as genes | miRNA-FFLs | +1.2 (0.0939) | +12.4 (0.0010) | vanishes |
| 18 | analysed network | published typing (157 TFs) | all three-node FFLs | +7.9 (0.0130) | +11.2 (0.0010) | persists (reference) |
| 19 | analysed network | published typing (157 TFs) | composite FFLs | +6.5 (0.0430) | +7.9 (0.0270) | persists (reference) |
| 20 | analysed network | published typing (157 TFs) | TF-FFLs | -45.8 (0.9650) | -43.0 (0.9550) | vanishes (reference) |
| 21 | analysed network | published typing (157 TFs) | miRNA-FFLs | -13.1 (0.9451) | +7.1 (0.2607) | vanishes (reference) |
| 22 | analysed network | 22 co-regulators typed as genes | all three-node FFLs | -1.5 (0.6783) | +1.4 (0.3407) | vanishes (reference) |
| 23 | analysed network | 22 co-regulators typed as genes | composite FFLs | +5.1 (0.0909) | +6.5 (0.0380) | vanishes (reference) |
| 24 | analysed network | 22 co-regulators typed as genes | TF-FFLs | -40.1 (0.9361) | -36.1 (0.9161) | vanishes (reference) |
| 25 | analysed network | 22 co-regulators typed as genes | miRNA-FFLs | -11.0 (0.9121) | +10.2 (0.1788) | vanishes (reference) |

