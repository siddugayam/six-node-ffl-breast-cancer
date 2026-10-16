# VALUES: reciprocal TF<->miRNA handling in the topological-definition census

Lines count the header as line 1. Tables are written by `summarise_reciprocal_variants.py`,
`six_node_pattern_variants.py` and `make_scc_file.py` from the counts in `counts/` (run by
`run_reciprocal_variants.sh`). The baseline is the published census (contracted census graph, 8,003 arcs).

## 1. Strongly connected components of the census graph (`scc_structure.csv`)

| line | graph | component (smallest vertex) | size | TF | miRNA | Gene | note |
|---|---|---|---|---|---|---|---|
| 2 | uncontracted | 0 | 452 | 157 | 219 | 76 |  |
| 3 | uncontracted | all | 587 | 157 | 223 | 207 | 9226 arcs; 136 SCCs (1 non-singleton); condensation DAG 136 vertices, 188 arcs |
| 4 | contracted | 0 | 435 | 149 | 210 | 76 |  |
| 5 | contracted | all | 587 | 157 | 223 | 207 | 8003 arcs; 153 SCCs (1 non-singleton); condensation DAG 153 vertices, 250 arcs |

## 2. Counts, overlap with the baseline, lost and gained (`census_reciprocal_variants.csv`)

Columns baseline_count, variant_count, shared, jaccard, lost_fraction (of the baseline), gained_fraction (of the
variant), variant_sd and baseline_sd (n = 7, SD of the 8 seed estimates).

| line | n | variant | method | baseline | variant | shared | Jaccard | lost | gained | SD (variant / baseline) |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 3 | U no contraction | exhaustive | 5,833 | 2,914 | 2,914 | 0.4996 | 0.5004 | 0.0000 |  |
| 3 | 3 | S graph-level SCC condensation | exhaustive (images) | 5,833 | 53 | 34 | 0.6415 | 0.9820 | 0.3585 |  |
| 4 | 3 | M module-level condensation | exhaustive | 5,833 | 2,914 | 2,914 | 0.4996 | 0.5004 | 0.0000 |  |
| 5 | 4 | U no contraction | exhaustive | 116,505 | 23,845 | 23,845 | 0.2047 | 0.7953 | 0.0000 |  |
| 6 | 4 | S graph-level SCC condensation | exhaustive (images) | 116,505 | 38 | 24 | 0.6316 | 0.9993 | 0.3684 |  |
| 7 | 4 | M module-level condensation | exhaustive | 116,505 | 164,199 | 56,889 | 0.2542 | 0.5117 | 0.6535 |  |
| 8 | 5 | U no contraction | exhaustive | 1,506,263 | 187,719 | 187,719 | 0.1246 | 0.8754 | 0.0000 |  |
| 9 | 5 | S graph-level SCC condensation | exhaustive (images) | 1,506,263 | 36 | 17 | 0.4722 | 1.0000 | 0.5278 |  |
| 10 | 5 | M module-level condensation | exhaustive | 1,506,263 | 7,255,939 | 1,039,673 | 0.1346 | 0.3098 | 0.8567 |  |
| 11 | 6 | U no contraction | exhaustive | 19,401,840 | 1,335,762 | 1,335,762 | 0.0688 | 0.9312 | 0.0000 |  |
| 12 | 6 | S graph-level SCC condensation | exhaustive (images) | 19,401,840 | 31 | 6 | 0.1935 | 1.0000 | 0.8065 |  |
| 13 | 6 | M module-level condensation | exhaustive | 19,401,840 | 271,782,976 | 13,982,004 | 0.0504 | 0.2793 | 0.9486 |  |
| 14 | 7 | U no contraction | RAND-ESU, 8 seeds (S: exhaustive on the condensation DAG) | 225,761,602 | 8,760,903 | 8,730,628 | 0.0387 | 0.9613 | 0.0000 | 161262.2 / 939186.5 |
| 15 | 7 | S graph-level SCC condensation | RAND-ESU, 8 seeds (S: exhaustive on the condensation DAG) | 225,761,602 | 19 | 0 |  | 1.0000 |  | 0.0 / 939186.5 |
| 16 | 7 | M module-level condensation | RAND-ESU, 8 seeds (S: exhaustive on the condensation DAG) | 225,761,602 | 9,019,061,186 | 162,850,228 | 0.0179 | 0.2787 | 0.9819 | 22657054.3 / 939186.5 |

Notes on the rows:

- line 3: baseline modules retained: 105; distinct images of baseline modules with n distinct SCCs: 34
- line 6: baseline modules retained: 83; distinct images of baseline modules with n distinct SCCs: 24
- line 9: baseline modules retained: 40; distinct images of baseline modules with n distinct SCCs: 17
- line 12: baseline modules retained: 25; distinct images of baseline modules with n distinct SCCs: 6
- line 14: shared = baseline estimate x fraction of sampled baseline modules valid under the variant
- line 15: Jaccard on images not estimable from the samples; shared = baseline estimate x retained fraction
- line 16: shared = baseline estimate x fraction of sampled baseline modules valid under the variant

## 3. Baseline modules lost, by first failing condition (`lost_under_module_condensation_first_condition.csv`)

| line | n | baseline | lost under M | of which D3 fails first | D4 fails first | lost under U (all fail D2) | lost under S: two nodes share an SCC | lost under S: image not a module |
|---|---|---|---|---|---|---|---|---|
| 2 | 3 | 5,833 | 2,919 | 1,544 | 1,375 | 2,919 | 5,728 | 0 |
| 3 | 4 | 116,505 | 59,616 | 12,108 | 47,508 | 92,660 | 116,422 | 0 |
| 4 | 5 | 1,506,263 | 466,590 | 101,344 | 365,246 | 1,318,544 | 1,506,223 | 0 |
| 5 | 6 | 19,401,840 | 5,419,836 | 836,339 | 4,583,497 | 18,066,078 | 19,401,815 | 0 |

## 4. Composite three-node FFLs under module-level condensation (`composite_three_node_under_module_condensation.csv`)

- line 2: composite three-node FFLs 1375; valid under M 0; fail D4 1375; fail D3 0.

## 5. The six-node composite pattern (`six_node_pattern_variants.csv`, columns instances and node_sets)

| line | treatment | instances | node sets | note |
|---|---|---|---|---|
| 2 | pattern recount (s1_null6 OBS, uncontracted census graph) | 3898 | 2498 | instance list identical to the deposited one: True |
| 3 | U no contraction | 3898 | 2498 | the pattern already reads both arcs of every reciprocal pair, so its graph is the uncontracted one |
| 4 | S graph-level SCC condensation: node sets with six distinct SCCs |  | 0 | of which topological-definition modules of the condensation DAG: 0 |
| 5 | M module-level condensation: node sets that are valid modules | 3161 | 2040 |  |
| 6 | B baseline census (contracted): node sets that are topological-definition modules |  | 0 | first failing condition: D2 2498 |
| 7 | U: node sets that are topological-definition modules of the uncontracted graph |  | 0 |  |

## 6. Independent check at small n (`check_small_n_output.txt`)

A separate Python implementation (undirected ESU over all weakly connected subgraphs; its own D2–D4 test) gives the
same counts as the C enumerator:
- `n=3 B=5833 U=2914 M=2914 BU=2914 BM=2914 S=53`
- `n=4 B=116505 U=23845 M=164199 BU=23845 BM=56889 S=38`

