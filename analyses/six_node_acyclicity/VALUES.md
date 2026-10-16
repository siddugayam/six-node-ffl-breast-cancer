# VALUES: the six-node pattern and acyclicity

Lines count the header as line 1. Part (a) is written by `pattern_first_failing.py`; part (b) by `summarise_d2prime.py`
from the counts of `six_node_d2prime_nulls.c` (run by `run_d2prime_nulls.sh`).

## (a) The 2,498 node sets of the six-node composite pattern against the topological definition (`pattern_first_failing_condition.csv`)

Graph: the census graph, contracted as in the census. Columns first_failing_condition, cycle_class,
with_D2_waived_literal (D3–D4 on the module as it is), with_D2_waived_condensed (D3–D4 with its cycles condensed),
node_sets.

| line | first failing | cycles present | D2 waived: literal | D2 waived: condensed | node sets |
|---|---|---|---|---|---|
| 2 | D2 | 2-cycles: TF-TF, miRNA-miRNA | D3 | pass | 2060 |
| 3 | D2 | 2-cycles: Gene-TF, TF-TF, miRNA-miRNA | D3 | pass | 187 |
| 4 | D2 | 2-cycles: Gene-TF, TF-TF, miRNA-miRNA | D3 | D3 | 82 |
| 5 | D2 | 2-cycles: Gene-TF, TF-TF, miRNA-miRNA | D3 | D4 | 50 |
| 6 | D2 | 2-cycles: TF-TF, miRNA-miRNA; a cycle of length >= 3 | D3 | D3 | 46 |
| 7 | D2 | 2-cycles: Gene-TF, TF-TF, miRNA-miRNA; a cycle of length >= 3 | D3 | D3 | 38 |
| 8 | D2 | 2-cycles: TF-TF, miRNA-miRNA; a cycle of length >= 3 | D3 | D4 | 29 |
| 9 | D2 | 2-cycles: TF-TF, miRNA-miRNA | D3 | D4 | 6 |

Total node sets: 2498. Every set holds both a co-transcribed miRNA pair and a reciprocal TF<->TF pair; the
composite pattern requires the TF1<->TF2 link and the 10 kb miRNA pair.

## (b) D2′ six-node modules (composite core plus a second TF)

### Observed (`d2prime_counts_observed.csv`, column value)

- line 2: D2' seeds (composite core + second TF, uncontracted census graph): 5,470
- line 3: D2' six-node modules, primary (D3-D4 on the module with its cycle): 67,229
- line 4: D2' six-node modules, secondary (D3-D4 with the cycle condensed): 57,217
- line 5: CENSUS6c seeds (contracted census graph, s1_null6.c): 1,339
- line 6: CENSUS6c six-node topological-definition modules with the composite core plus TF2: 210,115
- line 7: six-node composite pattern instances (check): 3,898

### Against the bipartite-layer null and the core-conditional null (`d2prime_null_summary.csv`)

Columns observed, null_mean, null_sd, fold, p_upper, p_lower, replicates.

| line | null | count | observed | null mean (SD) | fold | p upper | p lower | replicates |
|---|---|---|---|---|---|---|---|---|
| 2 | Bipartite-layer null | D2' seeds | 5,470 | 3,124.1 (341.0) | 1.751 | 0.0010 | 1.0000 | 1000 |
| 3 | Bipartite-layer null | D2' primary | 67,229 | 35,211.8 (7,683.0) | 1.909 | 0.0010 | 1.0000 | 1000 |
| 4 | Bipartite-layer null | D2' secondary | 57,217 | 42,389.9 (8,830.8) | 1.350 | 0.0649 | 0.9361 | 1000 |
| 5 | Bipartite-layer null | CENSUS6c (the topological definition, contracted; comparison) | 210,115 | 152,780.7 (38,493.7) | 1.375 | 0.0779 | 0.9231 | 1000 |
| 6 | Core-conditional null | D2' seeds | 5,470 | 3,462.8 (315.5) | 1.580 | 0.0010 | 1.0000 | 1000 |
| 7 | Core-conditional null | D2' primary | 67,229 | 67,649.1 (5,802.8) | 0.994 | 0.5185 | 0.4825 | 1000 |
| 8 | Core-conditional null | D2' secondary | 57,217 | 63,279.5 (5,103.4) | 0.904 | 0.8831 | 0.1179 | 1000 |
| 9 | Core-conditional null | CENSUS6c (the topological definition, contracted; comparison) | 210,115 | 252,463.2 (39,469.5) | 0.832 | 0.8462 | 0.1548 | 1000 |

### Decision (`d2prime_decision.csv`, column decision; rule: p_upper < 0.05 under both the bipartite-layer null and the core-conditional null)

- line 2: D2' primary: none of the three (not over-represented, not depleted) (the bipartite-layer null fold 1.909, p_upper 0.0010; the core-conditional null fold 0.994, p_upper 0.5185); rule applies: yes
- line 3: D2' secondary: none of the three (not over-represented, not depleted) (the bipartite-layer null fold 1.350, p_upper 0.0649; the core-conditional null fold 0.904, p_upper 0.8831); rule applies: reported
- line 4: CENSUS6c (the topological definition, contracted; comparison): none of the three (not over-represented, not depleted) (the bipartite-layer null fold 1.375, p_upper 0.0779; the core-conditional null fold 0.832, p_upper 0.8462); rule applies: reported

### Replicate check (`check_replicates_against_s1.txt`)

- the bipartite-layer null: 1000 of 1000 stored S1 replicate rows identical (all columns except the timing)
- the bipartite-layer null: replicates 1000; check passed
- the core-conditional null: 1000 of 1000 stored S1 replicate rows identical (all columns except the timing)
- the core-conditional null: replicates 1000; check passed

The number of replicates (1,000 per model) follows the rule in SETTINGS.md (section "Run decisions").

