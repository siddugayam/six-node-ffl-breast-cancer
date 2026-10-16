# SETTINGS: reciprocal handling in the census

Settings of the analysis.

## Question

How do the topological-definition census counts (n = 3–7) and the six-node composite pattern change when the reciprocal TF<->miRNA pairs
are handled differently? The paper's census replaces each of the 1,223 reciprocal pairs by its transcriptional arc
TF -> miRNA ("contracted"). This analysis compares that baseline with three alternatives.

## Definitions (paper numbering, Fig. 1)

- **D1.** Induced subgraph on n vertices, weakly connected.
- **D2.** Acyclic.
- **D3.** Exactly one source and one sink; every vertex lies on a source–sink path.
- **D4.** At least two internally vertex-disjoint source–sink paths.

The test code is the verbatim `acyclic()` + `check_module()` of
`analyses/census_and_motif_nulls/census/ffl_census_composition.c` (its own comments number the conditions differently).

## Graphs (all from the public repository; FFL_ROOT = a checkout of it)

| name | file | arcs |
|---|---|---|
| census graph, contracted (baseline) | `analyses/census_and_motif_nulls/graphs/census_noexemplar_fine.txt` | 8,003 |
| census graph, uncontracted | `analyses/census_and_motif_nulls/graphs/census_noexemplar_uncontracted.txt` | 9,226 |

The two files list the same 587 nodes in the same order.

## Variants

- **B, baseline.** the topological definition on the contracted census graph. This must reproduce the published counts:
  - exhaustive: 5,833 / 116,505 / 1,506,263 / 19,401,840 at n = 3–6;
  - n = 7: the published RAND-ESU estimates per seed.
  - If it does not reproduce them, nothing else is run.
- **U, no contraction.** the topological definition on the uncontracted census graph. Both arcs of every reciprocal pair are kept, so D2
  rejects any node set that holds both members of a pair.
  - Enumeration: the same enumerator with the uncontracted graph as its structure.
- **S, graph-level SCC condensation.** Every strongly connected component (SCC) of the uncontracted census graph is
  contracted to one node, giving the condensation DAG. the topological definition are then counted on that DAG, which is acyclic, so D2 always
  holds.
  - **Instance identity.** A baseline instance maps to the set of SCCs of its nodes (its image).
    - It is **retained** if the image has n distinct SCCs and is a topological-definition module of the condensation DAG.
    - Otherwise it is **lost**.
  - **Shared set:** the distinct images of retained baseline instances.
  - **Gained:** condensation modules that are not the image of any baseline instance.
  - **Jaccard** is computed on images: |shared| / (|distinct images of baseline instances with n distinct SCCs| +
    |condensation modules| − |shared|).
- **M, module-level condensation**. Take a node set S on n vertices whose induced
  subgraph in the uncontracted census graph satisfies D1.
  - Contract each SCC of that induced subgraph (not of the whole graph) to one node. The result is acyclic.
  - Apply D3 and D4 to it. n still counts the vertices before contraction.
  - **Enumeration.** Each valid S has one source component, from which every vertex of S is reachable. So S is
    enumerated from a root r: the smallest-index vertex of that component.
    - It is grown along out-arcs within S, with no acyclicity pruning, in the same ESU order.
    - It is counted only if r is that smallest vertex, so that each S is counted once.
  - **Instance identity** is the node set itself. Jaccard = |B ∩ M| / |B ∪ M|.
  - Also reported: the number of composite three-node FFLs that fail D4 under M, and the first failing condition for each
    baseline instance lost under M.

For B, U and M the instance is the node set, and Jaccard = |X ∩ B| / |X ∪ B|. The quantities reported are:
- **lost fraction** = |B \ X| / |B|;
- **gained fraction** = |X \ B| / |X|.

Instance sets are compared exactly: every n-set is packed into one 64-bit key (sorted 10-bit node indices; n ≤ 6), and
the key files are sorted and intersected.

## n = 7

- **B, U and M:** RAND-ESU with the published 8 seeds and keep-probabilities.
  - Seeds 11, 22, 33, 44, 55 use q = 1, 1, 0.2, 0.2, 0.2, 0.2.
  - Seeds 71, 72, 73 use q = 1, 1, 1, 0.1, 0.15, 0.15.
- **Estimates:** the count of each seed is divided by the product of its q; the estimate reported is the mean over the 8
  seeds, with their SD.
- **Overlap with the baseline at n = 7** is estimated from the baseline samples:
  - the fraction of sampled baseline instances that are valid under U, retained under S, or valid under M;
  - for M, also the fraction of sampled M instances that are not baseline instances.
- **The Jaccard estimate at n = 7:** est |X ∩ B| / (est |B| + est |X| − est |X ∩ B|), with est |X ∩ B| = est |B| × the
  retained fraction.
- **S at n = 7** is exhaustive, because the condensation DAG is small.

## Contingency for M at n = 6 (fixed now)

M enumerates without acyclicity pruning, so its cost is unknown until it runs.
- Timing is taken on n = 4 and n = 5 first. Then:
  - if the projected n = 6 CPU time is ≤ 12 CPU-hours, n = 6 is exhaustive (64 partitions);
  - otherwise n = 6 uses RAND-ESU with the n = 7 scheme of seeds 11–55, shortened to six nodes (q = 1, 1, 0.2, 0.2, 0.2),
    and this is stated.
- Comparing B with M on the baseline instances is exact at every n regardless, because the M test is applied to each
  baseline instance directly.

## Six-node composite pattern (3,898 instances on 2,498 node sets)

- The pattern is counted by `bhat6()` of `analyses/six_node_pattern/S1/s1_null6.c` (mode OBS) on the uncontracted census
  graph. It is recounted here with that program, unchanged, to confirm 3,898 / 2,498.
- **Under U.** The pattern reads both arc directions (TF<->miRNA as a link), so its graph is the same as in U and the count
  cannot change. This is confirmed by the recount.
- **Under S.** The pattern needs six distinct nodes in its roles. Its instances are mapped to SCC images, and the
  instances whose six nodes fall in six distinct SCCs are counted.
- **Under M.** The pattern is a literal arc pattern, not a topological-definition module. Reported: the number of its 2,498 node sets that
  are M-valid six-node modules, and the number that are baseline (B) modules.

## Seeds and software

- No random numbers are used except the RAND-ESU draws at n = 7: xorshift64 with the seeds above, in the published
  order of draws.
- C compiled with `gcc -O3 -march=native`. Python 3 for graph building (networkx for the SCCs) and the key comparisons
  (numpy).

## Outputs

- `census_reciprocal_variants.csv`, one row per n × variant, with columns:
  - n, variant, method, baseline_count, variant_count, shared, jaccard, lost_fraction, gained_fraction, the n = 7 SD;
- `six_node_pattern_variants.csv`;
- `composite_three_node_under_module_condensation.csv`;
- `lost_under_module_condensation_first_condition.csv`;
- `scc_structure.csv`: SCC sizes and the node types of the largest SCC.

