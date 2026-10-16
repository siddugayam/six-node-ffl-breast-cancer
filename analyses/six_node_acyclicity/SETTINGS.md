# SETTINGS: the six-node pattern and acyclicity

The design of the analysis.

Paper numbering (Fig. 1):
- D1: induced and weakly connected;
- D2: acyclic;
- D3: one source and one sink, every vertex on a source–sink path;
- D4: two internally vertex-disjoint source–sink paths.

## (a) The six-node composite pattern against the topological definition

**Pattern.** 3,898 instances on 2,498 node sets. They are counted by `bhat6()` of
`analyses/six_node_pattern/S1/s1_null6.c` on the census graph, and the deposited instance list must be reproduced (OBS
mode).

**First failing condition.** For each node set, the first of the topological definition that fails is determined on the census graph
(contracted, as the census uses).

**Failure classes.**
- The kind of directed cycle that makes D2 fail:
  - a co-transcribed (10 kb) miRNA pair, entered as two opposed arcs;
  - a reciprocal TF<->TF pair;
  - a reciprocal STRING/TRRUST gene pair;
  - a cycle of length ≥ 3;
  - and their combinations.
- For each class, it is also reported whether D3 and D4 would hold if D2 were waived:
  - literally, on the cyclic module;
  - with each strongly connected component of the module condensed.

## (b) D2′: one directed cycle allowed among {TF1, TF2, miRNA}

### The module

A six-node module that contains the composite core plus a second TF, the CENSUS6c seed of `s1_null6.c`:
- TF1 -> miRNA and miRNA -> TF1 (the reciprocal pair, both arcs present in the uncontracted census graph);
- TF1 -> TF2, TF2 -> miRNA;
- TF1 -> t, TF2 -> t and miRNA -> t, where t is a non-miRNA node;
- TF1 and TF2 are TF-typed;
- plus two further vertices, so that the six-set is weakly connected.

It is enumerated exactly as CENSUS6c does: every seed is extended by every pair of further vertices, and each six-set is
counted once.

### Conditions

- D1, D3 and D4 are unchanged.
- **D2′.** The induced subgraph has exactly one strongly connected component with more than one vertex, and that
  component lies within {TF1, TF2, miRNA} for at least one seed (TF1, TF2, miRNA, t) contained in the six-set.
  - This reads "exactly one directed cycle" as one cyclic strongly connected component: every seed contains two
    directed cycles, TF1<->miRNA and TF1->TF2->miRNA->TF1.
- **Primary.** D3 and D4 are tested on the induced subgraph as it is, with its cycle. Because TF1, TF2 and the miRNA
  all have an in-arc inside the cycle, the source must be one of the other three vertices.
- **Secondary.** D3 and D4 are tested after the cycle is condensed to one vertex (the census convention).

Both variants are counted on the observed census graph and on every null replicate.

### Null models

- the bipartite-layer null and the core-conditional null of `s1_null6.c`, unchanged: same xorshift64 stream, seed 20250908, 100 swaps or trades per
  rewired item, and each replicate restarting from the observed graph.
- The six-node counting consumes no random numbers. The replicates are therefore the same graphs as in the S1 runs,
  which is checked replicate by replicate on S1's stored columns (total3 … census6c) for the replicates of the stored
  runs (1,000 per model).

### Replicates

- **Rule.** 5,000 per model if the projected wall time is under 4 h; otherwise 1,000.
- **Projection.** The time per replicate of the counting is measured on the observed graph. It is multiplied by
  5,000 / (processes per model), plus the randomisation time of the S1 runs, with two models run at once on the
  processes available.
- With 5,000 replicates, the first 1,000 are the S1 replicates.

### Test and decision rule (as S1)

- p_upper = (1 + #{null ≥ observed}) / (R + 1); p_lower likewise.
- **Over-represented** if p_upper < 0.05 under both the bipartite-layer null and the core-conditional null.
- **Depleted** if p_lower < 0.05 under both.
- Otherwise "none of the three".
- Primary and secondary are reported separately. The rule is applied to the primary count, for instances (six-sets).

### Also reported

The topological-definition six-node count with the same seed (CENSUS6c, contracted graph) on the same replicates, for comparison.

## Outputs

- `pattern_first_failing_condition.csv`
- `d2prime_counts_observed.csv`
- `d2prime_null_summary.csv`
- `d2prime_decision.csv`
- `runs/`: per-replicate counts.
- `check_replicates_against_s1.txt`

## Replicates

- Counting the observed graph takes about 13 s, so 5,000 replicates per model would take about 6 h with four
  processes; the analysis uses **1,000 replicates per model**, the number used for the other null models.
