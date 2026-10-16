# SETTINGS: tier- and sign-aware null models for the three-node FFL test

The design of the analysis.

## Graphs (public repository; FFL_ROOT = a checkout of it)

| name | file | arcs |
|---|---|---|
| analysed network | `analyses/census_and_motif_nulls/graphs/null_dep_noexemplar.txt` | 6,829 |
| census graph | `analyses/census_and_motif_nulls/graphs/null_pub_noexemplar.txt` | 9,226 |

- Both files are uncontracted, with per-arc fine classes and node names in the `.meta.json`.
- Layer labels of the census-graph arcs come from `analyses/six_node_pattern/S1/graph_noexemplar_labels.txt`, which uses
  the same arc order (the md5 of the S1 graph equals that of `null_pub_noexemplar.txt`).

## Typings

- **Published.** 157 TFs.
- **Retyped.** The 22 co-regulators of `analyses/six_node_pattern/E/E5/e5bc_retype_graphs.py` (APEX1 … VHL) get type code
  1 -> 2 (gene). The fine classes (source type -> target type, sorted names) are recomputed, and the arc order is
  unchanged.
  - For the analysed network this must reproduce `E5/null_dep_noexemplar_retyped22.txt` byte for byte.
  - The census graph is retyped the same way.

## Arc attributes

- **Evidence tier**, for miRNA -> {Gene, TF} arcs: `data/network/edge_evidence_tier.tsv`, column tier (strong / weak /
  predicted_only). All 4,819 miRNA_target edges carry one.
- **Annotated sign**, for TF -> {TF, Gene} arcs whose source is TF-typed in the typing used:
  - analysed TF_target arcs: Table S2 trrust_mode;
  - TRRUST-layer arcs (label 3): `data/network/layer_TF_target.tsv`, column mode.
  - Activation = "+", Repression = "−", anything else (Unknown, Ambiguous, blank) = "unannotated".
  - STRING arcs leaving a TF (label 5) are "unannotated".

## Null models

Every model is degree-preserving and, as in the paper, starts every replicate from the observed graph. Shared settings:
- 100 swap attempts (or curveball trades) per arc of each group;
- xorshift64 seeded 20250908;
- 1,000 replicates;
- counting with `count3()` of `analyses/census_and_motif_nulls/nulls/motif_null_models.c`, verbatim. After randomisation,
  each reciprocal TF<->miRNA pair is contracted to its TF->miRNA arc and the induced transitive triangles are counted:
  - all three-node FFLs;
  - composite FFLs, counted once per reciprocal pair;
  - TF-FFLs;
  - miRNA-FFLs.

The models:
- **The degree-preserving null, the reciprocity-preserving null, the bipartite-layer null.** Verbatim from `motif_null_models.c`. With the published typing they must reproduce the
  published replicate files `nulls/null_{A,B,C}_{dep,pub}_noexemplar.tsv` line for line, and with the retyped analysed
  network `E5/null_{A,B,C}_dep_noexemplar_retyped22.tsv`. Nothing else is reported until these controls pass.
- **(i) the reciprocity-preserving null + tier.** the reciprocity-preserving null freezing (both arcs of every reciprocal TF<->miRNA pair). The free miRNA -> {Gene, TF}
  arcs are randomised by curveball trades as in the bipartite-layer null, but as three separate bipartite groups, one per tier. Every other
  class is swapped within class as in the reciprocity-preserving null.
- **(ii) the reciprocity-preserving null + sign.** the reciprocity-preserving null, except that the TF -> Gene and TF -> TF classes are each split by annotated sign (+, −,
  unannotated). The Maslov–Sneppen swaps then run within class × sign.
- **(iii) the reciprocity-preserving null + tier + sign.** (i) and (ii) together.
- **Waterfall steps beyond the bipartite-layer null:**
  - **The bipartite-layer null + tier:** the bipartite-layer null, with its miRNA-side curveball group split into the three tiers.
  - **The bipartite-layer null + tier + sign:** additionally, its TF-side curveball group (TF -> {TF, Gene}) is split into the three sign
    groups.
- **Group order** (it fixes the random stream): the fine-class groups in class order, then the split groups ordered by
  (layer side miRNA before TF, then tier strong / weak / predicted_only, or sign + / − / unannotated).

## Runs

8 models (A, B, C, i, ii, iii, C + tier, C + tier + sign) × 2 graphs × 2 typings = 32 runs of 1,000 replicates.

## Reported

For each run and each metric (all, composite, TF-FFL, miRNA-FFL):
- observed; expected (the mean of the 1,000) with its SD;
- fold (observed / expected); excess % (100 × (fold − 1));
- empirical p_upper = (1 + #{random ≥ observed}) / 1,001 and p_lower = (1 + #{random ≤ observed}) / 1,001.

**Waterfall table:** for each graph × typing × metric, the excess % at the degree-preserving null -> the reciprocity-preserving null (+ reciprocity) -> the bipartite-layer null
(+ bipartite curveball) -> the bipartite-layer null + tier -> the bipartite-layer null + tier + sign, with the step-to-step change in percentage points.

No decision rule is attached; the numbers are reported whichever way they fall.

## Software

C (`gcc -O2`, as in the published run) and Python 3 (numpy) for the summaries.

