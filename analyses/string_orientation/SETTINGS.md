# SETTINGS: orientation of the STRING associations in the census graph

Settings of the analysis.

## Background

- STRING associations are undirected. The census graph enters each one as a single arc in the orientation listed in
  `data/network/layer_gene_gene.tsv` (source < target alphabetically).
- That file lists 941 STRING associations (score ≥ 0.900), all with both nodes in the network.
- The census builder (`build_graph` of `scripts/03_ffl_census/03_ffl_census_noexemplarmirna.py`) adds an arc only if it
  is absent. Its order is:
  - the analysed network (6,829 arcs; the 30 exemplar miRNA–miRNA edges left out);
  - then the TRRUST layer (`layer_TF_target.tsv`, 604 arcs added);
  - then the gene–gene file: 833 STRING arcs added. The other 108 STRING arcs already exist in the listed orientation
    (91 TRRUST layer, 17 analysed network). The file's 697 TRRUST rows add none;
  - then the 10 kb miRNA pairs (960 arcs, both directions);
  - and finally it replaces each of the 1,223 reciprocal TF<->miRNA pairs by its TF->miRNA arc (8,003 arcs).
- The four conditions of the topological definition are those of the paper (Fig. 1). The census code is the verbatim enumerator and test of
  `analyses/census_and_motif_nulls/census/ffl_census_composition.c` (exhaustive at n = 3–6; RAND-ESU at n = 7 with the
  published 8 seeds).

## Graph construction (one script, `string_orientation_graphs.py`)

The builder above is re-implemented from the public files, and then checked:
- with every STRING association in its listed orientation, it must reproduce `census_noexemplar_fine.txt` arc for arc;
- without STRING, `census_noexemplar_nostring` (7,170 arcs).
Nothing is counted unless both checks pass.

## (a) Orientation bootstrap

- **Replicates.** In each, every one of the 941 associations is oriented a->b or b->a with probability 1/2. The builder
  is then applied unchanged, so an arc is added only if absent, and the reciprocal TF<->miRNA contraction runs last as in
  the census.
  - A flipped association can therefore add an arc where the listed one added none, and the reverse.
- **Random numbers.** numpy `default_rng(20261015 + replicate)`, one Bernoulli(0.5) vector of length 941 per replicate,
  in the file's row order.
- **Size.** 500 replicates; each is counted exhaustively at n = 3, 4, 5 and 6.
- **Reported per n:** the median count, the 2.5 % and 97.5 % percentiles, the mean, the SD, the coefficient of variation
  (SD / mean), and where the published count lies in the bootstrap distribution (its percentile).
- **Also per replicate:** the number of STRING arcs that entered the graph.

## (b) Gene–gene layer variants (exhaustive, n = 3–6, and the n = 7 RAND-ESU estimate)

1. **TRRUST-only gene–gene links:** the census graph minus every STRING arc. This is the published "noexemplar_nostring"
   graph, with 3,432 / 81,104 / 791,085 / 8,751,860 modules, which must be reproduced.
2. **No gene–gene layer:** the census graph minus the STRING arcs and minus the 604 TRRUST-layer arcs, all of which lie
   between protein-coding nodes. That leaves the analysed network plus the 10 kb layer, contracted.

## (c) The orientation-robust core

For the baseline census instances at each n = 3–6 (exact, from the instance list) and at n = 7 (from the 8 RAND-ESU
samples), two numbers and fractions:
1. **No STRING arc:** instances whose induced subgraph contains no STRING arc, i.e. no arc that the gene–gene file added.
2. **No STRING pair:** instances with no pair of their nodes forming one of the 941 STRING associations. Only these are
   unchanged under every orientation, because a flip can add an arc to a node set that held none.

Also reported for the n = 3–6 bootstrap: per replicate, the share of the baseline instances found again (exact, by
instance keys). Its median and 95 % interval are given.

## Seeds and software

- RAND-ESU at n = 7: xorshift64 with the published seeds and keep-probabilities.
  - Seeds 11, 22, 33, 44, 55: q = 1, 1, 0.2, 0.2, 0.2, 0.2.
  - Seeds 71, 72, 73: q = 1, 1, 1, 0.1, 0.15, 0.15.
- Bootstrap seeds as above.
- C (`gcc -O3`) for counting; Python 3 for the graphs and the summaries.

## Outputs

- `string_bootstrap_counts.csv`: one row per replicate and n.
- `string_bootstrap_summary.csv`
- `gene_gene_layer_variants.csv`
- `orientation_robust_core.csv`

