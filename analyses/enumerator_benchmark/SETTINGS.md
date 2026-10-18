# SETTINGS: benchmark of the topological-definition enumerator on planted modules

The design of the analysis.

## Enumerator under test

The census enumerator of `analyses/census_and_motif_nulls/census/ffl_census_composition.c`: source-rooted ESU,
acyclicity pruning, and `check_module()`. It is copied with key output and with two relaxed tests switched by an
argument. Exhaustive at n = 3–6.

## Background graphs (200 replicates, numpy seed 20261015 + replicate)

1. **Skeleton.** The analysed network with each reciprocal TF<->miRNA pair contracted to its TF -> miRNA arc (the graph
   the census starts from, before the census layers), taken as an undirected simple graph. The degree of every vertex is
   kept.
2. **Rewiring.** The skeleton is rewired with degree-preserving double-edge swaps: 100 attempts per edge, no self-loops,
   no multi-edges.
3. **Orientation.** A uniformly random total order of the 587 vertices orients every edge from the vertex placed first to the
   vertex placed second, so the background is a DAG with the network's vertex count, edge count and degree sequence.

## Planted modules

- **Topologies.** For each n = 3–6, the list of all topological-definition module topologies on n vertices up to isomorphism: every DAG
  on n vertices whose arcs follow one fixed vertex order, kept if it satisfies the topological definition and deduplicated by canonical form
  over all vertex permutations.
- **Placement.** In each replicate, 10 modules per n (40 in all) are planted on disjoint sets of vertices drawn at
  random.
  - Each topology is drawn uniformly from the list for its n.
  - Its vertex order is mapped onto the drawn vertices in the background order, so the graph stays acyclic.
  - Background arcs among the drawn vertices are removed first, so the planted subgraph is induced exactly.

## Noise

After planting, noise arcs equal to 10 % of the background edge count are added between random vertex pairs not
already linked, oriented by the background order (the graph stays a DAG). Noise arcs can fall inside planted modules.

## Measures, per replicate and n

- **Intact planted modules.** The independent topological-definition test (Python/networkx) is applied to each planted module's induced
  subgraph in the final graph. "Destroyed" means noise made it fail.
- **Recall of planted modules.** The share of intact planted modules whose vertex set the enumerator reports; and,
  separately, the share of destroyed ones it reports (should be 0).
- **Precision.**
  - 1,000 reported modules per n and replicate (all if fewer) are sampled at random and re-tested with the independent
    test. Precision is the share that passes.
  - Exhaustive cross-check: the independent Python enumerator (undirected ESU over every weakly connected n-set) is run
    at n = 3 on all 200 replicates and at n = 4 on the first 5. Recall and precision against that complete truth set
    are reported.
- **Relaxations**, each on all 200 replicates:
  - **D3 relaxed:** exactly one source and one sink still required; "every vertex on a source–sink path" dropped.
  - **D4 relaxed:** one source–sink path suffices, instead of two internally vertex-disjoint paths.
  - Reported for each: the module count relative to the topological definition, the recall of the intact planted modules, and the share of
    reported modules that are planted.

## Software

C (`gcc -O3`) for the enumerator; Python 3 (numpy, networkx) for the graphs, the planting, the independent tests and the
summaries.

## Outputs

- `benchmark_per_replicate.csv`
- `benchmark_summary.csv`: medians and ranges over the 200 replicates
- `planted_topologies.csv`: the number of the topological definition topologies per n
