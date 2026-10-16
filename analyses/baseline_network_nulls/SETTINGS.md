# SETTINGS: a baseline network without FFL pre-selection

Settings of the analysis.

## Why a substitute

The 1,057 curated genes and 421 curated miRNAs from which the network was assembled were selected by hand from GeneCards, DisGeNET,
PhenomiR, HMDD and miR2Disease. The selected lists were not saved. The closest
substitute keeps the 587 nodes of the analysed network (157 TFs, 207 other genes, 223 miRNAs) and takes every edge
among them directly from the primary interaction tables, not from the exported FFL networks.

## Edges among the 587 nodes, from the primary tables

- **TF -> gene or TF.** TRRUST version 2 human (`data/db/trrust_human.tsv`): every row whose regulator and target are network
  protein-coding nodes, self-loops excluded.
  - Check: the arcs equal those of `data/network/layer_TF_target.tsv`, the project's TRRUST layer for the same node
    universe, without its self-loops.
- **TF -> miRNA.** TransmiR version 2.0 human (`data/db/transmir_hsa.tsv`), miRNA names canonicalised with
  `scripts/01_network_assembly/01_mirna_canon.R`.
  - Check: the arcs equal those of `data/network/layer_TF_miRNA.tsv`.
- **miRNA -> gene or TF, validated:** the multiMiR validated tables (miRTarBase, TarBase, miRecords; the cached query
  of the evidence-tier script). Rows marked non-functional by miRTarBase are excluded, as in the tier script; names are
  canonicalised the same way.
- **miRNA -> gene or TF, predicted:** TargetScan 8.0 default predictions (`Predicted_Targets_Info.default_predictions.txt`,
  human).
  - Families are mapped to their human mature miRNAs with `miR_Family_Info.txt`, then canonicalised.
  - Gene symbols must equal a network node.

## Networks

| network | miRNA -> target layer |
|---|---|
| validated | validated only |
| validated plus predicted | validated, plus the TargetScan predictions |

Each network is run with two typings: published (157 TFs), and the 22 co-regulators of the analysed network typed as
genes (the list of the tier- and sign-aware null analysis). Reciprocal TF<->miRNA pairs arise wherever both directions
are present.

## Null models and test

- `analyses/census_and_motif_nulls/nulls/motif_null_models.c`, unchanged:
  - the degree-preserving null: swaps within each source-type -> target-type class;
  - the reciprocity-preserving null: as above, with both arcs of every reciprocal TF<->miRNA pair frozen;
  - the bipartite-layer null: reciprocity-preserving, plus curveball randomisation of the miRNA -> target and
    TF -> target layers.
- 1,000 randomised networks each; seed 20250908; 100 swap attempts or trades per arc.
- **Counts:** all three-node FFLs; composite FFLs (once per reciprocal pair); TF-FFLs; miRNA-FFLs.
- **Reported:** observed; expected (mean); fold; excess %; p_upper = (1 + #{null ≥ observed}) / 1,001; p_lower likewise.
- **Analysed-network comparison**, shown beside the substitute:
  - published typing (`results/expression_validation/census_tables/motif_nulls_all_graphs.csv`, analysed network
    without the exemplar edges): +7.9 % reciprocity-preserving, +11.2 % bipartite-layer;
  - retyped (`analyses/six_node_pattern/E/E5/e5b_null_summary.csv`).

## Rule (per network × typing × metric)

- **"persists":** p_upper < 0.05 under both the reciprocity-preserving and the bipartite-layer null.
- **"reverses":** p_lower < 0.05 under both.
- **"vanishes":** otherwise.
- The degree-preserving null is reported, but it is not part of the rule (it destroys reciprocity).

## Also reported

- Arcs per class in each network.
- The overlap with the analysed network: its arcs that the primary tables reproduce, and the arcs added here.

## Outputs

- `network_<name>_<typing>.txt` (null format) and `network_arc_counts.csv`
- `analysed_network_overlap.csv`
- `null_model_summary.csv`
- `decision.csv`
- `runs/`

## Notes on the program

- `motif_null_models.c` holds at most 20,000 arcs (MAXE), and the baseline networks have
  30,604 and 32,868. The runs therefore use `motif_null_models_large.c`, a copy that differs only in MAXE (40,000).
