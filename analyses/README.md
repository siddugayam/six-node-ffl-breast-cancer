# Analyses

Each folder is one analysis added after the main rounds, with its code, the inputs it reads and the outputs it wrote,
in the layout in which it was run. Scripts refer to other analyses as `analyses/<name>/…` under the analysis root.

| Folder | Paper | What it holds |
|---|---|---|
| [`census_and_motif_nulls/`](census_and_motif_nulls) | Table 2, Figs 3a and S2, Note S2 | The census of 3- to 7-node FFL modules and the three-node motif null models, on four versions of the network |
| [`six_node_pattern/`](six_node_pattern) | Note S2 | Tests of the six-node composite pattern: over-representation, dynamics, signed circuits, robustness and survival |
| [`six_node_pattern_networks/`](six_node_pattern_networks) | Note S2 | The three- to six-node composite FFL networks rebuilt on the census graph, with topology and MCODE |
| [`bhat_pattern_analysis/`](bhat_pattern_analysis) | Note S2 | The twelve networks of Bhat et al. (2024), analysed with the settings of the paper: over-representation, dynamics with a TF2→TF1 arc, sensitivity to clusters and to the evidence filter |
| [`six_node_followups/`](six_node_followups) | Notes S1, S2 | Follow-up questions on the miRNA–TF pair filter and the six-node pattern |
| [`analysed_network_reruns/`](analysed_network_reruns) | Notes S4, S7 | Re-runs without the 30 exemplar miRNA–miRNA edges: architecture, hubs, Tables S2, S3 and S6, module detection |
| [`dynamics_controls/`](dynamics_controls) | Results 3.7 | Feedback-loop structure of the modelled topologies and the four-node TF–TF control |
| [`metabric_power_downsampling/`](metabric_power_downsampling) | Results 3.8 | METABRIC downsampled to the TCGA event count |
| [`perturbation_tests/`](perturbation_tests) | Note S1 | Tests with public perturbation data, purified-cell miRNA atlases and a second sequence model |
| [`perturbation_sensitivity/`](perturbation_sensitivity) | Notes S1, S2 | Sensitivity analyses of the perturbation tests |
| [`original_submission_code/`](original_submission_code) | — | The network files and R scripts of the original submission (MIT licence) |

## Inputs read from outside these folders

- `six_node_pattern/E/E6` and `S8` read the TCGA-BRCA R objects `data/brca_*.rds`, built by
  `scripts/07_expression_validation/tcga_differential_expression/01_tcga_brca_prep_de.R` and `scripts/10_survival_and_clinical/cox_models/10_survival_cox_hubs.R` from the downloads listed in the main README.
- `six_node_pattern/S7/s7b_restricted_gene_gene.py` reads the STRING v12 downloads, and `six_node_followups/Q7` the
  TargetScan 8.0 miRNA family file.
- `six_node_pattern/E/E5/e5a_hypergeometric_filter.py` and `six_node_followups/Q123`, `Q5` and `Q9` take the authors'
  pair-filter archive as an argument; its files are in `data/pair_filter/`.
- The two scripts of `six_node_pattern/S5` read the edge list of `results/v5/tables/TableS1_all_interactions.csv` (the
  published Table S2), not its sign columns.
- The perturbation-data tests read, besides public downloads, the network tables,
  `results/v3/seqreg_ext_occlusion_allsites.csv`, TRRUST v2 and, for `P2/p2_analyse.py`, the six-node composite
  instances of `six_node_pattern/S1`.
