# Analyses

Each folder is one self-contained analysis, with its code, the inputs it reads and the outputs it wrote,
in the layout in which it was run. Scripts refer to other analyses as `analyses/<name>/…` under the analysis root.

| Folder | Paper | What it holds |
|---|---|---|
| [`census_and_motif_nulls/`](census_and_motif_nulls) | Table 2, Figs 3a and S2, Note S2 | The census of 3- to 7-node FFL modules and the three-node motif null models, on four versions of the network |
| [`six_node_pattern/`](six_node_pattern) | Note S2 | Tests of the six-node composite pattern: over-representation, dynamics, signed circuits, robustness and survival |
| [`six_node_pattern_networks/`](six_node_pattern_networks) | Note S2 | The three- to six-node composite FFL networks rebuilt on the census graph, with topology and MCODE |
| [`bhat_pattern_analysis/`](bhat_pattern_analysis) | Note S2 | The twelve networks of Bhat et al. (2024), analysed with the settings of the paper: over-representation, dynamics with a TF2→TF1 arc, sensitivity to clusters and to the evidence filter |
| [`six_node_followups/`](six_node_followups) | Notes S1, S2 | Follow-up questions on the miRNA–TF pair filter and the six-node pattern |
| [`reciprocal_handling/`](reciprocal_handling) | Note S2 | The census under three alternative treatments of the reciprocal TF↔miRNA pairs |
| [`six_node_acyclicity/`](six_node_acyclicity) | Results 3.2, Note S2 | The six-node composite pattern tested against the topological definition, and the effect of allowing one cyclic component |
| [`tier_sign_nulls/`](tier_sign_nulls) | Note S2 | The three-node FFL test under null models that also preserve evidence tier and edge sign |
| [`baseline_network_nulls/`](baseline_network_nulls) | Results 3.2, Discussion 4.8, Note S2 | A substitute network of the same 587 nodes with every edge taken from the primary interaction tables, under the three null models |
| [`string_orientation/`](string_orientation) | Results 3.2, Note S2 | The census with the 941 STRING associations oriented at random, gene–gene layer variants and the orientation-robust modules |
| [`stroma_adjusted_concordance/`](stroma_adjusted_concordance) | Results 3.4, Note S1 | Tier-stratified miRNA–target sign concordance against the matched null after stromal adjustment |
| [`concordance_combined_tests/`](concordance_combined_tests) | Results 3.4, Note S1 | The stromal-adjusted concordance tests under five rules for combining the matched-null draws |
| [`mediation_common_sample/`](mediation_common_sample) | Table 3, Note S1 | Stromal mediation of the twelve paths refitted on the 734 tumours shared by all 42 stromal estimates |
| [`multivariable_survival/`](multivariable_survival) | Results 3.6, Note S1 | Unadjusted and covariate-adjusted Cox models for the 11 miRNAs of the survival meta-analysis |
| [`er_pam50_adjustment/`](er_pam50_adjustment) | Results 3.6 and 3.8, Note S1 | Six-node circuit survival and the miRNA meta-analysis with oestrogen-receptor status and PAM50 subtype added |
| [`analysed_vs_full_network/`](analysed_vs_full_network) | Notes S4, S6 | Architecture, hubs, Tables S2, S3 and S6 and module detection on the analysed network, without the 30 exemplar miRNA–miRNA edges and with them |
| [`dynamics_controls/`](dynamics_controls) | Results 3.7 | Feedback-loop structure of the modelled topologies and the four-node TF–TF control |
| [`metabric_power_downsampling/`](metabric_power_downsampling) | Results 3.8 | METABRIC downsampled to the TCGA event count |
| [`perturbation_tests/`](perturbation_tests) | Note S1 | Tests with public perturbation data, purified-cell miRNA atlases and a second sequence model |
| [`perturbation_sensitivity/`](perturbation_sensitivity) | Notes S1, S2 | Sensitivity analyses of the perturbation tests |
| [`exemplar_circuits/`](exemplar_circuits) | — | The six FFL network files (three-node networks and four- to six-node exemplar circuits), their node tables and R scripts (MIT licence) |

## Inputs read from outside these folders

- `six_node_pattern/E/E6` and `S8` read the TCGA-BRCA R objects `data/brca_*.rds`, built by
  `scripts/07_expression_validation/tcga_differential_expression/01_tcga_brca_prep_de.R` and `scripts/10_survival_and_clinical/cox_models/10_survival_cox_hubs.R` from the downloads listed in the main README.
- `six_node_pattern/S7/s7b_restricted_gene_gene.py` reads the STRING v12 downloads, and `six_node_followups/Q7` the
  TargetScan 8.0 miRNA family file.
- `six_node_pattern/E/E5/e5a_hypergeometric_filter.py` and `six_node_followups/Q123`, `Q5` and `Q9` take the authors'
  pair-filter archive as an argument; its files are in `data/pair_filter/`.
- The two scripts of `six_node_pattern/S5` read the edge list of `results/prioritisation_and_stroma/tables/all_interactions.csv` (the
  published Table S2), not its sign columns.
- The perturbation-data tests read, besides public downloads, the network tables,
  `results/systems_and_dynamics/seqreg_ext_occlusion_allsites.csv`, TRRUST v2 and, for `P2/p2_analyse.py`, the six-node composite
  instances of `six_node_pattern/S1`.
