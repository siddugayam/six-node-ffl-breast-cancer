# Paper map

For each figure and table of the paper: the script that draws or writes it and the files that hold its values. For the
Results sections, `ANALYSIS_GUIDE.ipynb` lists the scripts of each analysis in the order they ran and loads the files
behind the reported numbers; the last table below gives its section for each part of the paper.

## Main figures

| Figure | Script | Data |
|---|---|---|
| Fig. 1 Study design, feed-forward loop classes and the topological definition of the *n*-node FFL | `scripts/15_figures_and_tables/74_fig1_concept.py` | none (schematic) |
| Fig. 2 The exemplar circuits, drawn as directed signed graphs | `scripts/15_figures_and_tables/72_fig2_exemplar_circuits.R` | `analyses/original_submission_code/SIF_files/`, `data/network/canonical_nodes.tsv`, `data/network/name_map.json` |
| Fig. 3 Over-representation of three-node FFLs and of the six-node composite pattern under null models | `scripts/15_figures_and_tables/71_fig3_nulls_sixnode.R` | (a) `results/v2/census_rerun/motif_nulls_all_graphs.csv`; (b) `results/v6/fig3b_six_node_pattern.csv` |
| Fig. 4 Stromal mediation, compartment occupancy and the hsa-miR-29a association with survival | `scripts/15_figures_and_tables/73_fig4_compartment_mir29a.R` | (a) `results/v3/deconv_mediation_extended.csv`; (b) `results/v6/spatial_compartment_marker_enrichment_meta.csv`; (c) `results/v6/pdx_matched_pairs_rank_change.csv`; (d) `results/v6/spatial_tf_collagen_correlation_meta.csv`; (e) `results/v6/mirna_meta_forest_table.csv` |
| Fig. 5 Behaviours gained and lost by higher-order modules relative to their three-node cores | `scripts/15_figures_and_tables/62_fig5_gain_loss.R` | `results/v3/dynamics_higher_order_gain_loss.csv` |

## Supplementary figures

| Figure | Script | Data |
|---|---|---|
| Fig. S1 The regulatory network and the interaction layers added for the census | `scripts/15_figures_and_tables/20_figures_core.R` | `data/network/canonical_nodes.tsv`, `canonical_edges.tsv` and the `layer_*.tsv` files |
| Fig. S2 Census of *n*-node feed-forward loops | `scripts/15_figures_and_tables/70_fig_supp_census_coherence_concordance.R` | `results/v2/census_rerun/census_all_graphs.csv` |
| Fig. S3 Coherence typing of the typed three-node cores | `scripts/15_figures_and_tables/70_fig_supp_census_coherence_concordance.R` | `results/v2/ffl_cores_coherence_corrected.csv` |
| Fig. S4 Prioritisation of regulatory nodes | `scripts/15_figures_and_tables/13_fig_prioritisation.R` | `results/v5/node_prioritisation_full.csv` |
| Fig. S5 Undirected and directed centrality of network nodes | `analyses/analysed_network_reruns/figure_fixes/fig_s5_information_flow_notitles.py` | the analysed-network outputs of `analyses/analysed_network_reruns/v3/` |
| Fig. S6 Sign concordance of predicted edges with co-expression | `scripts/15_figures_and_tables/70_fig_supp_census_coherence_concordance.R` | `results/sign_concordance_summary.csv` |
| Fig. S7 FFL module scores, the Farmer stroma-related signature and pathological complete response | `scripts/15_figures_and_tables/61_figS7_reactive_stroma.R` | `results/v5/farmer_module_score_correlation.csv`, `farmer_correlations_tcga.csv`, `farmer_neoadjuvant_meta.csv` |

## Tables

| Table | Scripts | Data |
|---|---|---|
| Table 1 Composition of the analysed network by edge class | counted from Table S2 (`ANALYSIS_GUIDE.ipynb`, section 1) | `supplementary_tables/TableS2_all_interactions.csv` (rows with `in_analysed_network` TRUE) |
| Table 2 Over-representation of three-node FFLs | `analyses/census_and_motif_nulls/run_all.sh` (`nulls/v2_null.c`) → `nulls/summarise_nulls.py` | `results/v2/census_rerun/motif_nulls_all_graphs.csv` |
| Table 3 Mediation by stromal content | `scripts/08_stroma_and_cell_types/deconvolution_and_mediation/46_deconv2_mediation.R` → `analyses/six_node_pattern/E/E7/e7_mediation_intervals.py` | `results/v3/deconv_mediation_extended.csv` |
| Table 4 Behaviours of higher-order modules | `scripts/11_dynamics/03_higher_order_sweep.py`, `03d_higher_order_compI1.py` | `results/v3/dynamics_higher_order_gain_loss.csv` |
| Table S1 Node prioritisation | `scripts/06_node_prioritisation/01_prioritisation.R` | `supplementary_tables/TableS1_node_prioritisation_full.csv` |
| Table S2 Complete interaction list | `scripts/15_figures_and_tables/21_tables.R` → `analyses/analysed_network_reruns/tables/flag_tableS1.py` | `supplementary_tables/TableS2_all_interactions.csv` |
| Table S3 Hub bootstrap stability | `analyses/analysed_network_reruns/hubboot/31_hub_bootstrap_nolegacy.py` | `supplementary_tables/TableS3_hub_bootstrap_stability.csv` |
| Table S4 Typed three-node FFL cores | `scripts/03_ffl_census/01b_coherence_corrected.py` → `scripts/15_figures_and_tables/21_tables.R` | `supplementary_tables/TableS4_ffl_cores.csv` |
| Table S5 Network partners of the prioritised nodes | `scripts/15_figures_and_tables/21_tables.R` | `supplementary_tables/TableS5_prioritised_node_partners.csv` |
| Table S6 Direction-aware hub table | `analyses/analysed_network_reruns/tables/rebuild_tableS6.py` | `supplementary_tables/TableS6_directed_hub_table.csv` |
| Table S7 Cohort inventory | `scripts/07_expression_validation/external_cohorts/N8_inventory.R` | `supplementary_tables/TableS7_cohort_inventory.csv` |
| Table S8 miRNA survival meta-analysis | `scripts/10_survival_and_clinical/mirna_survival_meta_analysis/c11_mirna_cox_meta.R` → `analyses/six_node_pattern/G/g1_bh_eleven_mirnas.py` | `supplementary_tables/TableS8_mirna_survival_meta.csv` |
| Table S9 Module-detection outcomes | summarised from the outputs of `scripts/14_module_detection/` and `analyses/analysed_network_reruns/` (`v7/`, `v7b/`) | `supplementary_tables/TableS9_module_detection_outcomes.csv` |

## Results sections and Supplementary Notes

| Part of the paper | Section of `ANALYSIS_GUIDE.ipynb` |
|---|---|
| Methods 2.1, Results 3.1 (network) | 1 Network assembly; 2 The miRNA–TF pair filter |
| Methods 2.2, Results 3.2 (cores, census, null models, six-node pattern) | 3 Three-node FFLs; 4 Census; 5 Motif significance; 6 The six-node composite pattern |
| Methods 2.3, Results 3.3 (prioritisation) | 7 Node prioritisation |
| Methods 2.4, Results 3.4 (expression validation) | 8 Expression validation and replication |
| Methods 2.6, Results 3.5 (stromal mediation, compartments, regulatory DNA) | 9 Stromal mediation and compartment tests; 10 Regulatory DNA |
| Results 3.6 (miR-29, reactive stroma, survival) | 11 The miR-29–collagen axis and clinical association |
| Methods 2.5, Results 3.7 (dynamics) | 12 Dynamics of higher-order modules |
| Results 3.8 (prognosis, essentiality, enrichment) | 13 Prognosis, essentiality and enrichment |
| Supplementary Note S7 (module detection) | 14 Module-detection methods |
| Supplementary Note S1 (public perturbation data) | 15 Tests with public perturbation data |
