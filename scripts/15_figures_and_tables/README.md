# Figures and tables

`10_theme.R` holds the shared theme. Each figure of the paper is drawn by one script:

| Script | Figure |
|---|---|
| `74_fig1_concept.py` | Fig. 1 |
| `72_fig2_exemplar_circuits.R` | Fig. 2 |
| `71_fig3_nulls_sixnode.R` | Fig. 3 |
| `73_fig4_compartment_mir29a.R` | Fig. 4 |
| `62_fig5_gain_loss.R` | Fig. 5 |
| `20_figures_core.R` | Fig. S1 |
| `70_fig_supp_census_coherence_concordance.R` | Figs S2, S3 and S6 |
| `13_fig_prioritisation.R` | Fig. S4 |
| `analyses/analysed_network_reruns/figure_fixes/fig_s5_information_flow_notitles.py` | Fig. S5 |
| `61_figS7_reactive_stroma.R` | Fig. S7 |

`21_tables.R` writes the main-text and supplementary tables. Three steps after it give the published Table S2
(`supplementary_tables/TableS2_all_interactions.csv`): edges without an annotated sign have an empty `sign`; the 521
repression edges (136 TRRUST, 385 TransmiR) have `sign` −1, as their `trrust_mode` and `transmir_mode` state; and the
`in_analysed_network` column, written by `analyses/analysed_network_reruns/tables/flag_tableS1.py`, is FALSE for the 30
exemplar miRNA–miRNA edges. `10_master_table.py` (Human Protein Atlas data for the 20 protein-coding prioritised nodes)
and `06_crosscheck_table.py` (the Atlas's prognostic calls against the TCGA-BRCA and METABRIC Cox results) write the
protein-atlas tables.
