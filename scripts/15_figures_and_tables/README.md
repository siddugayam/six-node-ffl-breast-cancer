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

## Scripts in this folder

| Script | What it does |
|---|---|
| `06_crosscheck_table.py` | HPA Pathology-Atlas breast prognostic calls vs our own TCGA-BRCA and METABRIC Cox results. |
| `10_master_table.py` | One master table: the 20 protein-coding prioritised nodes (10 TFs + 10 genes), HPA v25.1 - breast-cancer IHC, normal-breast IHC by cell type, pathology-atlas prognostics, subcellular location and antibody reliability. |
| `10_theme.R` | Shared publication theme for the figures. |
| `13_fig_prioritisation.R` | Fig. S4: the 30 prioritised nodes, their evidence-domain profiles, and the divergence between ranking by FFL topology and ranking by evidence. |
| `20_figures_core.R` | Fig. S1 (the global regulatory network, FigureS1_global_network), with a census panel and directed diagrams of the exemplar modules; Figs S2 and 2 of the paper are drawn by scripts 70 and 72. |
| `21_tables.R` | Complete interaction and supplementary tables. |
| `61_figS7_reactive_stroma.R` | Fig. S7: FFL module scores against the Farmer stroma-related signature, and the pCR meta-analysis. |
| `62_fig5_gain_loss.R` | Fig. 5: behaviours gained and lost by higher-order modules relative to their embedded three-node cores, for the composite circuit (a) and the incoherent type-1 circuit I1 (b). |
| `70_fig_supp_census_coherence_concordance.R` | Supplementary Figs S2, S3 and S6: S2 = census and edge-class saturation; S3 = coherence typology; S6 = sign concordance by edge class and evidence tier. |
| `71_fig3_nulls_sixnode.R` | Fig. 3: (a) over-representation of three-node FFLs and of reciprocal TF-miRNA pairs against the three null models; (b) the six-node composite pattern of Bhat et al. |
| `72_fig2_exemplar_circuits.R` | Fig. 2: the exemplar four-, five- and six-node circuits as directed signed graphs (drawn at print width with 8 pt labels, italic gene symbols, one shared key). |
| `73_fig4_compartment_mir29a.R` | Fig. 4: (a) mediation of the TF-collagen and miR-29-collagen associations by stromal content across 42 estimates; (b) spatial compartment occupancy, (c) patient-to-xenograft change and (d) bulk versus within-stroma TF-COL1A1 … |
| `74_fig1_concept.py` | Fig. 1: workflow (a), the three FFL classes and the edge-addition series from three to six nodes (b), and the topological definition of an n-node FFL (c). |
