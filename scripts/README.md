# Scripts

The code of the main analysis, grouped by analysis in the order of the paper. Within a folder, scripts run in the order
of their numeric prefixes. `ANALYSIS_GUIDE.ipynb` (repository root) names the scripts behind each reported number, and
`docs/script_index.tsv` gives every script's earlier path. The results keep the folders the scripts write to
(`results/`, `results/v2/` … `results/v7/`).

| Folder | Paper | What it holds |
|---|---|---|
| [`01_network_assembly/`](01_network_assembly) | Methods 2.1, Results 3.1 | Network assembly from TRRUST, TransmiR, multiMiR and STRING; evidence tiers |
| [`02_pair_filter/`](02_pair_filter) | Methods 2.1 | The hypergeometric miRNA–TF pair filter, recomputed from `data/pair_filter/` |
| [`03_ffl_census/`](03_ffl_census) | Methods 2.2, Results 3.2, Fig. S2, Note S2 | Three-node cores and coherence; census of 3- to 7-node modules |
| [`04_motif_significance/`](04_motif_significance) | Results 3.2 | Null models and motif over-representation; node sets |
| [`05_network_architecture/`](05_network_architecture) | Note S4 | Degree distribution, controllability, information flow, knockouts |
| [`06_node_prioritisation/`](06_node_prioritisation) | Results 3.3, Notes S3, S5 | Hubs, ExIR, prioritisation, node compendium and literature |
| [`07_expression_validation/`](07_expression_validation) | Results 3.4 | Differential expression, edge correlations, external cohorts, CPTAC and RPPA |
| [`08_stroma_and_cell_types/`](08_stroma_and_cell_types) | Results 3.5 | Deconvolution and mediation, single-cell, spatial, xenograft, protein atlas, cell lines |
| [`09_regulatory_evidence/`](09_regulatory_evidence) | Results 3.5, Note S1 | Sequence models, ENCODE accessibility, ChIP, ReMap and TargetScan evidence |
| [`10_survival_and_clinical/`](10_survival_and_clinical) | Results 3.6, 3.8 | Cox models, miRNA survival meta-analysis, reactive stroma and neoadjuvant response |
| [`11_dynamics/`](11_dynamics) | Methods 2.5, Results 3.7 | Dynamical models of three-node and higher-order modules |
| [`12_enrichment_and_metabolism/`](12_enrichment_and_metabolism) | Results 3.8, Note S1 | ORA, GSEA and the metabolic hypotheses |
| [`13_multiomics_and_screens/`](13_multiomics_and_screens) | Results 3.8 | Mutation, copy number, methylation, DepMap, drug and CRISPR screens, GWAS |
| [`14_module_detection/`](14_module_detection) | Note S7 | Community detection, MCODE, BioNet, jActiveModules |
| [`15_figures_and_tables/`](15_figures_and_tables) | — | Figure and table scripts and the shared plotting theme |

Later, self-contained analyses (with their own inputs and outputs) are in `analyses/`. C programs are given as source;
see the repository README for compile commands.
