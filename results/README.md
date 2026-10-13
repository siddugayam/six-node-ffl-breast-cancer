# Results

The outputs of the scripts in `scripts/`. The top level holds the core network analyses, and each subfolder holds one
group of analyses. `ANALYSIS_GUIDE.ipynb` names the file behind each reported number; `docs/files_cited_in_supplementary_notes.tsv`
gives the path of every file named in the Supplementary Notes.

| Folder | What it holds |
|---|---|
| top level | Network topology and hubs, three-node cores and coherence, the census of the full graph, motif significance, enrichment, edge correlations and survival |
| `expression_validation/` | Coherence, the census graphs, differential expression, CAF-adjusted correlations, mediation, cell lines and external cohorts; `census_tables/` holds the census and motif tables behind Table 2, Fig. 3a and Fig. S2 |
| `systems_and_dynamics/` | Dynamics, architecture and information flow, sequence-level regulation, deconvolution and mediation, single-cell and screen results |
| `metabolism/` | Gene-set enrichment of the FFL classes and the metabolic analyses |
| `prioritisation_and_stroma/` | Prioritisation, node compendium, literature counts, the Farmer reactive-stroma analyses and copies of the tables (`tables/`) |
| `tissue_and_cohorts/` | Spatial, xenograft, protein-atlas and ENCODE accessibility results, the miRNA survival meta-analysis and liquid biopsy |
| `module_detection/` | Module detection |
| `multiomics/` | Mutation, copy number, methylation, CPTAC, RPPA, DepMap, stromal adjustment, external cohorts and metabolomics |
| `figures/` | Plots kept with the results |

Files above 50 MB are gzip-compressed (`.gz`).

## Files behind the reported values

- Census: `expression_validation/census_tables/census_all_graphs.csv`. The top-level `ffl_census_summary.csv`, `ffl_census.csv` and
  `ffl_higher_order.csv.gz` count the graph that also holds the 30 exemplar miRNA–miRNA edges.
- Node ranks: Table S1 (`supplementary_tables/TableS1_node_prioritisation_full.csv`), ordered by its `priority` column;
  `ANALYSIS_GUIDE.ipynb` recomputes them. In `prioritisation_and_stroma/node_prioritisation_full.csv`, rank by the same column.
- `prioritisation_and_stroma/tables/` holds the tables under the analysis numbering; its README maps them to the published tables.
- `module_detection/mcode_19_headline_summary.csv` is compiled from the files `mcode_00` to `mcode_18` of the same folder, which hold
  its values.
- `tissue_and_cohorts/atac_encode_compartment_specificity_tests.csv` is rebuilt byte for byte from `tissue_and_cohorts/atac_encode_region_accessible.csv`
  by `scripts/09_regulatory_evidence/encode_accessibility/atac_06c_compartment_specificity_tests.py`, which writes the rebuilt copy and a summary next to itself.
