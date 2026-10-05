# Results

The outputs of the scripts in `scripts/`, in the same layout: the top level holds the first round and `v2/` … `v7/` the
later rounds. `ANALYSIS_GUIDE.ipynb` names the file behind each reported number; `docs/files_cited_in_supplementary_notes.tsv`
gives the path of every file named in the Supplementary Notes.

| Folder | What it holds |
|---|---|
| top level | Network topology and hubs, three-node cores and coherence, the first census, motif significance, enrichment, edge correlations and survival |
| `v2/` | Coherence, the census graphs, differential expression, CAF-adjusted correlations, mediation, cell lines and external cohorts; `census_rerun/` holds the census and motif tables behind Table 2, Fig. 3a and Fig. S2 |
| `v3/` | Dynamics, architecture and information flow, sequence-level regulation, deconvolution and mediation, single-cell and screen results |
| `v4/` | Gene-set enrichment of the FFL classes and the metabolic analyses |
| `v5/` | Prioritisation, node compendium, literature counts, the Farmer reactive-stroma analyses and copies of the tables (`tables/`) |
| `v6/` | Spatial, xenograft, protein-atlas and ENCODE accessibility results, the miRNA survival meta-analysis and liquid biopsy |
| `v7/` | Module detection |
| `multiomics/` | Mutation, copy number, methylation, CPTAC, RPPA, DepMap, stromal adjustment, external cohorts and metabolomics |
| `figures/` | Plots kept with the results |

Files above 50 MB are gzip-compressed (`.gz`).

## Files behind the reported values

- Census: `v2/census_rerun/census_all_graphs.csv`. The top-level `ffl_census_summary.csv`, `ffl_census.csv` and
  `ffl_higher_order.csv.gz` count the graph that also holds the 30 exemplar miRNA–miRNA edges.
- Node ranks: Table S1 (`supplementary_tables/TableS1_node_prioritisation_full.csv`), ordered by its `priority` column;
  `ANALYSIS_GUIDE.ipynb` recomputes them. In `v5/node_prioritisation_full.csv`, rank by the same column.
- `v5/tables/` holds the tables under the analysis numbering; its README maps them to the published tables.
- `v7/mcode_19_headline_summary.csv` is compiled from the files `mcode_00` to `mcode_18` of the same folder, which hold
  its values.
- `v6/atac_encode_compartment_specificity_tests.csv` is rebuilt byte for byte from `v6/atac_encode_region_accessible.csv`
  by `scripts/09_regulatory_evidence/encode_accessibility/atac_compartment_specificity_tests_rebuild.py`, which writes the rebuilt copy and a summary next to itself.
