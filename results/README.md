# Results

The outputs of the scripts in `scripts/`, in the same layout: the top level holds the first round and `v2/` … `v7/` the
later rounds. `ANALYSIS_GUIDE.ipynb` names the file behind each reported number; `docs/files_cited_in_supplementary_notes.tsv`
gives the path of every file named in the Supplementary Notes.

| Folder | What it holds |
|---|---|
| top level | Network topology and hubs, three-node cores and coherence, the first census, motif significance, enrichment, edge correlations and survival |
| `v2/` | Coherence, the census graphs, differential expression, CAF-adjusted correlations, mediation, cell lines and external cohorts; `census_rerun/` holds the census and motif tables of Table 3 and Fig. 1 |
| `v3/` | Dynamics, architecture and information flow, sequence-level regulation, deconvolution and mediation, single-cell and screen results |
| `v4/` | Gene-set enrichment of the FFL classes and the metabolic analyses |
| `v5/` | Prioritisation, node compendium, literature counts, the Farmer reactive-stroma analyses and copies of the tables (`tables/`) |
| `v6/` | Spatial, xenograft, protein-atlas and ENCODE accessibility results, the miRNA survival meta-analysis and liquid biopsy |
| `v7/` | Module detection |
| `multiomics/` | Mutation, copy number, methylation, CPTAC, RPPA, DepMap, stromal adjustment, external cohorts and metabolomics |
| `figures/` | Plots kept with the results |

Files above 50 MB are gzip-compressed (`.gz`).
