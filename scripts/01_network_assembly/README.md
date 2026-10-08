# Network assembly (Methods 2.1)

`01_mirna_canon.R` (miRNA name helper) → `02_build_canonical_network.py` → the layers `03_layer_TF_target.R` (TRRUST v2),
`04_layer_TF_miRNA.R` (TransmiR v2.0), `05_layer_miRNA_miRNA.R` (co-transcribed miRNA pairs), `06a_string_map_and_fetch.R`
and `06_layer_gene_gene.R` (STRING v12) → `07_edge_evidence_tier.R` (multiMiR evidence tiers) → `08_layer_summary.R`.
Outputs: `data/network/`. `02_build_canonical_network.py` reads the six networks of the original submission and their
node tables (`analyses/original_submission_code/SIF_files/` and `node_attributes/`). `06a_string_map_and_fetch.R` runs
before `06_layer_gene_gene.R` and maps the network genes to STRING v12.0 proteins. `06_layer_gene_gene.R` then reads the
STRING v12.0 links with a score of at least 900 among those proteins (`cache/string_physical_900_ensp.tsv`, 446 pairs;
`cache/string_functional_900_ensp.tsv`, 941 pairs). `analyses/six_node_pattern/S7/s7b_restricted_gene_gene.py` derives
the same two pair sets from the STRING v12.0 download files.

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_mirna_canon.R` | Shared helper: canonicalise a miRBase name (precursor or mature, any case) to the canonical mature stem form used by canonical_nodes.tsv: hsa-mir-155 -> hsa-miR-155 hsa-mir-29b-1 -> hsa-miR-29b (strip precursor copy-number … |
| `02_build_canonical_network.py` | Rebuild a single, canonical, DIRECTED, SIGNED, type-annotated regulatory network from the six Cytoscape exports in miRNA_Github_GPR/. |
| `03_layer_TF_target.R` | A) TF-TF and TF-gene layer with literature-derived signs, from TRRUST v2 (human). |
| `04_layer_TF_miRNA.R` | B) TF-miRNA layer with literature-derived signs, from TransmiR v2 (human). |
| `05_layer_miRNA_miRNA.R` | C) miRNA-miRNA layer defined as POLYCISTRONIC CO-TRANSCRIPTION (genomic clustering), NOT as invented direct regulatory interactions. |
| `06_layer_gene_gene.R` | D) Gene-gene layer: two clearly separated evidence tiers. |
| `06a_string_map_and_fetch.R` | Exact symbol -> STRING v12 protein-ID mapping using the official STRING info/aliases flat files (the /get_string_ids API does FUZZY text matching and mis-resolves e.g. |
| `07_edge_evidence_tier.R` | E) Evidence tier for every miRNA_target edge in the canonical network. |
| `08_layer_summary.R` | Summary of the interaction layers: edges by sign or source for each layer (TRRUST, TransmiR, miRNA clusters at 3, 10 and 50 kb, gene-gene sources and the miRNA-target evidence tiers), entered from the outputs of the layer scripts … |
