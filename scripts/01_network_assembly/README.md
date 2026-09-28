# Network assembly (Methods 2.1)

`00_mirna_canon.R` (miRNA name helper) → `01_build_canonical_network.py` → the layers `02_layer_TF_target.R` (TRRUST v2),
`03_layer_TF_miRNA.R` (TransmiR v2.0), `04_layer_miRNA_miRNA.R` (co-transcribed miRNA pairs), `05a_string_map_and_fetch.R`
and `05_layer_gene_gene.R` (STRING v12) → `06_edge_evidence_tier.R` (multiMiR evidence tiers) → `07_layer_summary.R`.
The `v2_*` scripts re-derive the network from the original deposit and check its claims. Outputs: `data/network/`.
