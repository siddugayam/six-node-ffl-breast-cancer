# Network assembly (Methods 2.1)

`01_mirna_canon.R` (miRNA name helper) → `02_build_canonical_network.py` → the layers `03_layer_TF_target.R` (TRRUST v2),
`04_layer_TF_miRNA.R` (TransmiR v2.0), `05_layer_miRNA_miRNA.R` (co-transcribed miRNA pairs), `06a_string_map_and_fetch.R`
and `06_layer_gene_gene.R` (STRING v12) → `07_edge_evidence_tier.R` (multiMiR evidence tiers) → `08_layer_summary.R`.
Outputs: `data/network/`.
