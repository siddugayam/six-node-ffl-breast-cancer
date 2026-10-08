# The regulatory network

The harmonised, directed, signed and evidence-tiered network of 587 nodes (157 TFs, 223 miRNAs and 207 genes). The
scripts read these files from `data/` under the analysis root.

| File | What it holds |
|---|---|
| `canonical_nodes.tsv` | Every node with its type |
| `canonical_edges.tsv` | Every edge (6,859; the analysed network excludes the 30 miRNA–miRNA edges of the exemplar circuits). Its `sign` column is set by regulator class (−1 for miRNA–target edges, +1 for all others); the literature signs are in the layer files and in `supplementary_tables/TableS2_all_interactions.csv` |
| `layer_TF_target.tsv`, `layer_TF_miRNA.tsv`, `layer_miRNA_miRNA.tsv`, `layer_gene_gene.tsv` | The four layers, from TRRUST v2, TransmiR v2.0, miRNA co-transcription and STRING v12, with their sources and signs |
| `edge_evidence_tier.tsv` | The multiMiR evidence tier of each miRNA–target edge (strong, weak, prediction only) |
| `mirna_id_map.tsv`, `name_map.json` | miRNA and gene identifier maps |
| `string_symbol_map_exact.tsv` | Gene symbol to STRING v12 protein identifier |

`scripts/01_network_assembly/02_build_canonical_network.py` builds `canonical_nodes.tsv`, `canonical_edges.tsv` and
`name_map.json` from the six networks of the original submission and their node tables
(`analyses/original_submission_code/`). The STRING tier of `layer_gene_gene.tsv` starts from the STRING v12.0 links with a
score of at least 900 among the 362 network proteins mapped to STRING: 446 pairs in the physical subnetwork and 941 in
the full (functional) network. `analyses/six_node_pattern/S7/s7b_restricted_gene_gene.py` derives both pair sets from
the STRING v12.0 download files and records the comparison in `s7b_restricted_gene_gene.txt`.

For the network with literature signs, evidence tiers and correlations in one table, use
`supplementary_tables/TableS2_all_interactions.csv`.
