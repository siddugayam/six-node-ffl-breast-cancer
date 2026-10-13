# Tcga differential expression

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_de_limma.R` | Independent re-derivation of the TCGA-BRCA differential expression tables with limma. |
| `01_tcga_brca_prep_de.R` | TCGA-BRCA expression matrices + differential expression (tumour vs normal) Shared substrate for downstream FFL validation analyses. |
| `01b_gene_de_node_match.R` | The pan-cancer gene matrix is keyed by HGNC SYMBOL (a minority of rows are bare Entrez IDs for loci with no symbol). |
