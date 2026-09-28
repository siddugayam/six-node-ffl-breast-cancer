#!/bin/bash
# Single streaming pass over the 178M-entry GSE176078 sparse matrix.
# Emits (a) per-cell total UMI counts, (b) counts for the genes of interest only.
set -euo pipefail
D=/path/to/revision/cache/celltype
M=$D/Wu_etal_2021_BRCA_scRNASeq
awk -v WANT="$D/scrna_want_genes.txt" -v TOT="$D/scrna_cell_totals.tsv" '
BEGIN{ while((getline g < WANT) > 0) want[g]=1 }
NR==FNR && FNR>0 { gi[FNR]=$1; if($1 in want) keep[FNR]=1; next }
FNR<=2 { next }
{ tot[$2]+=$3; if ($1 in keep) print gi[$1]"\t"$2"\t"$3 }
END{ for(c in tot) print c"\t"tot[c] > TOT }
' "$M/count_matrix_genes.tsv" "$M/count_matrix_sparse.mtx" > "$D/scrna_gene_counts.tsv"
wc -l "$D/scrna_gene_counts.tsv" "$D/scrna_cell_totals.tsv"
