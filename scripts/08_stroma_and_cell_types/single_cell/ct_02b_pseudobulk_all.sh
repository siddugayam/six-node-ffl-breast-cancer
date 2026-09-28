#!/bin/bash
# Second streaming pass: per-(gene x celltype_major) UMI sums and n-cells-expressing,
# for ALL 29,733 genes. Used to derive a data-driven, collagen-free CAF signature.
set -euo pipefail
D=/path/to/revision/cache/celltype
M=$D/Wu_etal_2021_BRCA_scRNASeq
mawk -v CT="$D/cell_to_ctcode.tsv" '
BEGIN{ while((getline line < CT) > 0){ split(line,a,"\t"); ct[a[1]+0]=a[2]+0 } }
NR==FNR { gi[FNR]=$1; next }
FNR<=2 { next }
{ c=ct[$2+0]; k=$1 SUBSEP c; s[k]+=$3; n[k]++ }
END{ print "gene\tct_code\tsum_umi\tn_cells_expressing";
     for(k in s){ split(k,b,SUBSEP); print gi[b[1]+0]"\t"b[2]"\t"s[k]"\t"n[k] } }
' "$M/count_matrix_genes.tsv" "$M/count_matrix_sparse.mtx" > "$D/scrna_pseudobulk_all.tsv"
wc -l "$D/scrna_pseudobulk_all.tsv"
