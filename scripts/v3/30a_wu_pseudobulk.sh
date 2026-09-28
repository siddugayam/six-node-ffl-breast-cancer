#!/bin/bash
# Stream the Wu et al. 2021 (GSE176078) 29,733 x 100,064 count matrix once and accumulate
# summed UMI per gene for three groupings: celltype_minor (CAF subsets), patient|celltype_major
# and patient|celltype_minor.  Labels have had whitespace replaced by "_" so that awk's default
# whitespace field splitting is safe for both the map and the MatrixMarket file.
set -e
WU=/path/to/revision/cache/celltype/Wu_etal_2021_BRCA_scRNASeq
OUT=/path/to/revision/cache/v7/wu
awk -v OFS='\t' -v O="$OUT" '
 NR==FNR { if (FNR>1) { mi[$1]=$2; pm[$1]=$4; pmi[$1]=$5 } ; next }
 FNR<=2 { next }
 { g=$1; c=$2; v=$3; A[g SUBSEP mi[c]]+=v; B[g SUBSEP pm[c]]+=v; C[g SUBSEP pmi[c]]+=v }
 END {
   for (k in A) { split(k,a,SUBSEP); print a[1], a[2], A[k] > (O "/pb_minor.tsv") }
   for (k in B) { split(k,b,SUBSEP); print b[1], b[2], B[k] > (O "/pb_patmajor.tsv") }
   for (k in C) { split(k,d,SUBSEP); print d[1], d[2], C[k] > (O "/pb_patminor.tsv") }
 }' $OUT/cellmap.tsv $WU/count_matrix_sparse.mtx
wc -l $OUT/pb_minor.tsv $OUT/pb_patmajor.tsv $OUT/pb_patminor.tsv
