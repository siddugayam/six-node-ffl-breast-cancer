#!/bin/bash
# Extract the rows we need from the 331MB pan-cancer Xena gene-expression matrix.
set -euo pipefail
ML=/path/to/home/Desktop/DD/R_GPR/ML
R=/path/to/revision
C=$R/cache/v6/pancancer
mkdir -p "$C"
# 1) build the entrez id list: COL1A1, COL3A1, EZH2 + CAF-A signature + controls
python3 - <<'PY'
import csv
R="/path/to/revision"
sym2id={}
with open(R+"/data/geneinfo_beta.txt") as f:
    r=csv.DictReader(f,delimiter="\t")
    for row in r:
        s=row["gene_symbol"]
        if s and s not in sym2id: sym2id[s]=row["gene_id"]
want=["COL1A1","COL3A1","EZH2","COL1A2","POSTN","FN1","PDGFRB","ACTA2","DCN","LUM",
      "PTPRC","PECAM1","EPCAM","KRT8","KRT18","VIM","MKI67","ESR1","ERBB2","FAP","THY1",
      "CDH1","MYC","GAPDH","ACTB"]
caf=[l.strip() for l in open(R+"/results/v2/cafA_gene_list_v3.txt") if l.strip()]
allw=list(dict.fromkeys(want+caf))
out=[]; miss=[]
for s in allw:
    if s in sym2id: out.append((sym2id[s],s))
    else: miss.append(s)
with open(R+"/cache/v6/pancancer/wanted_genes.tsv","w") as o:
    o.write("entrez\tsymbol\n")
    for i,s in out: o.write("%s\t%s\n"%(i,s))
print("wanted", len(allw), "mapped", len(out), "unmapped", miss)
PY
cut -f1 "$C/wanted_genes.tsv" | tail -n +2 > "$C/ids.txt"
zcat "$ML/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz" \
  | awk -F'\t' 'NR==FNR{ids[$1]=1;next} FNR==1||($1 in ids)' "$C/ids.txt" - \
  > "$C/pancan_gene_subset.tsv"
wc -l "$C/pancan_gene_subset.tsv"
