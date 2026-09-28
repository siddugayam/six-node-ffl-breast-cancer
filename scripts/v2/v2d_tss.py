#!/usr/bin/env python3
"""Build hg38 TSS table for every gene/TF node of the canonical network from
UCSC ncbiRefSeqSelect (one MANE/representative transcript per gene)."""
import csv, gzip, json, os
BASE="/path/to/revision"
nodes=list(csv.DictReader(open(BASE+"/data/canonical_nodes.tsv"),delimiter="\t"))
want={n["name"] for n in nodes if n["type"] in ("Gene","TF")}
print("gene/TF nodes:",len(want))
MAIN={"chr%s"%c for c in list(range(1,23))+["X","Y","M"]}
tss={}
with gzip.open(BASE+"/cache/direct/ncbiRefSeqSelect.txt.gz","rt") as f:
    for line in f:
        p=line.rstrip("\n").split("\t")
        name2=p[12]; chrom=p[2]; strand=p[3]
        if chrom not in MAIN: continue
        if name2 not in want: continue
        s=int(p[4]); e=int(p[5])
        t = e if strand=="-" else s
        tss[name2]=(chrom,t,strand,p[1])
print("nodes with TSS:",len(tss),"missing:",len(want-set(tss)))
miss=sorted(want-set(tss)); print("missing list:",miss)
json.dump(tss,open(BASE+"/cache/direct/node_tss_hg38.json","w"),indent=0)
