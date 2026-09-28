#!/usr/bin/env python3
"""Stream the 592 MB SCAN-B (GSE96058) expression CSV and keep only the rows
needed for the replication analysis (hubs + CAF signatures + collagens/TFs).
Writes cache/newcohorts/GSE96058_subset.csv"""
import gzip, sys, csv, os
BASE="/path/to/revision"
CA=os.path.join(BASE,"cache/newcohorts")

want=set()
with open(os.path.join(BASE,"results/network_topology_hubs.csv")) as f:
    for r in csv.DictReader(f):
        if r["is_hub"]=="TRUE" and r["type"]!="miRNA": want.add(r["name"])
n_hub=len(want)
with open(os.path.join(CA,"cafA_signature_genes.txt")) as f:
    caf=[l.strip() for l in f if l.strip()]
want|=set(caf)
want|={"DCN","LUM","FAP","THY1","COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1"}
print(f"hubs={n_hub} cafA={len(caf)} total wanted symbols={len(want)}", flush=True)

src=os.path.join(CA,"GSE96058_expr.csv.gz")
out=os.path.join(CA,"GSE96058_subset.csv")
found=set(); nrow=0
with gzip.open(src,"rt") as fi, open(out,"w") as fo:
    hdr=fi.readline()
    fo.write(hdr)
    ncol=hdr.count(",")
    print("header columns (samples) =", ncol, flush=True)
    for line in fi:
        nrow+=1
        i=line.find(",")
        g=line[:i].strip().strip('"')
        if g in want:
            found.add(g); fo.write(line)
print(f"rows scanned={nrow}  rows kept={len(found)}", flush=True)
missing=sorted(want-found)
print(f"symbols not found in SCAN-B ({len(missing)}): {','.join(missing)}", flush=True)
