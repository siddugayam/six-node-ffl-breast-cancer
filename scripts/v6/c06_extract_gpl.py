#!/usr/bin/env python3
"""Extract probe -> miRNA-name maps from GEO platform SOFT files."""
import gzip, os, sys, csv, re
G="/path/to/revision/cache/v6/cohorts/gpl"
OUT=os.path.join(G,"maps"); os.makedirs(OUT, exist_ok=True)
NAMECOLS=["miRNA_ID","Gene Symbol","SYMBOL","miRNA_ID_LIST","TargetMatureName","ILMN_Gene","SPOT_ID","detector","ACCESSION_LIST","Human Mature Name","Probe Name"]
for f in sorted(os.listdir(G)):
    if not f.endswith("_family.soft.gz"): continue
    gpl=f.split("_")[0]
    hdr=None; rows=[]
    with gzip.open(os.path.join(G,f),"rt",encoding="utf-8",errors="replace") as fh:
        intab=False
        for line in fh:
            if line.startswith("!platform_table_begin"): intab=True; continue
            if line.startswith("!platform_table_end"): break
            if intab:
                p=line.rstrip("\n").split("\t")
                if hdr is None: hdr=p; continue
                rows.append(p)
    if hdr is None: print(gpl,"NO TABLE"); continue
    idx={c:i for i,c in enumerate(hdr)}
    have=[c for c in NAMECOLS if c in idx]
    op=os.path.join(OUT,gpl+"_map.tsv")
    with open(op,"w",newline="") as o:
        w=csv.writer(o,delimiter="\t"); w.writerow(["probe"]+have)
        for r in rows:
            w.writerow([r[0]]+[ (r[idx[c]] if idx[c]<len(r) else "") for c in have])
    ids=[r[0] for r in rows]
    print("%s  nprobes=%d  cols=%s  id_examples=%s" % (gpl,len(rows),have,ids[:3]))
