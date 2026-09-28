#!/usr/bin/env python3
"""Genome-wide background: ReMap 2022 promoter occupancy at 600 randomly chosen
protein-coding RefSeqSelect genes (excluding the network's own nodes)."""
import gzip, json, os, random, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
BASE="/path/to/revision"
BB=BASE+"/cache/direct/bigBedToBed"; URL="https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMap2022.bb"
OUT=BASE+"/cache/direct/remap/bg"; os.makedirs(OUT,exist_ok=True)
MAIN={"chr%s"%c for c in list(range(1,23))+["X"]}
net=set(json.load(open(BASE+"/cache/direct/node_tss_hg38.json")))
cand={}
with gzip.open(BASE+"/cache/direct/ncbiRefSeqSelect.txt.gz","rt") as f:
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[2] not in MAIN or not p[1].startswith("NM_"): continue
        if p[12] in net or p[12] in cand: continue
        s,e=int(p[4]),int(p[5])
        cand[p[12]]=(p[2], e if p[3]=="-" else s)
print("candidate background genes:",len(cand))
random.seed(20260909)
sel=random.sample(sorted(cand),600)
json.dump({g:cand[g] for g in sel},open(BASE+"/cache/direct/bg_tss.json","w"))
def q(g):
    chrom,t=cand[g]; path=os.path.join(OUT,g+".bed")
    if os.path.exists(path+".done"): return (g,"cached")
    for _ in range(3):
        try:
            r=subprocess.run([BB,"-chrom=%s"%chrom,"-start=%d"%max(0,t-5000),"-end=%d"%(t+5000),URL,path],
                             capture_output=True,timeout=600)
            if r.returncode==0:
                open(path+".done","w").write("1"); return (g,"ok")
        except Exception: pass
    return (g,"fail")
from collections import Counter
res=[]
with ThreadPoolExecutor(max_workers=8) as ex:
    for i,r in enumerate(ex.map(q,sel)):
        res.append(r)
        if (i+1)%100==0: print(" %d/600"%(i+1)); sys.stdout.flush()
print(Counter(x[1] for x in res))
