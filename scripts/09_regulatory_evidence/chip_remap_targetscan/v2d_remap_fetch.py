#!/usr/bin/env python3
"""Query ReMap 2022 (hg38, all non-merged datasets) for every ChIP-seq peak whose
summit falls within TSS +/-5 kb of each gene/TF node of the canonical network.
Uses UCSC bigBedToBed against the remote reMap2022.bb (region-restricted reads;
no bulk download of the 1.45 GB catalogue)."""
import json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
BASE="/path/to/revision"
BB=BASE+"/cache/direct/bigBedToBed"
URL="https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMap2022.bb"
OUT=BASE+"/cache/direct/remap/regions"
os.makedirs(OUT,exist_ok=True)
WIN=5000
tss=json.load(open(BASE+"/cache/direct/node_tss_hg38.json"))
print("genes to query:",len(tss)); sys.stdout.flush()

def q(item):
    g,(chrom,t,strand,tx)=item
    path=os.path.join(OUT,g+".bed")
    if os.path.exists(path) and os.path.getsize(path)>=0 and os.path.exists(path+".done"):
        return (g,"cached",sum(1 for _ in open(path)))
    s=max(0,t-WIN); e=t+WIN
    for attempt in range(3):
        try:
            r=subprocess.run([BB,"-chrom=%s"%chrom,"-start=%d"%s,"-end=%d"%e,URL,path],
                             capture_output=True,timeout=600)
            if r.returncode==0:
                open(path+".done","w").write("%s\t%d\t%d\n"%(chrom,s,e))
                return (g,"ok",sum(1 for _ in open(path)))
        except Exception as ex:
            pass
    return (g,"fail",0)

res=[]
with ThreadPoolExecutor(max_workers=8) as ex:
    for i,r in enumerate(ex.map(q,sorted(tss.items()))):
        res.append(r)
        if (i+1)%40==0: print("  %d/%d"%(i+1,len(tss))); sys.stdout.flush()
from collections import Counter
print("status:",Counter(x[1] for x in res))
print("failed:",[x[0] for x in res if x[1]=="fail"])
print("total peak lines:",sum(x[2] for x in res))
