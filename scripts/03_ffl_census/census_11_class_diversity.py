#!/usr/bin/env python3
"""Edge-class diversity per module for n=3..7 from the class-mask histograms."""
import re, sys, json, os
from collections import Counter

FINE = ['Gene->Gene','Gene->TF','TF->Gene','TF->TF','TF->miRNA','miRNA->Gene','miRNA->TF','miRNA->miRNA']
# taxonomy maps: index -> group id
TAX = {
 "fine8 (source_type -> target_type)": {i:i for i in range(8)},
 "7-class (all Gene-> edges pooled as gene_gene)": {0:0,1:0,2:1,3:2,4:3,5:4,6:5,7:6},
 "6-class coarse edge_type (miRNA_target pools miRNA->Gene/TF)": {0:0,1:1,2:2,3:2,4:3,5:4,6:4,7:5},
 "5-class deposited edge_type (Gene->TF pooled into gene_gene)": {0:0,1:0,2:1,3:1,4:2,5:3,6:3,7:4},
}

def parse(path):
    """returns list of (K, seed, found, est, maskhist dict)"""
    out=[]
    cur=None
    for line in open(path):
        m=re.match(r"K=(\d+) q=([\d.,]*) seed=(\d+) found=(\d+) est=([\d.]+) visited=(\d+) maxclasses=(\d+)", line)
        if m:
            cur=dict(K=int(m.group(1)), q=m.group(2), seed=m.group(3), found=int(m.group(4)),
                     est=float(m.group(5)), visited=int(m.group(6)), maxcls=int(m.group(7)), mask={})
            out.append(cur); continue
        if line.startswith("maskhist:") and cur is not None:
            for tok in line.split()[1:]:
                k,v=tok.split(":"); cur["mask"][int(k)]=int(v)
    return out

def diversity(mask, tax):
    """max distinct groups, and count of modules attaining it"""
    best=0; cnt=0; hist=Counter()
    for mk,n in mask.items():
        g=set()
        for i in range(8):
            if mk>>i & 1: g.add(tax[i])
        d=len(g); hist[d]+=n
        if d>best: best=d; cnt=n
        elif d==best: cnt+=n
    return best,cnt,hist

if __name__=="__main__":
    logs="/path/to/revision/logs/v2"
    files={3:"census2_n3.log",4:"census2_n4.log",5:"census2_n5.log",
           6:"census_n6_sampled.log",7:"census_n7_sampled.log"}
    for tname,tax in TAX.items():
        print("=== taxonomy:", tname, " (%d possible classes)"%len(set(tax.values())))
        for n in (3,4,5,6,7):
            p=os.path.join(logs,files[n])
            if not os.path.exists(p): print("  n=%d  MISSING %s"%(n,p)); continue
            runs=parse(p)
            if not runs: print("  n=%d  no runs parsed"%n); continue
            allmask=Counter()
            for r in runs:
                for k,v in r["mask"].items(): allmask[k]+=v
            best,cnt,hist=diversity(allmask,tax)
            exhaustive = all(r["q"].strip(",")=="" or set(r["q"].split(","))<= {"1",""} for r in runs)
            print("  n=%d  %-11s max_classes=%d  modules_at_max(raw counts over %d run(s))=%d   hist=%s"
                  %(n,"exhaustive" if exhaustive else "sampled",best,len(runs),cnt,dict(sorted(hist.items()))))
        print()
