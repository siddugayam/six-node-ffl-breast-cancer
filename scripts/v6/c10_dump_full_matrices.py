#!/usr/bin/env python3
"""Dump collapsed mature-miRNA matrices + pheno for the cohorts new to v6."""
import sys, os, csv, numpy as np, re
sys.path.insert(0,"/path/to/revision/scripts/v6")
from c07_build_mirna_cohorts import read_matrix, chars_to_dict, load_map, MX
OUT="/path/to/revision/cache/v6/cohorts/built"

def collapse(probes, X, gplmap):
    """probe -> best hsa- mature name; keep highest-mean probe per name."""
    name={}
    for i,p in enumerate(probes):
        cands=([p]+ (gplmap.get(p,[]) if gplmap else []))
        pick=None
        for c in cands:
            c=c.strip()
            if c.lower().startswith("hsa-"): pick=c; break
        if pick: name[i]=pick
    by={}
    for i,nm in name.items():
        mu=np.nanmean(X[i,:]) if np.any(np.isfinite(X[i,:])) else -1e18
        if nm not in by or mu>by[nm][1]: by[nm]=(i,mu)
    rows=sorted(by.keys())
    M=np.vstack([X[by[r][0],:] for r in rows])
    return rows, M

JOBS=[("GSE59829","GSE59829_series_matrix.txt.gz","GPL8179"),
      ("GSE78870","GSE78870_series_matrix.txt.gz","GPL20662"),
      ("GSE57897","GSE57897_series_matrix.txt.gz","GPL18722"),
      ("GSE103161","GSE103161_series_matrix.txt.gz","GPL23960"),
      ("GSE97811","GSE97811_series_matrix.txt.gz","GPL21263"),
      ("GSE28321","GSE28321_series_matrix.txt.gz","GPL5106"),
      ("GSE40267","GSE40267_series_matrix.txt.gz","GPL10850"),
      ("GSE26666","GSE26666-GPL8227_series_matrix.txt.gz","GPL8227"),
      ("GSE73002","GSE73002_series_matrix.txt.gz","GPL18941")]
for cohort,f,gpl in JOBS:
    probes,gsm,X,chars,meta=read_matrix(os.path.join(MX,f))
    n=len(gsm); cd=chars_to_dict(chars,n); gm=load_map(gpl)
    rows,M=collapse(probes,X,gm)
    with open(os.path.join(OUT,cohort+"_expr.csv"),"w",newline="") as o:
        w=csv.writer(o); w.writerow(["miRNA"]+list(gsm))
        for i,r in enumerate(rows):
            w.writerow([r]+["" if not np.isfinite(v) else "%.6f"%v for v in M[i,:]])
    keys=[k for k in cd if k!="_free"]
    with open(os.path.join(OUT,cohort+"_pheno.csv"),"w",newline="") as o:
        w=csv.writer(o); w.writerow(["gsm","title","source"]+keys)
        for j in range(n):
            w.writerow([gsm[j], meta.get("title",[""]*n)[j], meta.get("source",[""]*n)[j]]
                       +[(cd[k][j] if cd[k][j] is not None else "") for k in keys])
    print("%-11s n=%-5d probes=%d -> mature=%d"%(cohort,n,len(probes),len(rows)), flush=True)
