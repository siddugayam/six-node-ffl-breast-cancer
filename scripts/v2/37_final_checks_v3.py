#!/usr/bin/env python3
import csv, os, math, collections
import numpy as np
from scipy import stats
OUT="/path/to/revision/results/v2"
def rd(f):
    with open(os.path.join(OUT,f)) as fh: return list(csv.DictReader(fh))
def fl(v):
    try:
        x=float(v); return x if math.isfinite(x) else float('nan')
    except: return float('nan')

pr=rd("pancancer_ecological_product_v3.csv")
by=collections.defaultdict(dict)
for r in pr: by[(r["regulator"],r["target"])][r["cohort"]]=r
cohorts=sorted({r["cohort"] for r in pr})
print("PAIRED TEST ACROSS 31 COHORTS: |residual| (observed rho minus stroma-path product)")
print(f"{'regulator':12s}{'target':8s}{'mean_resid':>12s}{'mean_|resid|':>14s}{'wilcoxon p (resid=0)':>22s}")
resid={}
for (reg,tgt),d in by.items():
    v=np.array([fl(d[c]["residual"]) for c in cohorts if c in d])
    resid[(reg,tgt)]=v
    w=stats.wilcoxon(v)
    print(f"{reg:12s}{tgt:8s}{v.mean():>12.3f}{np.abs(v).mean():>14.3f}{w.pvalue:>22.2e}")
print("\nETS1 vs miR-29c residual (paired across cohorts):")
for tgt in ["COL1A1","COL3A1"]:
    a=np.abs(resid[("ETS1",tgt)]); b=np.abs(resid[("hsa-miR-29c",tgt)])
    w=stats.wilcoxon(a,b)
    print(f"  {tgt}: mean|resid| ETS1={a.mean():.3f}  miR-29c={b.mean():.3f}  paired Wilcoxon p={w.pvalue:.2e}")
print("\nMean rho(regulator, CAF) across 31 cohorts:")
for reg in ["ETS1","NFKB1","SP1","RELA","MYC","TP53","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-101"]:
    v=np.array([fl(by[(reg,"COL1A1")][c]["rho_reg_CAF"]) for c in cohorts if c in by[(reg,"COL1A1")]])
    print(f"  {reg:12s} mean {v.mean():+.3f}  median {np.median(v):+.3f}  range [{v.min():+.3f},{v.max():+.3f}]  "
          f"positive in {int((v>0).sum())}/{len(v)}")
