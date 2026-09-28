#!/usr/bin/env python3
"""Dose-response: does the human collagen signal fall progressively with PDX passage,
as residual human stroma is replaced by mouse?"""
import json
import numpy as np, pandas as pd
from scipy import stats
BASE="/path/to/revision"; RES=f"{BASE}/results/v6"
P=pd.read_csv(f"{RES}/pdx_percentile_ranks_persample.csv", index_col=0)
P=P[P.passage.notna()].copy()
def pnum(p):
    if p=="Originator": return -1
    if isinstance(p,str) and p.startswith("P") and p[1:].isdigit(): return int(p[1:])
    return np.nan
P["pnum"]=P.passage.map(pnum)
D=P[P.pnum.notna() & (P.pnum<=5)].copy()
genes=["COL1A1","COL3A1","COL1A2","FN1","POSTN","PDGFRB","CXCL12","DCN","LUM","PTPRC",
       "ACTA2","THY1","EPCAM","KRT8","E2F1","EZH2","MYBL2","NFKB1","ETS1","GATA3"]
rows=[]
for g in genes:
    if g not in D.columns: continue
    sub=D[["pnum",g]].dropna()
    pdxonly=sub[sub.pnum>=0]
    r,p=stats.spearmanr(pdxonly.pnum, pdxonly[g])
    med=sub.groupby("pnum")[g].median()
    rows.append(dict(gene=g, n=len(sub),
        median_Originator=med.get(-1,np.nan), median_P0=med.get(0,np.nan),
        median_P1=med.get(1,np.nan), median_P2=med.get(2,np.nan),
        median_P3=med.get(3,np.nan), median_P4=med.get(4,np.nan), median_P5=med.get(5,np.nan),
        spearman_rho_passage_within_PDX=r, p_value=p, n_pdx=len(pdxonly)))
T=pd.DataFrame(rows)
T.to_csv(f"{RES}/pdx_passage_trend.csv", index=False)
print(T.round(3).to_string(index=False))
print("\nn per passage:", D.pnum.value_counts().sort_index().to_dict())
