#!/usr/bin/env python3
"""
(1) Spot-level mediation: TF -> stromal content -> COL1A1/COL3A1/FN1, run inside
    each Visium section and pooled, mirroring the bulk mediation analysis but
    with the mediator measured in situ rather than deconvolved.
(2) Summary of how much of the bulk TF-collagen correlation survives
    compartment stratification.
"""
import os
import numpy as np, pandas as pd
from scipy import stats

BASE="/path/to/revision"; RES=f"{BASE}/results/v6"
S=pd.read_csv(f"{RES}/spatial_spot_level_values.csv.gz", index_col=0)
OUT=["COL1A1","COL3A1","FN1"]
TFS=["NFKB1","ETS1","E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3"]
rng=np.random.default_rng(20260909)

def z(v):
    v=np.asarray(v,float); s=np.nanstd(v)
    return (v-np.nanmean(v))/s if s>0 else np.zeros_like(v)

def paths(x,m,y):
    a=np.polyfit(x,m,1)[0]
    A=np.column_stack([np.ones_like(x),x,m])
    beta=np.linalg.lstsq(A,y,rcond=None)[0]
    cprime=beta[1]; b=beta[2]
    c=np.polyfit(x,y,1)[0]
    return a,b,cprime,c

rows=[]
for sec,d in S.groupby("section"):
    if not set(OUT+TFS).issubset(d.columns): continue
    for o in OUT:
        for t in TFS:
            sub=d[[o,t,"stroma_score"]].dropna()
            if len(sub)<200: continue
            x=z(sub[t]); m=z(sub["stroma_score"]); y=z(sub[o])
            a,b,cp,c=paths(x,m,y)
            bs=np.empty(500)
            n=len(x)
            for i in range(500):
                idx=rng.integers(0,n,n)
                aa,bb,_,_=paths(x[idx],m[idx],y[idx]); bs[i]=aa*bb
            lo,hi=np.percentile(bs,[2.5,97.5])
            sa=np.std(x)  # standardised, so SEs from regression:
            # Sobel
            resm=m-(np.mean(m)+a*(x-np.mean(x)))
            sea=np.sqrt(np.sum(resm**2)/(n-2)/np.sum((x-np.mean(x))**2))
            A=np.column_stack([np.ones(n),x,m]); beta=np.linalg.lstsq(A,y,rcond=None)[0]
            resy=y-A@beta; s2=np.sum(resy**2)/(n-3)
            cov=s2*np.linalg.inv(A.T@A); seb=np.sqrt(cov[2,2])
            sob=a*b/np.sqrt(b*b*sea*sea+a*a*seb*seb)
            rows.append(dict(section=sec, tf=t, outcome=o, n_spots=n,
                total_c=c, a_path=a, b_path=b, direct_cprime=cp, mediated_ab=a*b,
                mediated_lo=lo, mediated_hi=hi,
                prop_mediated=(a*b)/c if abs(c)>1e-6 else np.nan,
                sobel_z=sob, sobel_p=2*stats.norm.sf(abs(sob))))
M=pd.DataFrame(rows)
M.to_csv(f"{RES}/spatial_spotlevel_mediation_bysection.csv", index=False)

agg=M.groupby(["outcome","tf"]).agg(n_sections=("section","size"),
      median_total_c=("total_c","median"), median_mediated_ab=("mediated_ab","median"),
      median_direct_cprime=("direct_cprime","median"),
      median_prop_mediated=("prop_mediated","median"),
      n_sections_mediated_positive=("mediated_ab",lambda v:int((v>0).sum())),
      n_sections_direct_opposite_sign=("direct_cprime",lambda v:0)).reset_index()
sgn=M.assign(opp=np.sign(M.direct_cprime)!=np.sign(M.mediated_ab)).groupby(["outcome","tf"]).opp.sum()
agg["n_sections_direct_opposite_sign"]=agg.set_index(["outcome","tf"]).index.map(sgn).astype(int)
agg.to_csv(f"{RES}/spatial_spotlevel_mediation_summary.csv", index=False)
print(agg[agg.outcome=="COL1A1"].round(3).to_string(index=False))

# ---- how much of the bulk correlation survives stratification --------------
meta=pd.read_csv(f"{RES}/spatial_tf_collagen_correlation_meta.csv")
srows=[]
for (sch,o),d in meta.groupby(["scheme","outcome"]):
    d=d.dropna(subset=["rho_tcga_bulk"])
    if len(d)<5: continue
    w,p=stats.wilcoxon(d.rho_stromal_meta.abs(), d.rho_tcga_bulk.abs())
    w2,p2=stats.wilcoxon(d.rho_stromal_meta.abs(), d.rho_all_meta.abs())
    srows.append(dict(scheme=sch, outcome=o, n_tfs=len(d),
        mean_abs_rho_tcga_bulk=d.rho_tcga_bulk.abs().mean(),
        mean_abs_rho_allspots=d.rho_all_meta.abs().mean(),
        mean_abs_rho_stromal=d.rho_stromal_meta.abs().mean(),
        mean_abs_rho_epithelial=d.rho_epithelial_meta.abs().mean(),
        pct_of_bulk_retained_in_stroma=100*d.rho_stromal_meta.abs().mean()/d.rho_tcga_bulk.abs().mean(),
        wilcoxon_p_stromal_vs_bulk=p, wilcoxon_p_stromal_vs_allspots=p2,
        n_tfs_sign_flip_vs_bulk=int(d.sign_flip_vs_tcga_bulk.sum())))
SS=pd.DataFrame(srows); SS.to_csv(f"{RES}/spatial_bulk_vs_compartment_summary.csv", index=False)
print(SS.round(4).to_string(index=False))
print("DONE")
