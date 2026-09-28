#!/usr/bin/env python3
import numpy as np, pandas as pd
from scipy.stats import spearmanr
REV="/path/to/revision"; OUT=f"{REV}/results/v3"; CA=f"{REV}/cache/v7"
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=np.empty(n)
    r[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]; return np.minimum(r,1)
S=pd.read_csv(f"{OUT}/screens_celline_module_scores.csv"); scores=["COLLAGEN","MIR29_TARGET","MIR130A_TARGET","MIR29_STRONG","MIR130A_STRONG"]
bset=set(S.ModelID)
chunks=[]
for ch in pd.read_csv(f"{CA}/prism_secondary_dose_response.csv",usecols=["broad_id","depmap_id","auc","name","moa","target","phase"],chunksize=1_000_000):
    chunks.append(ch[ch.depmap_id.isin(bset)])
P=pd.concat(chunks); print("PRISM secondary rows in breast lines with a score:",len(P),"| lines",P.depmap_id.nunique(),"| compounds",P.broad_id.nunique())
P=P[np.isfinite(P.auc)]
agg=P.groupby(["broad_id","depmap_id"]).agg(auc=("auc","median"),name=("name","first"),moa=("moa","first"),target=("target","first"),phase=("phase","first")).reset_index()
cnt=agg.groupby("broad_id").size(); keep=set(cnt[cnt>=15].index)
print("compounds in >=15 breast lines:",len(keep),"of",len(cnt))
agg=agg[agg.broad_id.isin(keep)]
W=agg.pivot(index="broad_id",columns="depmap_id",values="auc")
meta=agg.groupby("broad_id")[["name","moa","target","phase"]].first()
rows=[]
for sc in scores:
    x=S.set_index("ModelID")[sc]
    common=[c for c in W.columns if c in x.index]; xx=x[common].values.astype(float)
    sub=W[common].values.astype(float)
    for i,bid in enumerate(W.index):
        yy=sub[i]; ok=np.isfinite(yy)&np.isfinite(xx)
        if ok.sum()<15: continue
        r=spearmanr(xx[ok],yy[ok]); rows.append(dict(score=sc,broad_id=bid,n_lines=int(ok.sum()),rho=r.statistic,p=r.pvalue))
PR=pd.DataFrame(rows)
PR["fdr_within_score"]=PR.groupby("score")["p"].transform(lambda v: bh(v.values)); PR["fdr_global"]=bh(PR.p.values)
PR=PR.merge(meta,on="broad_id",how="left").sort_values("p")
print("\nPRISM tests:",len(PR))
print(PR.groupby("score").agg(n=("p","size"),n_fdr05=("fdr_within_score",lambda v:(v<0.05).sum()),min_fdr=("fdr_within_score","min")).to_string())
print("\nglobal-BH survivors:",int((PR.fdr_global<0.05).sum()))
print(PR.head(6)[["score","name","moa","n_lines","rho","p","fdr_within_score","fdr_global"]].to_string(index=False))
PR.to_csv(f"{OUT}/screens_verify_prism_associations.csv",index=False)
dep=pd.read_csv(f"{OUT}/screens_prism_associations.csv")
m=PR.merge(dep[["score","id","rho","p","fdr"]],left_on=["score","broad_id"],right_on=["score","id"],suffixes=("_new","_dep"))
print("\nmerged with deposited: %d of %d rows | max|d rho|=%.2e max|d fdr|=%.2e"%(len(m),len(PR),np.nanmax(np.abs(m.rho_new-m.rho_dep)),np.nanmax(np.abs(m.fdr_within_score-m.fdr))))
