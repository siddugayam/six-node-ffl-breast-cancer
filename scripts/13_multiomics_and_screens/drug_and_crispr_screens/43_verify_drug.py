#!/usr/bin/env python3
"""Independent verification of PART 1C (PRISM + GDSC2 drug-sensitivity associations)."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
REV="/path/to/revision"; OUT=f"{REV}/results/v3"; CA=f"{REV}/cache/v7"
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=np.empty(n)
    r[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]
    return np.minimum(r,1)
S=pd.read_csv(f"{OUT}/screens_celline_module_scores.csv")
scores=["COLLAGEN","MIR29_TARGET","MIR130A_TARGET","MIR29_STRONG","MIR130A_STRONG"]
print("breast lines with module scores:",len(S))
# --- independent re-derivation of the COLLAGEN score from DepMap expression ---
E=pd.read_csv(f"{REV}/data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1_24Q4.csv",index_col=0)
E.columns=[c.split(" (")[0] for c in E.columns]
mod=pd.read_csv(f"{REV}/data/depmap/Model_24Q4.csv")
bre=set(mod.loc[mod.OncotreeLineage=="Breast","ModelID"])
Eb=E.loc[E.index.isin(bre)]
keep=(Eb>0.5).mean(axis=0)>=0.25
Xb=Eb.loc[:,keep]
Z=(Xb-Xb.mean(axis=0))/Xb.std(axis=0,ddof=1)
COL=["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL11A1"]
g=[x for x in COL if x in Z.columns]
mycol=Z[g].mean(axis=1)
chk=S.set_index("ModelID").reindex(mycol.index)["COLLAGEN"]
print("collagen genes retained:",len(g),"of",len(COL),"| max |deposited - reindependent| = %.3e"%np.nanmax(np.abs(mycol-chk)))
# --- GDSC2 ---
G=pd.read_excel(f"{CA}/GDSC2_fitted_dose_response.xlsx")
print("GDSC2 rows:",G.shape)
map_=mod[["SangerModelID","ModelID"]].dropna(); map_=map_[map_.SangerModelID!=""]
G=G.merge(map_,left_on="SANGER_MODEL_ID",right_on="SangerModelID")
Gb=G[G.ModelID.isin(set(S.ModelID))].copy()
print("GDSC2 breast rows:",len(Gb),"lines",Gb.ModelID.nunique(),"drugs",Gb.DRUG_ID.nunique())
Gb=Gb[np.isfinite(Gb.LN_IC50)]
agg=Gb.groupby(["DRUG_ID","ModelID"]).agg(LN_IC50=("LN_IC50","median"),DRUG_NAME=("DRUG_NAME","first"),
      PUTATIVE_TARGET=("PUTATIVE_TARGET","first"),PATHWAY_NAME=("PATHWAY_NAME","first")).reset_index()
cnt=agg.groupby("DRUG_ID").size(); keepd=set(cnt[cnt>=15].index)
agg=agg[agg.DRUG_ID.isin(keepd)]
print("GDSC2 drugs in >=15 breast lines:",len(keepd),"of",len(cnt))
W=agg.pivot(index="DRUG_ID",columns="ModelID",values="LN_IC50")
meta=agg.groupby("DRUG_ID")[["DRUG_NAME","PUTATIVE_TARGET","PATHWAY_NAME"]].first()
rows=[]
for sc in scores:
    x=S.set_index("ModelID")[sc]
    for did,y in W.iterrows():
        common=[c for c in W.columns if c in x.index]
        yy=y[common].values.astype(float); xx=x[common].values.astype(float)
        ok=np.isfinite(yy)&np.isfinite(xx)
        if ok.sum()<15: continue
        r=spearmanr(xx[ok],yy[ok])
        rows.append(dict(score=sc,DRUG_ID=did,n_lines=int(ok.sum()),rho=r.statistic,p=r.pvalue))
GR=pd.DataFrame(rows)
GR["fdr_within_score"]=GR.groupby("score")["p"].transform(lambda v: bh(v.values))
GR["fdr_global"]=bh(GR.p.values)
GR=GR.merge(meta,on="DRUG_ID",how="left").sort_values("p")
print("\nGDSC2 tests:",len(GR))
print(GR.groupby("score").agg(n=("p","size"),n_fdr05=("fdr_within_score",lambda v:(v<0.05).sum()),min_fdr=("fdr_within_score","min")).to_string())
print("\nGDSC2 survivors at within-score BH<0.05:")
print(GR[GR.fdr_within_score<0.05].to_string(index=False))
print("\nGDSC2 survivors at GLOBAL BH<0.05 across all %d tests: %d"%(len(GR),(GR.fdr_global<0.05).sum()))
print(GR.head(5)[["score","DRUG_NAME","PUTATIVE_TARGET","n_lines","rho","p","fdr_within_score","fdr_global"]].to_string(index=False))
GR.to_csv(f"{OUT}/screens_verify_gdsc2_associations.csv",index=False)
dep=pd.read_csv(f"{OUT}/screens_gdsc2_associations.csv")
m=GR.merge(dep[["score","id","rho","p","fdr"]],left_on=["score","DRUG_ID"],right_on=["score","id"],suffixes=("_new","_dep"))
print("\nmerged with deposited: %d rows | max |d rho| = %.2e | max |d p| = %.2e | max |d fdr| = %.2e"%(
   len(m),np.nanmax(np.abs(m.rho_new-m.rho_dep)),np.nanmax(np.abs(m.p_new-m.p_dep)),np.nanmax(np.abs(m.fdr_within_score-m.fdr))))
