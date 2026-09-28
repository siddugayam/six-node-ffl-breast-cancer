#!/usr/bin/env python3
"""Independent verification of PART 1B (Project SCORE / Sanger)."""
import numpy as np, pandas as pd
REV="/path/to/revision"; OUT=f"{REV}/results/v3"; CA=f"{REV}/cache/v7"
P=f"{CA}/ps2/Project_score_combined_Sanger_v2_Broad_21Q2_fitness_scores_scaled_bayesian_factors_20250624.tsv"
B=f"{CA}/ps2/Project_score_combined_Sanger_v2_Broad_21Q2_fitness_scores_binary_matrix_20250624.tsv"
hdr=pd.read_csv(P,sep="\t",nrows=5,header=None,dtype=str)
print("header rows:"); print(hdr.iloc[:,:6].to_string())
model_name=hdr.iloc[0].tolist(); model_id=hdr.iloc[1].tolist(); src=hdr.iloc[2].tolist(); qc=hdr.iloc[3].tolist()
M=pd.read_csv(P,sep="\t",skiprows=5,header=None)
print("matrix:",M.shape)
sym=M.iloc[:,1].astype(str)
mid=model_id[3:]; qcv=qc[3:]; srcv=src[3:]; mnm=model_name[3:]
assert len(mid)==M.shape[1]-3, (len(mid),M.shape)
X=M.iloc[:,3:]; X.columns=mid; X.index=sym
print("dup symbols:",X.index.duplicated().sum()); X=X[~X.index.duplicated()]
print("sources:",pd.Series(srcv).value_counts().to_dict(),"| qc TRUE:",sum(np.array(qcv)=="TRUE"))
print("frac index starting with a letter: %.4f"%np.mean([bool(s[:1].isalpha()) for s in X.index]))
cmp=pd.read_csv(f"{CA}/cmp_model_list.csv",low_memory=False)
bre=set(cmp.loc[cmp.cancer_type.astype(str).str.contains("Breast",case=False,na=False),"model_id"])
qcok=set(np.array(mid)[np.array(qcv)=="TRUE"])
bcols=[c for c in X.columns if c in bre and c in qcok]
print("CMP breast models:",len(bre),"| breast cols passing QC:",len(bcols))
print("breast model names:",", ".join(sorted(set(np.array(mnm)[[mid.index(c) for c in bcols]]))))
XB=X[bcols].astype(float)
scb=XB.mean(axis=1)
print("\nvalue range %.3f .. %.3f | median %.4f"%(np.nanmin(XB.values),np.nanmax(XB.values),np.nanmedian(XB.values)))
print("SIGN CHECK (positive should = fitness gene):")
for g in ["RPL7","PSMA1","POLR2A","EIF3B","RPL5","MYC","PTEN","CDKN2A","KRAS","GAPDH"]:
    if g in scb.index: print("   %-8s mean scaled BF = %+8.3f"%(g,scb[g]))
# binary
bh=pd.read_csv(B,sep="\t",nrows=5,header=None,dtype=str); bmid=bh.iloc[1].tolist()[3:]
BM=pd.read_csv(B,sep="\t",skiprows=5,header=None); bsym=BM.iloc[:,1].astype(str)
BX=BM.iloc[:,3:]; BX.columns=bmid; BX.index=bsym; BX=BX[~BX.index.duplicated()]
bcb=[c for c in bcols if c in BX.columns]
print("binary breast cols:",len(bcb))
fsig=(BX[bcb]==1).sum(axis=1)/BX[bcb].notna().sum(axis=1)
print("genome-wide mean frac models a gene is a significant fitness gene: %.4f"%fsig.mean())
FOC=["COL1A1","COL3A1","ESR1","ETS1","HIF1A","MKL1","MYC","NFKB1","RELA","SP1","STAT3","TP53"]
alias={"MKL1":"MRTFA"}
rows=[]
for g in FOC:
    k=g if g in scb.index else alias.get(g,g)
    rows.append(dict(gene=g,used_symbol=k,mean_scaledBF_breast=scb.get(k,np.nan),
                     frac_models_significant=fsig.get(k,np.nan),n_models=int(XB.loc[k].notna().sum()) if k in XB.index else 0))
F=pd.DataFrame(rows); print("\n--- PROJECT SCORE FOCUS GENES (independent) ---")
print(F.to_string(index=False,float_format=lambda x:"%.4f"%x))
F.to_csv(f"{OUT}/screens_verify_projectscore_focus.csv",index=False)
scb.to_csv("/path/to/scratch/ps_breast_means.csv")
