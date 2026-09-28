#!/usr/bin/env python3
"""Project SCORE: handle replicate screens of the same model explicitly."""
import numpy as np, pandas as pd
REV="/path/to/revision"; OUT=f"{REV}/results/v3"; CA=f"{REV}/cache/v7"
P=f"{CA}/ps2/Project_score_combined_Sanger_v2_Broad_21Q2_fitness_scores_scaled_bayesian_factors_20250624.tsv"
B=f"{CA}/ps2/Project_score_combined_Sanger_v2_Broad_21Q2_fitness_scores_binary_matrix_20250624.tsv"
hdr=pd.read_csv(P,sep="\t",nrows=5,header=None,dtype=str)
mid=hdr.iloc[1].tolist()[3:]; qcv=np.array(hdr.iloc[3].tolist()[3:]); srcv=np.array(hdr.iloc[2].tolist()[3:]); mnm=np.array(hdr.iloc[0].tolist()[3:])
M=pd.read_csv(P,sep="\t",skiprows=5,header=None); sym=M.iloc[:,1].astype(str)
X=M.iloc[:,3:].astype(float); X.index=sym; X.columns=range(X.shape[1])   # positional columns, no name collisions
cmp=pd.read_csv(f"{CA}/cmp_model_list.csv",low_memory=False)
bre=set(cmp.loc[cmp.cancer_type.astype(str).str.contains("Breast",case=False,na=False),"model_id"])
midA=np.array(mid)
sel=np.where(np.isin(midA,list(bre)) & (qcv=="TRUE"))[0]
print("breast QC-passing SCREEN columns:",len(sel),"| unique models:",len(set(midA[sel])))
print("source split of those screens:",pd.Series(srcv[sel]).value_counts().to_dict())
dupmods=pd.Series(midA[sel]).value_counts(); print("models screened twice:",int((dupmods>1).sum()),"| max screens per model:",int(dupmods.max()))
XB=X.iloc[:,sel]; XB.columns=midA[sel]
per_model=XB.T.groupby(level=0).mean().T          # average replicates within a model
scb=per_model.mean(axis=1)                        # then across the 43 models
print("per-model matrix:",per_model.shape)
# comparison: first-occurrence-only (what the deposited R intersect+.. selection does)
first_idx=[np.where(midA==m)[0][0] for m in sorted(set(midA[sel]))]
Xfirst=X.iloc[:,first_idx]; Xfirst.columns=sorted(set(midA[sel])); scb_first=Xfirst.mean(axis=1)
# binary
bh=pd.read_csv(B,sep="\t",nrows=5,header=None,dtype=str); bmid=np.array(bh.iloc[1].tolist()[3:])
BM=pd.read_csv(B,sep="\t",skiprows=5,header=None); BX=BM.iloc[:,3:]; BX.index=BM.iloc[:,1].astype(str); BX.columns=range(BX.shape[1])
bsel=np.where(np.isin(bmid,list(set(midA[sel]))))[0]
BB=BX.iloc[:,bsel]; BB.columns=bmid[bsel]
bin_per_model=(BB==1).T.groupby(level=0).max().T   # significant in either replicate
fsig=bin_per_model.mean(axis=1)
bfirst_idx=[np.where(bmid==m)[0][0] for m in sorted(set(midA[sel]))]
BF=(BX.iloc[:,bfirst_idx]==1); BF.columns=sorted(set(midA[sel])); fsig_first=BF.mean(axis=1)
dep=pd.read_csv(f"{OUT}/screens_projectscore_focus_genes.csv")
print("\ndeposited focus columns:",list(dep.columns))
FOC=["COL1A1","COL3A1","ESR1","ETS1","HIF1A","MKL1","MYC","NFKB1","RELA","SP1","STAT3","TP53"]
alias={"MKL1":"MRTFA"}
rows=[]
for g in FOC:
    k=g if g in scb.index else alias.get(g,g)
    rows.append(dict(gene=g,BF_perModelAvg=scb.get(k,np.nan),BF_firstScreenOnly=scb_first.get(k,np.nan),
                     fracSig_perModel=fsig.get(k,np.nan),fracSig_first=fsig_first.get(k,np.nan)))
F=pd.DataFrame(rows)
dep=dep.rename(columns={"node":"gene"});F=F.merge(dep[["gene","score_breast","frac_significant","n_lines"]],on="gene",how="left")
print("\n--- PROJECT SCORE FOCUS: replicate handling ---")
print(F.to_string(index=False,float_format=lambda x:"%.4f"%x))
F.to_csv(f"{OUT}/screens_verify_projectscore_focus.csv",index=False)
scb.to_csv("/path/to/scratch/ps_breast_permodel.csv")
scb_first.to_csv("/path/to/scratch/ps_breast_first.csv")
from scipy.stats import spearmanr
ok=scb.notna()&scb_first.notna()
print("\nagreement between the two replicate policies over %d genes: Spearman %.4f, Pearson %.4f"%(ok.sum(),spearmanr(scb[ok],scb_first[ok]).statistic,np.corrcoef(scb[ok],scb_first[ok])[0,1]))
