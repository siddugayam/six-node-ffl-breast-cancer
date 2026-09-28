#!/usr/bin/env python3
"""INDEPENDENT verification of PART 1A (DEMETER2 RNAi).
Deliberately does not import netlib or reuse the deposited R code path."""
import numpy as np, pandas as pd, sys, re
REV="/path/to/revision"
CA=f"{REV}/cache/v7"; OUT=f"{REV}/results/v3"
RAW="/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
rng=np.random.default_rng(20260910)

D=pd.read_csv(f"{CA}/D2_combined_gene_dep_scores.csv",index_col=0)
print("D2 matrix raw:",D.shape)
D.index=[re.sub(r" \(\d+\)$","",g) for g in D.index]
print("duplicated gene symbols in D2:",pd.Index(D.index).duplicated().sum())
bcl=[c for c in D.columns if c.endswith("_BREAST")]
print("breast lines (col suffix _BREAST):",len(bcl))
si=pd.read_csv(f"{CA}/D2_sample_info.csv")
sb=si[si.CCLE_ID.isin(bcl)]
print("  in_Achilles",int(sb.in_Achilles.sum()),"in_DRIVE",int(sb.in_DRIVE.sum()),"in_Marcotte",int(sb.in_Marcotte.sum()))
# cross-check via sample_info disease field rather than the name suffix
sb2=si[si.disease=="breast"]
print("  sample_info disease=='breast':",len(sb2),"| set equal to suffix method:",set(sb2.CCLE_ID)==set(bcl))

M=D[bcl]
d2b=M.mean(axis=1,skipna=True)
d2o=D[[c for c in D.columns if c not in bcl]].mean(axis=1,skipna=True)
d2fe=(M<-0.5).sum(axis=1)/M.notna().sum(axis=1)
d2n=M.notna().sum(axis=1)
print("genome-wide mean D2 in breast lines: %.4f | frac genes mean<-0.5: %.4f"%(d2b.mean(),(d2b<-0.5).mean()))
for g in ["RPL7","PSMA1","POLR2A","EIF3B","KRAS","GAPDH","PTEN","CDKN2A"]:
    if g in d2b.index: print("   control %-8s D2=%.3f n=%d"%(g,d2b[g],d2n[g]))

# CRISPR
mod=pd.read_csv(f"{REV}/data/depmap/Model_24Q4.csv")
bid=set(mod.loc[mod.OncotreeLineage=="Breast","ModelID"])
ce=pd.read_csv(f"{RAW}/CRISPRGeneEffect.csv",index_col=0)
ce.columns=[re.sub(r" \(\d+\)$","",c) for c in ce.columns]
print("CRISPR dup gene cols:",pd.Index(ce.columns).duplicated().sum())
CB=ce.loc[ce.index.isin(bid)]
print("CRISPR breast lines:",CB.shape[0],"genes:",CB.shape[1])
chb=CB.mean(axis=0,skipna=True); chfe=(CB<-0.5).sum(axis=0)/CB.notna().sum(axis=0)

FOC=["NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","MYC"]
rows=[]
for g in FOC:
    rows.append(dict(gene=g,
        d2_breast=d2b.get(g,np.nan), d2_n_lines=int(d2n.get(g,0)),
        d2_frac_lines_essential=d2fe.get(g,np.nan),
        d2_percentile_breast=float((d2b<=d2b.get(g,np.nan)).mean()) if g in d2b.index else np.nan,
        chronos_breast=chb.get(g,np.nan),
        chronos_percentile_breast=float((chb<=chb.get(g,np.nan)).mean()) if g in chb.index else np.nan))
F=pd.DataFrame(rows)
print("\n--- FOCUS GENES (independent) ---"); print(F.to_string(index=False,float_format=lambda x:"%.4f"%x))

common=[g for g in d2b.index if g in chb.index]
from scipy.stats import spearmanr,pearsonr
a=d2b[common].values; b=chb[common].values; ok=np.isfinite(a)&np.isfinite(b)
print("\nCRISPR-vs-RNAi over %d shared genes: Spearman %.4f Pearson %.4f"%(ok.sum(),spearmanr(a[ok],b[ok]).statistic,pearsonr(a[ok],b[ok]).statistic))
F.to_csv(f"{OUT}/screens_verify_demeter2_focus.csv",index=False)
np.save("/path/to/scratch/d2b.npy",d2b.values)
d2b.to_csv("/path/to/scratch/d2_breast_means.csv")
chb.to_csv("/path/to/scratch/chronos_breast_means.csv")
