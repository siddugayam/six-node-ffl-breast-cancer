#!/usr/bin/env python3
"""(a) Locus-level test: does aggregate ATAC accessibility at a gene locus track
RNA-derived stromal content in the same 74 TCGA-BRCA patients?

Locus score = mean of per-peak z-scored accessibility over peaks within a window
of the TSS.  Null = 10,000 random sets of peaks matched on peak number, GC decile
and accessibility decile.  BH across genes."""
import numpy as np, pandas as pd
from scipy import stats
ROOT="/path/to/revision"; C=f"{ROOT}/cache/v6/atac"; R=f"{ROOT}/results/v6"
rng=np.random.default_rng(11)
meta=pd.read_csv(f"{C}/peak_stroma_correlation.tsv.gz",sep="\t")
S=np.load(f"{C}/brca_sample_matrix.npy"); samp=pd.read_csv(f"{C}/brca_samples.tsv",sep="\t")
tss=pd.read_csv(f"{C}/refseq_select_tss.tsv",sep="\t").drop_duplicates("symbol").set_index("symbol")
sc=pd.read_csv(f"{ROOT}/results/v5/farmer_tcga_persample_scores.csv")
pat=samp.tcga_patient.values; upat=pd.unique(pat)
P=np.column_stack([S[:,pat==p].mean(1) for p in upat])
sc["patient"]=sc["sample"].str.slice(0,12)
scm=sc[sc.patient.isin(upat)].drop_duplicates("patient").set_index("patient").reindex(upat)
STR=scm["meanZ_ESTIMATE_STROMAL_noCOL"].values; CAF=scm["meanZ_CAF_scRNA_50"].values
Z=(P-P.mean(1,keepdims=True))/(P.std(1,keepdims=True)+1e-9)
chrom=meta.seqnames.values; mid=((meta.start+meta.end)//2).values
gcq=pd.qcut(meta.percentGC,10,labels=False,duplicates="drop").values
acq=pd.qcut(meta.mean_acc,10,labels=False,duplicates="drop").values
strat=gcq*10+acq
GENES=["COL1A1","COL3A1","COL1A2","COL5A1","COL6A3","FN1","DCN","LUM","POSTN","FAP","PDGFRB","THY1",
       "ACTA2","VIM","EPCAM","KRT8","KRT18","ESR1","GATA3","FOXA1","PTPRC","CD3E","GAPDH","ACTB",
       "ETS1","NFKB1","RELA","SP1"]
def idx(g,win):
    if g not in tss.index: return np.array([],dtype=int)
    t=tss.loc[g]; return np.where((chrom==t.chrom)&(np.abs(mid-t.tss)<=win))[0]
rows=[]
for win,tag in [(5000,"prom5kb"),(100000,"locus100kb")]:
    for g in GENES:
        i=idx(g,win)
        if len(i)==0:
            rows.append(dict(gene=g,window=tag,n_peaks=0)); continue
        v=Z[i].mean(0)
        for sname,sv in [("stromal_noCOL",STR),("CAF_scRNA50",CAF)]:
            rho,p=stats.spearmanr(v,sv)
            null=np.zeros(10000)
            pools={s:np.where(strat==s)[0] for s in np.unique(strat[i])}
            rv=stats.rankdata(sv); rv=(rv-rv.mean())/rv.std()
            for b in range(10000):
                pick=np.concatenate([rng.choice(pools[s],int((strat[i]==s).sum())) for s in np.unique(strat[i])])
                vv=Z[pick].mean(0); rr=stats.rankdata(vv); rr=(rr-rr.mean())/rr.std()
                null[b]=(rr@rv)/len(rv)
            phi=(null>=rho).mean(); plo=(null<=rho).mean()
            rows.append(dict(gene=g,window=tag,n_peaks=len(i),score=sname,rho=rho,p_param=p,
                             null_mean=null.mean(),null_sd=null.std(),
                             p_emp=min(1.0,2*min(phi,plo)),pctile_vs_null=(null<rho).mean()*100))
        print(f"{tag} {g} n={len(i)}",flush=True)
D=pd.DataFrame(rows)
m=D.p_emp.notna(); p=D.loc[m,"p_emp"].values; o=np.argsort(p); q=np.empty(len(p))
q[o]=np.minimum.accumulate((p[o]*len(p)/(np.arange(len(p))+1))[::-1])[::-1]
D.loc[m,"fdr_BH"]=np.clip(q,0,1)
D.to_csv(f"{R}/atac_locus_stroma_tracking.csv",index=False)
pd.set_option("display.width",220)
print(D[D.score=="stromal_noCOL"].sort_values(["window","rho"],ascending=[True,False]).to_string(index=False))
print("DONE")
