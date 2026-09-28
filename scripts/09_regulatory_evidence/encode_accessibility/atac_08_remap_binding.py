#!/usr/bin/env python3
"""(b)/(c) Measured genome-wide TF occupancy: ReMap 2022 non-redundant hg38 peaks
for ETS1, NFKB1, RELA, SP1 at promoters (TSS +/-1 kb).

Null 1: all RefSeq-Select promoters.
Null 2: promoters that carry a TCGA-BRCA ATAC peak (i.e. accessible in breast
        tumours) - the conditioning that makes the comparison fair.
"""
import os, subprocess, tempfile, numpy as np, pandas as pd
from scipy.stats import fisher_exact
ROOT="/path/to/revision"; C=f"{ROOT}/cache/v6/atac"; R=f"{ROOT}/results/v6"
TFS=["ETS1","NFKB1","RELA","SP1"]
tss=pd.read_csv(f"{C}/refseq_select_tss.tsv",sep="\t").drop_duplicates("symbol")
prom=tss.assign(start=(tss.tss-1000).clip(lower=0),end=tss.tss+1000)[["chrom","start","end","symbol"]]
prom=prom.sort_values(["chrom","start"])
tmp=tempfile.mkdtemp(); pb=f"{tmp}/prom.bed"; prom.to_csv(pb,sep="\t",header=False,index=False)
B=prom[["symbol"]].copy().set_index("symbol")
for tf in TFS:
    f=f"{C}/remap/{tf}_nr.bed.gz"
    out=subprocess.run(f"zcat {f} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > {tmp}/{tf}.bed",shell=True)
    r=subprocess.run(["bedtools","intersect","-a",pb,"-b",f"{tmp}/{tf}.bed","-c"],capture_output=True,text=True).stdout
    d={l.split('\t')[3]:int(l.split('\t')[4]) for l in r.strip().split('\n')}
    B[f"remap_{tf}"]=(pd.Series(d).reindex(B.index)>0).astype(int)
    print(tf,"promoters bound:",B[f"remap_{tf}"].sum(),f"({B[f'remap_{tf}'].mean()*100:.1f}%)")
# accessible-in-BRCA flag
brca=pd.read_csv(f"{C}/brca_peaks.bed",sep="\t",header=None,names=["chrom","start","end","name"]).sort_values(["chrom","start"])
bb=f"{tmp}/brca.bed"; brca.to_csv(bb,sep="\t",header=False,index=False)
r=subprocess.run(["bedtools","intersect","-a",pb,"-b",bb,"-c"],capture_output=True,text=True).stdout
B["brca_atac_promoter_peak"]=(pd.Series({l.split('\t')[3]:int(l.split('\t')[4]) for l in r.strip().split('\n')}).reindex(B.index)>0).astype(int)
print("promoters with a TCGA-BRCA ATAC peak:",B.brca_atac_promoter_peak.sum(),f"({B.brca_atac_promoter_peak.mean()*100:.1f}%)")
B.to_csv(f"{R}/atac_remap_promoter_binding_allgenes.csv")

pc=pd.read_csv(f"{R}/atac_tf_target_compartment_adjusted.csv")
rows=[]
for tf in TFS:
    tg=[g for g in pc.loc[(pc.TF==tf)&(pc.poscontrol),"target"] if g in B.index]
    col=[g for g in ["COL1A1","COL3A1"] if g in B.index]
    for label,genes in [("positive_control_targets",tg),("COL1A1",["COL1A1"]),("COL3A1",["COL3A1"])]:
        if not genes: continue
        k=int(B.loc[genes,f"remap_{tf}"].sum()); n=len(genes)
        for bg_name,bgset in [("all_promoters",B),("accessible_promoters",B[B.brca_atac_promoter_peak==1])]:
            bg=bgset.drop(index=[g for g in genes if g in bgset.index],errors="ignore")
            K=int(bg[f"remap_{tf}"].sum()); N=len(bg)
            odds,p=fisher_exact([[k,n-k],[K,N-K]])
            rows.append(dict(TF=tf,set=label,n_genes=n,n_bound=k,frac_bound=k/n,background=bg_name,
                             bg_frac=K/N,odds_ratio=odds,p=p))
D=pd.DataFrame(rows)
o=np.argsort(D.p.values); q=np.empty(len(D)); m=len(D)
q[o]=np.minimum.accumulate((D.p.values[o]*m/(np.arange(m)+1))[::-1])[::-1]
D["fdr_BH"]=np.clip(q,0,1)
D.to_csv(f"{R}/atac_remap_binding_tests.csv",index=False)
pd.set_option("display.width",220)
print(D.to_string(index=False))
