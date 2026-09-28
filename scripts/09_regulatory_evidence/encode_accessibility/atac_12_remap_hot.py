#!/usr/bin/env python3
"""(b, supporting) Genome-wide HOT-promoter index from ReMap 2022 non-redundant
hg38 (all TFs, all cell types): how many distinct TFs have a ChIP-seq peak in
each RefSeq-Select promoter (TSS +/-1 kb)?  Places the ETS1/NFKB1/RELA/SP1 peaks
at COL1A1 in context."""
import subprocess, sys, numpy as np, pandas as pd
ROOT="/path/to/revision"; C=f"{ROOT}/cache/v6/atac"; R=f"{ROOT}/results/v6"
tss=pd.read_csv(f"{C}/refseq_select_tss.tsv",sep="\t").drop_duplicates("symbol")
prom=tss.assign(start=(tss.tss-1000).clip(lower=0),end=tss.tss+1000)[["chrom","start","end","symbol"]].sort_values(["chrom","start"])
pb=f"{C}/prom1kb.bed"; prom.to_csv(pb,sep="\t",header=False,index=False)
pidx={s:i for i,s in enumerate(prom.symbol)}
M=np.zeros((len(prom),4000),dtype=bool); tfidx={}
cmd=(f"zcat {C}/remap/remap2022_nr_hg38.bed.gz | cut -f1-4 | "
     f"bedtools intersect -a stdin -b {pb} -wa -wb")
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True,bufsize=1<<20)
n=0
for line in p.stdout:
    f=line.rstrip("\n").split("\t")
    tf=f[3].split(":")[0]; sym=f[7]
    j=tfidx.get(tf)
    if j is None: j=len(tfidx); tfidx[tf]=j
    M[pidx[sym],j]=True; n+=1
    if n%20_000_000==0: print("lines",n,"TFs",len(tfidx),flush=True)
p.wait()
print("overlapping lines",n,"distinct TFs",len(tfidx),flush=True)
D=prom.copy(); D["n_distinct_TF"]=M[:,:len(tfidx)].sum(1)
D["pctile"]=D.n_distinct_TF.rank(pct=True)*100
for tf in ["ETS1","NFKB1","RELA","SP1"]:
    D[f"has_{tf}"]=M[:,tfidx[tf]] if tf in tfidx else False
D.to_csv(f"{R}/atac_remap_hot_promoters.csv",index=False)
pd.set_option("display.width",200)
print(D[D.symbol.isin(["COL1A1","COL3A1","COL1A2","GAPDH","ACTB","EPCAM","DCN","POSTN","FAP","LUM"])].to_string(index=False))
print("\nmedian distinct TFs per promoter:",D.n_distinct_TF.median(),"mean:",round(D.n_distinct_TF.mean(),1))
for tf in ["ETS1","NFKB1","RELA","SP1"]:
    print(f"{tf}: bound at {D[f'has_{tf}'].sum()}/{len(D)} promoters ({D[f'has_{tf}'].mean()*100:.1f}%)")
print("DONE")
