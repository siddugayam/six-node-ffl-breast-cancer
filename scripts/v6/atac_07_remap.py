#!/usr/bin/env python3
"""(b, supporting) Measured TF occupancy (ReMap 2022, hg38, all cell types) at the
COL1A1 and COL3A1 promoters, placed against the local background of all RefSeq
promoters inside the same cached ReMap intervals.  ReMap peaks were cached for
five loci only (chr1:207.50-208.12, chr2:188.67-189.27, chr7:130.58-131.18,
chr11:57.34-57.94, chr17:49.90-50.50 Mb), so the background is those loci."""
import pandas as pd, numpy as np
ROOT="/path/to/revision"; C=f"{ROOT}/cache/v6/atac"; R=f"{ROOT}/results/v6"
rm=pd.read_csv(f"{ROOT}/cache/seqreg/remap2022_hg38_loci.bed",sep="\t",header=None,
               names=["chrom","start","end","id","score","strand","ts","te","rgb","tf","biotype"])
cov=rm.groupby("chrom").agg(lo=("start","min"),hi=("end","max"))
tss=pd.read_csv(f"{C}/refseq_select_tss.tsv",sep="\t").drop_duplicates("symbol")
tss=tss[tss.chrom.isin(cov.index)]
tss=tss[[cov.loc[r.chrom,"lo"]+2000 < r.tss < cov.loc[r.chrom,"hi"]-2000 for r in tss.itertuples()]]
print("promoters inside cached ReMap intervals:",len(tss))
FOC=["ETS1","NFKB1","RELA","SP1"]
rows=[]
for r in tss.itertuples():
    s,e=r.tss-2000,r.tss+2000
    sub=rm[(rm.chrom==r.chrom)&(rm.end>s)&(rm.start<e)]
    tfs=set(sub.tf)
    rows.append(dict(symbol=r.symbol,chrom=r.chrom,tss=r.tss,n_remap_peaks=len(sub),n_distinct_TF=len(tfs),
                     **{f"has_{f}":int(f in tfs) for f in FOC}))
D=pd.DataFrame(rows).sort_values("n_distinct_TF",ascending=False)
D["pct_n_distinct_TF"]=D.n_distinct_TF.rank(pct=True)*100
D.to_csv(f"{R}/atac_remap_promoter_occupancy.csv",index=False)
print(D.head(12).to_string(index=False))
print("\n== focal genes ==")
print(D[D.symbol.isin(["COL1A1","COL3A1"])].to_string(index=False))
for f in FOC:
    q=D[f"has_{f}"].mean()
    print(f"{f}: present at {D[f'has_{f}'].sum()}/{len(D)} ({q*100:.0f}%) of background promoters in these intervals")
print("\nHOT-region check: Spearman(n_distinct_TF, has_focal_TF)")
from scipy.stats import spearmanr
for f in FOC:
    rho,p=spearmanr(D.n_distinct_TF,D[f"has_{f}"]); print(f"  {f}: rho={rho:.2f} p={p:.2g}")
