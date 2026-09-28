#!/usr/bin/env python3
"""29_seqreg_cpg_architecture.py -- promoter sequence architecture at the COL1A1/COL3A1 and
miR-29 / miR-130a loci: GC content, CpG observed/expected, CpG-island overlap (UCSC cpgIslandExt),
ENCODE cCRE class, and the position of the Illumina 450k probes used in the project's
methylation analysis. Answers the "is the miR-130a promoter CpG-dense and silenceable?" part
of the task from sequence, independently of any deep-learning model.
Uses only cached hg38 sequence (UCSC REST) - no model, no key.
"""
import os
import numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
HALF=300_000

reg=pd.read_csv(f"{RES}/seqreg_regions.csv")
def seq(chrom,anchor):
    return open(f"{CACHE}/seq/{chrom}_{anchor-HALF}_{anchor+HALF}.txt").read().strip()
def revcomp(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

def stats(s):
    n=len(s); g=s.count('G'); c=s.count('C'); cg=s.count('CG')
    gc=(g+c)/n
    oe=(cg*n/(g*c)) if g*c>0 else np.nan
    return dict(length=n, gc=round(gc,4), n_CpG=cg, cpg_oe=round(oe,4) if oe==oe else np.nan,
                cpg_per_kb=round(1000*cg/n,3))

rows=[]
for _,r in reg.iterrows():
    s=seq(r.chrom,int(r.anchor)); c=HALF
    for label,(up,dn) in {"TSS_-1000_+500":(1000,500),"TSS_-2000_+2000":(2000,2000),
                          "TSS_-500_+500":(500,500),"TSS_-5000_+5000":(5000,5000)}.items():
        sub = s[c-up:c+dn] if r.strand=='+' else revcomp(s[c-dn:c+up])
        d=stats(sub); d.update(region=r.region, window=label, chrom=r.chrom,
                               anchor=int(r.anchor), strand=r.strand)
        rows.append(d)
arch=pd.DataFrame(rows)

# genome-wide promoter reference for the same window, so GC/CpG are percentiles not raw numbers
import pyfaidx
fa=pyfaidx.Fasta(f"{CACHE}/genome/hg38.fa", as_raw=True, sequence_always_upper=True)
rs=pd.read_csv(f"{CACHE}/genome/ncbiRefSeqSelect.txt.gz", sep="\t", header=None,
               names=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd",
                      "exonCount","exonStarts","exonEnds","score","name2","cdsStartStat",
                      "cdsEndStat","exonFrames"])
rs=rs[rs.name.str.startswith("NM_")]
MAIN=[f"chr{i}" for i in list(range(1,23))+["X"]]
rs=rs[rs.chrom.isin(MAIN)].drop_duplicates("name2")
rs["tss"]=np.where(rs.strand=="-", rs.txEnd, rs.txStart)
bg=rs.sample(n=3000, random_state=7)
brows=[]
for _,r in bg.iterrows():
    t=int(r.tss)
    if t-2000<0 or t+2000>len(fa[r.chrom]): continue
    sub=str(fa[r.chrom][t-1000:t+500]) if r.strand=='+' else revcomp(str(fa[r.chrom][t-500:t+1000]))
    if 'N' in sub: continue
    d=stats(sub); d['gene']=r.name2; brows.append(d)
B=pd.DataFrame(brows)
B.to_csv(f"{CACHE}/promoter_background_cpg.csv", index=False)
print("background promoters:", len(B), "median GC", B.gc.median(), "median CpG o/e", B.cpg_oe.median())

sub=arch[arch.window=="TSS_-1000_+500"].copy()
sub["gc_percentile"]=[round(float((B.gc<v).mean()*100),2) for v in sub.gc]
sub["cpg_oe_percentile"]=[round(float((B.cpg_oe<v).mean()*100),2) for v in sub.cpg_oe]
sub["cpg_per_kb_percentile"]=[round(float((B.cpg_per_kb<v).mean()*100),2) for v in sub.cpg_per_kb]
sub["bg_median_gc"]=round(B.gc.median(),4); sub["bg_median_cpg_oe"]=round(B.cpg_oe.median(),4)

# CpG island and cCRE overlap of the promoter window
cpg=pd.read_csv(f"{RES}/seqreg_cpg_islands.csv")
ccre=pd.read_csv(f"{RES}/seqreg_encode_ccre.csv")
def overlap(df,r,up=1000,dn=500):
    a=int(r.anchor); s,e=(a-up,a+dn) if r.strand=='+' else (a-dn,a+up)
    d=df[(df.region==r.region)&(df.start<e)&(df.end>s)]
    return d
isl=[];cc=[]
for _,r in sub.iterrows():
    d=overlap(cpg,r); isl.append(";".join(f"{x.chrom}:{x.start}-{x.end}(obsExp={x.obsExp},perGc={x.perGc},len={x.length})" for x in d.itertuples()))
    d2=overlap(ccre,r); cc.append(";".join(f"{x.accession}:{x.ccre_class}" for x in d2.itertuples()))
sub["cpg_island_overlap"]=isl; sub["ccre_overlap"]=cc
arch.to_csv(f"{RES}/seqreg_promoter_architecture_all_windows.csv", index=False)
sub.to_csv(f"{RES}/seqreg_promoter_architecture.csv", index=False)
pd.set_option("display.width",250)
print(sub[["region","gc","gc_percentile","cpg_oe","cpg_oe_percentile","cpg_per_kb",
           "cpg_per_kb_percentile","cpg_island_overlap","ccre_overlap"]].to_string(index=False))
