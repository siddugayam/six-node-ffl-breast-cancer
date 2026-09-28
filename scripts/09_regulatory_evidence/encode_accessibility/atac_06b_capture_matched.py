#!/usr/bin/env python3
"""Reproducibility-matched version of the compartment-capture comparison:
elements called in ALL 5 samples of one compartment and in NONE of the other."""
import json, os, subprocess, tempfile, pandas as pd, numpy as np
ROOT="/path/to/revision"; C=f"{ROOT}/cache/v6/atac"; E=f"{C}/encode"; R=f"{ROOT}/results/v6"
man=pd.DataFrame(json.load(open(f"{E}/manifest.json")))
EPI=["MCF-7","MCF 10A","breast epithelium"]; FIB=["fibroblast of mammary gland","fibroblast of dermis","fibroblast of lung","IMR-90"]
man["compartment"]=np.where(man.biosample.isin(EPI),"epithelial","fibroblast")
tmp=tempfile.mkdtemp()
def prep(f):
    d=pd.read_csv(f"{E}/{f}.bed.gz",sep="\t",header=None,usecols=[0,1,2],names=["chrom","start","end"])
    d=d[d.chrom.str.match(r"^chr[0-9XY]+$")].sort_values(["chrom","start"])
    p=f"{tmp}/{f}.bed"; d.to_csv(p,sep="\t",header=False,index=False); return p
files={r.file:prep(r.file) for _,r in man.iterrows()}
def run(cmd,out):
    with open(out,"w") as fh: subprocess.run(cmd,stdout=fh,check=True)
    return out
def core(comp):
    fs=[files[f] for f in man[man.compartment==comp].file]
    cur=fs[0]
    for i,f in enumerate(fs[1:]):
        cur=run(["bedtools","intersect","-a",cur,"-b",f,"-u"],f"{tmp}/{comp}_core{i}.bed")
    return run(["bedtools","merge","-i",cur],f"{tmp}/{comp}_core.bed")
def anyof(comp):
    fs=[files[f] for f in man[man.compartment==comp].file]
    cat=f"{tmp}/{comp}_cat.bed"
    with open(cat,"w") as fh: subprocess.run(["cat"]+fs,stdout=fh)
    srt=run(["bedtools","sort","-i",cat],f"{tmp}/{comp}_srt.bed")
    return run(["bedtools","merge","-i",srt],f"{tmp}/{comp}_any.bed")
fc,ec=core("fibroblast"),core("epithelial")
fa,ea=anyof("fibroblast"),anyof("epithelial")
fspec=run(["bedtools","intersect","-a",fc,"-b",ea,"-v"],f"{tmp}/fspec.bed")
espec=run(["bedtools","intersect","-a",ec,"-b",fa,"-v"],f"{tmp}/espec.bed")
brca=pd.read_csv(f"{C}/brca_peaks.bed",sep="\t",header=None,names=["chrom","start","end","name"]).sort_values(["chrom","start"])
bb=f"{tmp}/brca.bed"; brca.to_csv(bb,sep="\t",header=False,index=False)
rows=[]
for tag,f in [("fibroblast_specific_all5",fspec),("epithelial_specific_all5",espec)]:
    n=sum(1 for _ in open(f))
    h=subprocess.run(["bedtools","intersect","-a",f,"-b",bb,"-u"],capture_output=True,text=True).stdout
    nh=len([l for l in h.strip().split("\n") if l])
    # median width, to check the two sets are comparable
    w=pd.read_csv(f,sep="\t",header=None,names=["c","s","e"]); rows.append(dict(set=tag,n=n,captured=nh,frac=nh/n,median_width=int((w.e-w.s).median())))
D=pd.DataFrame(rows); D.to_csv(f"{R}/atac_encode_compartment_capture_matched.csv",index=False)
print(D.to_string(index=False))
from scipy.stats import fisher_exact
a,b=D.captured.values; na,nb=D.n.values
odds,p=fisher_exact([[a,na-a],[b,nb-b]])
print(f"Fisher OR (fib vs epi capture) = {odds:.3f}, p = {p:.3g}")
open(f"{R}/atac_encode_compartment_capture_matched.csv","a").write(f"# Fisher OR={odds:.4f} p={p:.3g}\n")
