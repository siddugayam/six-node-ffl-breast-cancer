#!/usr/bin/env python3
"""24_seqreg_chipatlas.py -- ChIP-Atlas (Zou et al., NAR 2022) evidence that a TF is bound
within +/-1 kb and +/-5 kb of the COL1A1 and COL3A1 TSS, resolved by cell type, with
particular attention to fibroblast / mesenchymal cell types.

Two sources are used:
 (1) ChIP-Atlas 'target genes' tables https://chip-atlas.dbcls.jp/data/hg38/target/{TF}.{d}.tsv
     (rows = genes, columns = individual SRX experiments annotated by cell type, values =
     MACS2 -log10(Q) of the strongest peak in the window). Downloaded per TF and streamed.
 (2) ChIP-Atlas experimentList.tab, used to quantify how much fibroblast ChIP-seq exists at all,
     so that an absence of fibroblast evidence can be interpreted.
"""
import io, os, re, sys, time, urllib.request, gzip
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
UA={"User-Agent":"Mozilla/5.0 (research; contact your.email@example.org)"}
GENES=["COL1A1","COL3A1"]

# ---- (2) experiment inventory ----------------------------------------------
# the full experimentList.tab has a ragged number of trailing attribute columns, so the
# first seven fields for hg38 were extracted with awk into chipatlas_hg38.tsv beforehand:
#   awk -F'\t' 'BEGIN{OFS="\t"} $2=="hg38"{print $1,$2,$3,$4,$5,$6,$7}' chipatlas_experimentList.tab
hg=pd.read_csv(f"{CACHE}/chipatlas_hg38.tsv", sep="\t", header=None,
               names=["srx","genome","antigen_class","antigen","celltype_class",
                      "celltype","celltype_desc"], dtype=str, quoting=3)
tf=hg[hg.antigen_class=="TFs and others"]
print("hg38 experiments:", len(hg), " TF experiments:", len(tf), " unique antigens:", tf.antigen.nunique())
FIB=r"(?i)(fibroblast|IMR-?90|MRC-?5|WI-?38|^BJ$|HFF|NHDF|NHLF|myofibroblast|mesenchymal|hMSC|MSC)"
tf_fib=tf[tf.celltype.str.contains(FIB, na=False, regex=True)]
print("TF ChIP experiments in fibroblast/mesenchymal cell types:", len(tf_fib),
      " unique antigens:", tf_fib.antigen.nunique())
inv=pd.DataFrame({
 "metric":["hg38_experiments_total","hg38_TF_experiments","hg38_TF_unique_antigens",
           "hg38_TF_experiments_fibroblast_mesenchymal","hg38_TF_unique_antigens_fibroblast_mesenchymal"],
 "value":[len(hg),len(tf),tf.antigen.nunique(),len(tf_fib),tf_fib.antigen.nunique()]})
inv.to_csv(f"{RES}/seqreg_chipatlas_inventory.csv", index=False)
tf_fib.groupby(["celltype","antigen"]).size().reset_index(name="n_experiments")\
     .to_csv(f"{RES}/seqreg_chipatlas_fibroblast_experiments.csv", index=False)

ANTIGENS=sorted(tf.antigen.dropna().unique())
print("antigens to query:", len(ANTIGENS))

# ---- (1) per-TF target-gene tables -----------------------------------------
def fetch_tf(args):
    ag, dist = args
    url=f"https://chip-atlas.dbcls.jp/data/hg38/target/{urllib.request.quote(ag)}.{dist}.tsv"
    try:
        req=urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            raw=r.read().decode("utf-8", "replace")
    except Exception as e:
        return ag, dist, None, str(e)
    lines=raw.split("\n")
    if not lines: return ag,dist,None,"empty"
    header=lines[0].split("\t")
    out=[]
    for ln in lines[1:]:
        if not ln: continue
        f=ln.split("\t")
        if f[0] in GENES:
            for h,v in zip(header[1:], f[1:]):
                if v in ("","0","0.0"): continue
                try: val=float(v)
                except ValueError: continue
                if val<=0: continue
                srx, _, cell = h.partition("|")
                out.append(dict(tf=ag, distance_kb=dist, gene=f[0], srx=srx, celltype=cell, score=val))
    return ag, dist, out, None

jobs=[(a,d) for a in ANTIGENS for d in (1,5)]
res=[]; errs=[]
t0=time.time()
with ThreadPoolExecutor(max_workers=10) as ex:
    for i,(ag,dist,out,err) in enumerate(ex.map(fetch_tf, jobs)):
        if err: errs.append((ag,dist,err))
        elif out: res.extend(out)
        if (i+1)%200==0:
            print(f"  {i+1}/{len(jobs)} rows={len(res)} errs={len(errs)} {time.time()-t0:.0f}s", flush=True)
df=pd.DataFrame(res)
df.to_csv(f"{RES}/seqreg_chipatlas_binding_at_collagens.csv", index=False)
pd.DataFrame(errs, columns=["tf","distance_kb","error"]).to_csv(
    f"{RES}/seqreg_chipatlas_fetch_errors.csv", index=False)
print("binding rows:", len(df), "errors:", len(errs))
if len(df):
    df["is_fibro"]=df.celltype.str.contains(FIB, na=False, regex=True)
    summ=(df.groupby(["gene","distance_kb","tf"])
            .agg(n_experiments=("srx","nunique"), n_celltypes=("celltype","nunique"),
                 max_score=("score","max"), n_fibro_exp=("is_fibro","sum"),
                 fibro_celltypes=("celltype", lambda s: ";".join(sorted(set(
                     x for x,f in zip(s, df.loc[s.index,"is_fibro"]) if f)))[:300]))
            .reset_index().sort_values(["gene","distance_kb","max_score"], ascending=[1,1,0]))
    summ.to_csv(f"{RES}/seqreg_chipatlas_tf_summary.csv", index=False)
    print(summ.head(20).to_string(index=False))
