#!/usr/bin/env python3
"""
Ecological decomposition: per TCGA cohort compute rho(TF, CAF), rho(CAF, COL1A1),
their product, and the observed rho(TF, COL1A1). If the TF->collagen correlation is
stroma-carried, observed ~ product. Same for miR-29 (expected NOT to be).
Also caches the extracted gene submatrix for reuse.
"""
import gzip, os, json, math
import numpy as np
from scipy import stats
BASE="/path/to/revision"; ML="/path/to/home/Desktop/DD/R_GPR/ML"
OUT=os.path.join(BASE,"results/v2")
GEXP=os.path.join(ML,"EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
MEXP=os.path.join(ML,"pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz")
PHENO=os.path.join(ML,"TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
cafA=[l.strip() for l in open(os.path.join(OUT,"cafA_gene_list_v3.txt")) if l.strip()]
FOCUS=["COL1A1","COL3A1","ETS1","NFKB1","SP1","RELA","MYC","TP53"]
WANT=sorted(set(cafA)|set(FOCUS))
ph={}
with gzip.open(PHENO,"rt") as fh:
    fh.readline()
    for line in fh:
        p=line.rstrip("\n").split("\t"); ph[p[0]]=(p[2],p[3])
with gzip.open(MEXP,"rt") as fh:
    msamp=fh.readline().rstrip("\n").split("\t")[1:]
    mn,mr=[],[]
    for line in fh:
        p=line.rstrip("\n").split("\t"); mn.append(p[0])
        mr.append(np.array([np.nan if v in("","NA") else float(v) for v in p[1:]],dtype=np.float32))
Mmat=np.vstack(mr); midx={n:i for i,n in enumerate(mn)}
want=set(WANT); gn,gr=[],[]
with gzip.open(GEXP,"rt") as fh:
    gsamp=fh.readline().rstrip("\n").split("\t")[1:]
    for line in fh:
        p=line.rstrip("\n").split("\t")
        if p[0] in want:
            gn.append(p[0]); gr.append(np.array([np.nan if v in("","NA") else float(v) for v in p[1:]],dtype=np.float32))
Gmat=np.vstack(gr)
np.savez_compressed(os.path.join(OUT,"pancan_gene_cache_v3.npz"),G=Gmat,genes=np.array(gn),samples=np.array(gsamp))
print("cached gene submatrix", Gmat.shape, flush=True)
gidx={}
for i,n in enumerate(gn): gidx.setdefault(n,i)
gpos={}
for i,s in enumerate(gsamp): gpos.setdefault(s,i)
mpos={}
for i,s in enumerate(msamp): mpos.setdefault(s,i)
common=[s for s in gsamp if s in mpos and s in ph and ph[s][0]=="Primary Tumor"]
coh={}
for s in common: coh.setdefault(ph[s][1],[]).append(s)
coh={k:v for k,v in sorted(coh.items(),key=lambda kv:-len(kv[1])) if len(v)>=30}
def sp(x,y):
    k=np.isfinite(x)&np.isfinite(y)
    if k.sum()<12: return np.nan
    return float(stats.spearmanr(x[k],y[k])[0])
MIR={"hsa-miR-29a":["hsa-miR-29a-3p","hsa-miR-29a-5p"],
     "hsa-miR-29b":["hsa-miR-29b-1-5p","hsa-miR-29b-2-5p","hsa-miR-29b-3p"],
     "hsa-miR-29c":["hsa-miR-29c-3p","hsa-miR-29c-5p"],
     "hsa-miR-101":["hsa-miR-101-3p","hsa-miR-101-5p"]}
rows=[]
for c,ss in coh.items():
    gi=np.array([gpos[s] for s in ss]); mi=np.array([mpos[s] for s in ss])
    Gc=Gmat[:,gi]; Mc=Mmat[:,mi]
    zs=[]
    for g in cafA:
        i=gidx.get(g)
        if i is None: continue
        v=Gc[i].astype(np.float64); sd=np.nanstd(v)
        if sd>0: zs.append((v-np.nanmean(v))/sd)
    CAF=np.mean(np.vstack(zs),axis=0)
    col1=Gc[gidx["COL1A1"]].astype(np.float64); col3=Gc[gidx["COL3A1"]].astype(np.float64)
    r_caf_c1=sp(CAF,col1); r_caf_c3=sp(CAF,col3)
    for reg in ["ETS1","NFKB1","SP1","RELA","MYC","TP53"]+list(MIR):
        if reg in MIR:
            arms=[a for a in MIR[reg] if a in midx]
            if not arms: continue
            x=np.nanmean(np.vstack([Mc[midx[a]] for a in arms]).astype(np.float64),axis=0)
        else:
            x=Gc[gidx[reg]].astype(np.float64)
        r_reg_caf=sp(x,CAF)
        for tgt,y,rcc in (("COL1A1",col1,r_caf_c1),("COL3A1",col3,r_caf_c3)):
            obs=sp(x,y)
            rows.append(dict(cohort=c,n=len(ss),regulator=reg,target=tgt,
                rho_reg_CAF=r_reg_caf, rho_CAF_target=rcc,
                product=r_reg_caf*rcc, rho_observed=obs,
                residual=obs-r_reg_caf*rcc))
    print("done",c,flush=True)
import csv
with open(os.path.join(OUT,"pancancer_ecological_product_v3.csv"),"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote pancancer_ecological_product_v3.csv",len(rows),flush=True)
# summary
import collections
print(f"\n{'regulator':12s}{'target':8s}{'k':>4s}{'mean_obs':>10s}{'mean_prod':>11s}{'mean_resid':>12s}"
      f"{'r(obs,prod)':>13s}{'p':>10s}{'mean_rho_reg_CAF':>18s}")
by=collections.defaultdict(list)
for r in rows: by[(r["regulator"],r["target"])].append(r)
for (reg,tgt),d in by.items():
    obs=np.array([x["rho_observed"] for x in d]); pr=np.array([x["product"] for x in d])
    rc=np.array([x["rho_reg_CAF"] for x in d]); k=np.isfinite(obs)&np.isfinite(pr)
    rr,pp=stats.spearmanr(obs[k],pr[k])
    print(f"{reg:12s}{tgt:8s}{k.sum():>4d}{obs[k].mean():>10.3f}{pr[k].mean():>11.3f}"
          f"{(obs[k]-pr[k]).mean():>12.3f}{rr:>13.3f}{pp:>10.2e}{rc[k].mean():>18.3f}")
