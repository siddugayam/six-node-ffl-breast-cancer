#!/usr/bin/env python3
"""28b_seqreg_motif_calibration.py -- higher-resolution GC-matched calibration of the
focus-TF motif counts at all ten promoters.

28_ used 1,000 background promoters, which left only 74 GC-matched controls for the very
AT-rich COL3A1 promoter, so its empirical p could not fall below 0.054. This repeats the
calibration for the NF-kB / SP / ETS focus motifs against 5,000 random RefSeq-Select
promoters and reports the matched-set size explicitly, plus a Poisson-model p-value that
does not depend on the matched-set size.
"""
import os, re, time
import numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
PROM_UP,PROM_DN=1000,500; PTHRESH=1e-4; SCALE=100
NBG=int(os.environ.get("NBG_MOTIF2","5000"))

def read_jaspar(path):
    ms={};name=None;rows=[]
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name: ms[name]=np.array(rows,dtype=float)
            p=line[1:].split(); name=f"{p[0]}|{p[1]}"; rows=[]
        elif line.strip():
            rows.append([float(x) for x in re.findall(r"[-\d.]+", line.split("[")[-1])])
    if name: ms[name]=np.array(rows,dtype=float)
    return {k:v.T for k,v in ms.items()}
def read_hocomoco(path):
    ms={};name=None;rows=[]
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name: ms[name]=np.array(rows,dtype=float)
            name=line[1:].strip(); rows=[]
        elif line.strip(): rows.append([float(x) for x in line.split()])
    if name: ms[name]=np.array(rows,dtype=float)
    return ms
MAPI={'A':0,'C':1,'G':2,'T':3}
def encode(s):
    a=np.full(len(s),4,dtype=np.int8); b=np.frombuffer(s.encode(),dtype=np.uint8)
    for c,i in MAPI.items(): a[b==ord(c)]=i
    return a
def revcomp(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def pfm_to_pwm(pfm,bg,pseudo=0.1):
    n=pfm.sum(axis=1,keepdims=True); n[n==0]=1
    p=(pfm+pseudo*bg[None,:])/(n+pseudo); return np.log2(p/bg[None,:])
def hoc_to_pwm(lw,bg):
    p=np.exp(lw)*0.25; p=p/p.sum(axis=1,keepdims=True)
    p=np.clip(p,1e-6,None); p=p/p.sum(axis=1,keepdims=True); return np.log2(p/bg[None,:])
def score_threshold(pwm,bg,pval):
    q=np.round(pwm*SCALE).astype(int); off=q.min(axis=1); q=q-off[:,None]
    maxs=int(q.max(axis=1).sum()); dist=np.zeros(maxs+1); dist[0]=1.0
    for l in range(q.shape[0]):
        nd=np.zeros(maxs+1)
        for j in range(4):
            nd[q[l,j]:q[l,j]+len(dist)-q[l,j]] += bg[j]*dist[:len(dist)-q[l,j]]
        dist=nd
    tail=np.cumsum(dist[::-1])[::-1]
    ok=np.nonzero(tail<=pval)[0]
    cut=int(ok[0]) if len(ok) else int(maxs)
    return (cut+off.sum())/SCALE
def count_sites(idx,pwm,cut):
    L=pwm.shape[0]; M=len(idx)
    if M<L: return 0
    lut=np.concatenate([pwm,np.full((L,1),-1e6)],axis=1)
    sc=np.zeros(M-L+1)
    for l in range(L): sc+=lut[l][idx[l:M-L+1+l]]
    pr=pwm[::-1,::-1]; lut2=np.concatenate([pr,np.full((L,1),-1e6)],axis=1)
    sc2=np.zeros(M-L+1)
    for l in range(L): sc2+=lut2[l][idx[l:M-L+1+l]]
    return int(((sc>=cut)|(sc2>=cut)).sum())

JAS=read_jaspar(f"{CACHE}/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
HOC=read_hocomoco(f"{CACHE}/HOCOMOCOv11_core_pwms_HUMAN_mono.txt")

reg=pd.read_csv(f"{RES}/seqreg_regions.csv"); HALF=300_000
targets={}
for _,r in reg.iterrows():
    s=open(f"{CACHE}/seq/{r.chrom}_{int(r.anchor)-HALF}_{int(r.anchor)+HALF}.txt").read().strip()
    c=HALF
    targets[r.region]= s[c-PROM_UP:c+PROM_DN] if r.strand=='+' else revcomp(s[c-PROM_DN:c+PROM_UP])

import pyfaidx
fa=pyfaidx.Fasta(f"{CACHE}/genome/hg38.fa", as_raw=True, sequence_always_upper=True)
rs=pd.read_csv(f"{CACHE}/genome/ncbiRefSeqSelect.txt.gz", sep="\t", header=None,
    names=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount",
           "exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"])
rs=rs[rs.name.str.startswith("NM_")]
MAIN=[f"chr{i}" for i in list(range(1,23))+["X"]]
rs=rs[rs.chrom.isin(MAIN)].drop_duplicates("name2")
rs["tss"]=np.where(rs.strand=="-",rs.txEnd,rs.txStart)
bgseqs=[]
for _,r in rs.sample(n=min(len(rs),NBG*2), random_state=11).iterrows():
    if len(bgseqs)>=NBG: break
    t=int(r.tss)
    if t-PROM_UP<0 or t+PROM_UP>len(fa[r.chrom]): continue
    s=str(fa[r.chrom][t-PROM_UP:t+PROM_DN]) if r.strand=='+' else revcomp(str(fa[r.chrom][t-PROM_DN:t+PROM_UP]))
    if 'N' in s or len(s)!=PROM_UP+PROM_DN: continue
    bgseqs.append(s)
print("background promoters:",len(bgseqs), flush=True)
allbg="".join(bgseqs); cnt=np.array([allbg.count(x) for x in "ACGT"],float); BG0=cnt/cnt.sum()
bg_gc=np.array([(s.count('G')+s.count('C'))/len(s) for s in bgseqs])
bg_idx=[encode(s) for s in bgseqs]

FOC=r"(?i)^(NFKB1|NFKB2|RELA|RELB|REL|SP1|SP2|SP3|ETS1|ETS2|ELK1|ELK4|GABPA)($|_)"
motifs={}
for k,v in JAS.items():
    if re.match(FOC,k.split("|")[-1]): motifs[("JASPAR",k)]=pfm_to_pwm(v,BG0)
for k,v in HOC.items():
    if re.match(FOC,k.split("|")[-1]): motifs[("HOCOMOCO",k)]=hoc_to_pwm(v,BG0)
print("focus motifs:",len(motifs), flush=True)

rows=[]; t0=time.time()
for i,((db,name),pwm) in enumerate(motifs.items()):
    if pwm.shape[0]>30: pwm=pwm[:30]
    cut=score_threshold(pwm,BG0,PTHRESH)
    counts=np.array([count_sites(x,pwm,cut) for x in bg_idx])
    for regname,s in targets.items():
        n=count_sites(encode(s),pwm,cut)
        gc=(s.count('G')+s.count('C'))/len(s)
        for tol in (0.03,0.05):
            sel=np.abs(bg_gc-gc)<=tol
            if sel.sum()>=100: break
        lam=float(counts[sel].mean())
        from scipy.stats import poisson
        rows.append(dict(db=db, motif=name, tf=name.split("|")[-1].split("_HUMAN")[0], region=regname,
                         n_sites=n, promoter_gc=round(gc,4), gc_tolerance=tol,
                         n_gc_matched=int(sel.sum()),
                         bg_gcmatched_mean=round(lam,4),
                         bg_all_mean=round(float(counts.mean()),4),
                         pct_gcmatched=round(float((counts[sel]<n).mean()*100),2),
                         emp_p_gcmatched=round(float((counts[sel]>=n).mean()),5),
                         emp_p_resolution=round(1/max(1,int(sel.sum())),5),
                         poisson_p_ge=round(float(poisson.sf(n-1,lam)) if n>0 else 1.0,5)))
    if (i+1)%5==0: print(f"  {i+1}/{len(motifs)} {time.time()-t0:.0f}s", flush=True)
out=pd.DataFrame(rows)
out.to_csv(f"{RES}/seqreg_motif_focus_calibrated_5000.csv", index=False)
pd.set_option("display.width",250)
for reg in ["COL1A1","COL3A1"]:
    print(f"\n===== {reg}")
    print(out[out.region==reg][["db","tf","n_sites","promoter_gc","n_gc_matched","bg_gcmatched_mean",
                                "pct_gcmatched","emp_p_gcmatched","emp_p_resolution","poisson_p_ge"]].to_string(index=False))
