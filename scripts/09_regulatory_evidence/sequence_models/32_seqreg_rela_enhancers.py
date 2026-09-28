#!/usr/bin/env python3
"""32_seqreg_rela_enhancers.py -- the only fibroblast-context TF-binding evidence at the
collagen loci in ReMap 2022 is RELA in TNF-stimulated IMR-90 lung fibroblasts (GSE43070) and
nutlin-treated human dermal fibroblasts (GSE77225). This script (a) tabulates those peaks,
(b) scans each peak sequence for NF-kB-family motifs with the same FIMO-equivalent procedure
used in 28_, and (c) reports which ENCODE cCRE each peak falls in and the Enformer-predicted
fibroblast DNase signal there.
"""
import re, numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
COLS=["chrom","start","end","name","score","strand","tstart","tend","rgb","tf","biotype"]
FIBRO={"BJ","BJ1-hTERT","GM23248","HDF","HFF","IMR-90","MRC-5","WI-38","WI-38VA13",
       "dermal-fibroblast","fibroblast","hMSC","hMSC-TERT","hMSC-TERT4","mesenchymal",
       "myofibroblast","primary-dermal-fibroblasts","primary-lung-fibroblast",
       "proliferating-human-fibroblast"}
SCALE=100; PTHRESH=1e-4
def read_jaspar(path):
    ms={};name=None;rows=[]
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name: ms[name]=np.array(rows,dtype=float)
            p=line[1:].split(); name=f"{p[0]}|{p[1]}"; rows=[]
        elif line.strip(): rows.append([float(x) for x in re.findall(r"[-\d.]+", line.split("[")[-1])])
    if name: ms[name]=np.array(rows,dtype=float)
    return {k:v.T for k,v in ms.items()}
def read_hoc(path):
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
def pwm_from_pfm(pfm,bg,pseudo=0.1):
    n=pfm.sum(axis=1,keepdims=True); n[n==0]=1
    return np.log2((pfm+pseudo*bg[None,:])/(n+pseudo)/bg[None,:])
def pwm_from_hoc(lw,bg):
    p=np.exp(lw)*0.25; p=p/p.sum(axis=1,keepdims=True); p=np.clip(p,1e-6,None); p=p/p.sum(axis=1,keepdims=True)
    return np.log2(p/bg[None,:])
def thr(pwm,bg,pv):
    q=np.round(pwm*SCALE).astype(int); off=q.min(axis=1); q=q-off[:,None]
    maxs=int(q.max(axis=1).sum()); d=np.zeros(maxs+1); d[0]=1.0
    for l in range(q.shape[0]):
        nd=np.zeros(maxs+1)
        for j in range(4): nd[q[l,j]:q[l,j]+len(d)-q[l,j]] += bg[j]*d[:len(d)-q[l,j]]
        d=nd
    tail=np.cumsum(d[::-1])[::-1]; ok=np.nonzero(tail<=pv)[0]
    cut=int(ok[0]) if len(ok) else maxs
    return (cut+off.sum())/SCALE, tail, int(off.sum())
def best_hits(idx,pwm,cut,tail,offsum,seq):
    L=pwm.shape[0]; M=len(idx)
    if M<L: return []
    lut=np.concatenate([pwm,np.full((L,1),-1e6)],axis=1); sc=np.zeros(M-L+1)
    for l in range(L): sc+=lut[l][idx[l:M-L+1+l]]
    pr=pwm[::-1,::-1]; lut2=np.concatenate([pr,np.full((L,1),-1e6)],axis=1); sc2=np.zeros(M-L+1)
    for l in range(L): sc2+=lut2[l][idx[l:M-L+1+l]]
    out=[]
    for pos in np.where((sc>=cut)|(sc2>=cut))[0]:
        s=float(max(sc[pos],sc2[pos]))
        pv=float(tail[int(np.clip(round(s*SCALE)-offsum,0,len(tail)-1))])
        out.append(dict(offset=int(pos), strand='+' if sc[pos]>=sc2[pos] else '-',
                        score=round(s,3), pvalue=pv, match=seq[pos:pos+L]))
    return out

JAS=read_jaspar(f"{CACHE}/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
HOC=read_hoc(f"{CACHE}/HOCOMOCOv11_core_pwms_HUMAN_mono.txt")
NK={}
for k,v in JAS.items():
    if re.match(r"^(NFKB1|NFKB2|RELA|RELB|REL)$", k.split("|")[-1]): NK[("JASPAR",k)]=v
HK={}
for k,v in HOC.items():
    if re.match(r"^(NFKB1|NFKB2|RELA|RELB|REL)_HUMAN", k): HK[("HOCOMOCO",k)]=v

import pyfaidx
fa=pyfaidx.Fasta(f"{CACHE}/genome/hg38.fa", as_raw=True, sequence_always_upper=True)
bgprom=pd.read_csv(f"{CACHE}/promoter_background_cpg.csv") if False else None
# 0-order background: the same 1,500 bp promoter set composition used in 28_ is not
# appropriate for enhancers; use the composition of the peak set itself.
p=pd.read_csv(f"{CACHE}/remap2022_hg38_loci.bed", sep="\t", header=None, names=COLS)
reg=pd.read_csv(f"{RES}/seqreg_regions.csv").set_index("region")
ccre=pd.read_csv(f"{RES}/seqreg_encode_ccre.csv").drop_duplicates("accession")
targets=pd.read_csv(f"{CACHE}/enformer_targets_human.txt", sep="\t")
targets['assay']=targets.description.str.split(':').str[0]
FIB_DNASE=targets.index[(targets.assay=='DNASE')&targets.description.str.contains('ibroblast')].to_numpy()
SEQ_LEN=196_608; N_BINS=896; BIN=128; CROP=(SEQ_LEN-N_BINS*BIN)//2

rows=[]; seqs=[]
for g in ["COL1A1","COL3A1"]:
    r=reg.loc[g]; a=int(r.anchor)
    sub=p[(p.chrom==r.chrom)&(p.start<a+100000)&(p.end>a-100000)&(p.biotype.isin(FIBRO))&
          (p.tf.isin(["NFKB1","NFKB2","RELA","RELB","REL"]))]
    pred=np.load(f"{CACHE}/enformer/{g}_full.npy"); fib=pred[:,FIB_DNASE].mean(axis=1)
    win_start=a-SEQ_LEN//2+CROP
    for _,q in sub.iterrows():
        s=str(fa[q.chrom][int(q.start):int(q.end)])
        mid=(int(q.start)+int(q.end))//2
        b=int((mid-win_start)//BIN)
        rows.append(dict(gene=g, tf=q.tf, biotype=q.biotype, dataset=q["name"],
                         chrom=q.chrom, start=int(q.start), end=int(q.end), width=int(q.end-q.start),
                         dist_to_TSS=(mid-a)*(1 if r.strand=='+' else -1),
                         enformer_fibroblast_dnase=round(float(fib[b]),4) if 0<=b<N_BINS else np.nan,
                         ccre=";".join(ccre[(ccre.chrom==q.chrom)&(ccre.start<q.end)&(ccre.end>q.start)]
                                       .apply(lambda x: f"{x.accession}({x.ccre_class})", axis=1))))
        seqs.append(s)
D=pd.DataFrame(rows)
allseq="".join(seqs); cnt=np.array([allseq.count(c) for c in "ACGT"],float); BG0=cnt/cnt.sum()
print("peak-set base composition ACGT:", np.round(BG0,4))
mot={}
for k,v in NK.items(): mot[k]=pwm_from_pfm(v,BG0)
for k,v in HK.items(): mot[k]=pwm_from_hoc(v,BG0)
hits=[]
for (db,name),pwm in mot.items():
    cut,tail,offsum=thr(pwm,BG0,PTHRESH)
    for i,s in enumerate(seqs):
        for h in best_hits(encode(s),pwm,cut,tail,offsum,s):
            hits.append(dict(peak_index=i, db=db, motif=name, **h))
H=pd.DataFrame(hits)
if len(H):
    agg=H.groupby("peak_index").agg(n_kB_motifs=("motif","size"), n_matrices=("motif","nunique"),
                                    best_p=("pvalue","min"),
                                    best_match=("match", lambda x: list(x)[0])).reset_index()
    D=D.merge(agg, left_index=True, right_on="peak_index", how="left")
else:
    D["n_kB_motifs"]=0; D["n_matrices"]=0; D["best_p"]=np.nan; D["best_match"]=""
D["n_kB_motifs"]=D.n_kB_motifs.fillna(0).astype(int)
D.to_csv(f"{RES}/seqreg_rela_fibroblast_peaks.csv", index=False)
H.to_csv(f"{RES}/seqreg_rela_fibroblast_peak_motifs.csv", index=False)
pd.set_option("display.width",260)
print(D[["gene","tf","biotype","dataset","chrom","start","end","dist_to_TSS",
         "enformer_fibroblast_dnase","n_kB_motifs","n_matrices","best_p","best_match","ccre"]]
      .sort_values(["gene","dist_to_TSS"]).to_string(index=False))
print("\npeaks with >=1 NF-kB consensus motif at p<1e-4:", int((D.n_kB_motifs>0).sum()), "of", len(D))
