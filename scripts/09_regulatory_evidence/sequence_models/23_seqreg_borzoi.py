#!/usr/bin/env python3
"""23_seqreg_borzoi.py -- Borzoi (Linder et al., Nat Genet 2025) sequence-to-function
predictions at the COL1A1/COL3A1/miR-29 loci.

Borzoi is used in addition to Enformer because its human head contains RELA and RELB ChIP
tracks (Enformer's does not contain any NF-kB-family track at all) and 1,543 RNA-seq tracks
including GTEx fibroblast and breast, so the NFKB1/RELA -> collagen edges of the network can be
put to a sequence-based test.
Input 524,288 bp, output 6,144 bins x 32 bp (196,608 bp) x 7,611 human tracks.
Weights: johahi/borzoi-replicate-{0..3} on HuggingFace (Calico weights ported to PyTorch).
"""
import os, time
import numpy as np, pandas as pd, torch
from borzoi_pytorch import Borzoi

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
SEQ_LEN=524_288; N_BINS=6144; BIN=32
CROP=(SEQ_LEN-N_BINS*BIN)//2          # 163,840
CENTRE_BIN=(SEQ_LEN//2-CROP)//BIN     # 3072
NBG=int(os.environ.get("NBG_BORZOI","200"))
REPS=[int(x) for x in os.environ.get("BORZOI_REPS","0").split(",")]

targets=pd.read_csv(f"{CACHE}/borzoi_targets_human.txt", sep="\t")
targets['assay']=targets.description.str.split(':').str[0]
targets['tf']=np.where(targets.assay=='CHIP', targets.description.str.split(':').str[1], '')

MAP={'A':0,'C':1,'G':2,'T':3}
def onehot(s):
    a=np.zeros((4,len(s)),dtype=np.float32)
    idx=np.frombuffer(s.encode(),dtype=np.uint8)
    for b,i in MAP.items(): a[i, idx==ord(b)]=1.0
    return a

def load_rep(r):
    """transformers' from_pretrained leaves the 8 non-persistent relative-position
    buffers on the meta device; rebuild them from a freshly constructed model."""
    m=Borzoi.from_pretrained(f'johahi/borzoi-replicate-{r}', low_cpu_mem_usage=False)
    fresh=Borzoi(m.config); fb=dict(fresh.named_buffers()); nfix=0
    for name,buf in list(m.named_buffers()):
        if buf.is_meta:
            mod=m; parts=name.split('.')
            for pp in parts[:-1]: mod=getattr(mod,pp)
            setattr(mod, parts[-1], fb[name].clone()); nfix+=1
    assert not any(b.is_meta for _,b in m.named_buffers())
    del fresh
    print(f"replicate {r}: rebuilt {nfix} meta buffers", flush=True)
    return m.eval().cuda()

models=[load_rep(r) for r in REPS]
print("loaded replicates:", REPS, flush=True)

@torch.no_grad()
def predict(seq):
    x=torch.from_numpy(onehot(seq)).unsqueeze(0).cuda()
    outs=[]
    for m in models:
        with torch.autocast('cuda', dtype=torch.bfloat16):
            o=m(x)
        o=o.float().cpu().numpy()[0]
        if o.shape[0]==targets.shape[0]: o=o.T      # -> (bins, tracks)
        outs.append(o)
    return np.mean(outs, axis=0)

def promoter_value(pred, half_bins=8):   # +/- 256 bp
    return pred[CENTRE_BIN-half_bins:CENTRE_BIN+half_bins+1].mean(axis=0)

reg=pd.read_csv(f"{RES}/seqreg_regions.csv")
HALF=300_000
def locus_seq(chrom,anchor):
    s=open(f"{CACHE}/seq/{chrom}_{anchor-HALF}_{anchor+HALF}.txt").read().strip()
    c=HALF; return s[c-SEQ_LEN//2:c+SEQ_LEN//2]

os.makedirs(f"{CACHE}/borzoi", exist_ok=True)
loc={}
for _,r in reg.iterrows():
    s=locus_seq(r.chrom,int(r.anchor)); assert len(s)==SEQ_LEN
    t0=time.time(); p=predict(s)
    np.save(f"{CACHE}/borzoi/{r.region}_full.npy", p.astype(np.float32))
    loc[r.region]=promoter_value(p)
    print(f"{r.region} {time.time()-t0:.1f}s shape={p.shape} max={p.max():.2f}", flush=True)

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
bg=rs.sample(n=NBG*2, random_state=20260909)
bgv=[]; names=[]
t0=time.time()
for _,r in bg.iterrows():
    if len(bgv)>=NBG: break
    t=int(r.tss)
    if t-SEQ_LEN//2<0 or t+SEQ_LEN//2>len(fa[r.chrom]): continue
    s=str(fa[r.chrom][t-SEQ_LEN//2:t+SEQ_LEN//2])
    if len(s)!=SEQ_LEN or s.count('N')>SEQ_LEN*0.1: continue
    bgv.append(promoter_value(predict(s))); names.append(r.name2)
    if len(bgv)%25==0: print(f"  bg {len(bgv)}/{NBG} {time.time()-t0:.0f}s", flush=True)
BG=np.vstack(bgv)
np.save(f"{CACHE}/borzoi/background_promoter_values.npy", BG)
pd.Series(names).to_csv(f"{CACHE}/borzoi/background_promoter_genes.csv", index=False, header=["gene"])
print("BG", BG.shape, flush=True)

rows=[]
for regname,v in loc.items():
    pct=(BG<v[None,:]).mean(axis=0)*100
    d=targets[['identifier','description','assay','tf']].copy()
    d['region']=regname; d['value']=v; d['bg_percentile']=pct
    d['bg_mean']=BG.mean(axis=0); d['bg_sd']=BG.std(axis=0)
    rows.append(d)
out=pd.concat(rows, ignore_index=True)
out.to_csv(f"{RES}/seqreg_borzoi_promoter_tracks.csv.gz", index=False, compression="gzip")
foc=out[(out.tf.isin(['RELA','RELB','NFKB2','NFKBIZ','SP1','SP2','SP3','ETS1','ETS2','ELK1',
                      'GABPA','GABPB1','CTCF','EP300','POLR2A','JUND','FOSL2','SMAD3','TWIST1',
                      'RUNX1','TEAD1','SRF','MYC','STAT3','TP53','ESR1','EGR1','KLF4','YY1','MED1'])) |
        (out.assay.isin(['DNASE','ATAC','CAGE','RNA']))]
foc.to_csv(f"{RES}/seqreg_borzoi_focus_tracks.csv.gz", index=False, compression="gzip")
print("wrote", len(out), "rows; focus", len(foc))
