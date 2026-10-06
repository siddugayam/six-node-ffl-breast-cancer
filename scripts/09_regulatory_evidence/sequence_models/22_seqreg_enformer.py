#!/usr/bin/env python3
"""22_seqreg_enformer.py -- Enformer (Avsec et al. 2021) sequence-to-function predictions
at the COL1A1 / COL3A1 / miR-29 loci (hg38), with a genome-wide promoter
background so that every statement is a calibrated percentile rather than a raw number.

Model: EleutherAI/enformer-official-rough (the official Enformer weights ported to PyTorch).
Input 196,608 bp; output 896 bins x 128 bp (114,688 bp) x 5,313 human tracks.
Track metadata: Basenji/Enformer targets_human.txt (cached).
"""
import os, sys, json, time
import numpy as np, pandas as pd, torch
from enformer_pytorch import from_pretrained

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
SEQ_LEN=196_608; N_BINS=896; BIN=128
CROP=(SEQ_LEN-N_BINS*BIN)//2      # 40,960
CENTRE_BIN=(SEQ_LEN//2-CROP)//BIN # 448
NBG=int(os.environ.get("NBG","500"))

targets=pd.read_csv(f"{CACHE}/enformer_targets_human.txt", sep="\t")
targets['assay']=targets.description.str.split(':').str[0]
targets['tf']=np.where(targets.assay=='CHIP', targets.description.str.split(':').str[1], '')
targets['cell']=targets.description.str.split(':').str[-1]

MAP={'A':0,'C':1,'G':2,'T':3}
def onehot(s):
    a=np.zeros((len(s),4),dtype=np.float32)
    idx=np.frombuffer(s.encode(),dtype=np.uint8)
    for b,i in MAP.items():
        a[idx==ord(b),i]=1.0
    return a

model=from_pretrained('EleutherAI/enformer-official-rough').eval().cuda()

@torch.no_grad()
def predict(seq):
    x=torch.from_numpy(onehot(seq)).unsqueeze(0).cuda()
    o=model(x)['human'][0].float().cpu().numpy()   # (896, 5313)
    return o

def promoter_value(pred, half_bins=1):
    lo=CENTRE_BIN-half_bins; hi=CENTRE_BIN+half_bins+1
    return pred[lo:hi].mean(axis=0)

# ---------------- loci -------------------------------------------------------
reg=pd.read_csv(f"{RES}/seqreg_regions.csv")
HALF=300_000
def locus_seq(chrom, anchor):
    f=f"{CACHE}/seq/{chrom}_{anchor-HALF}_{anchor+HALF}.txt"
    s=open(f).read().strip()
    c=HALF
    return s[c-SEQ_LEN//2 : c+SEQ_LEN//2]

os.makedirs(f"{CACHE}/enformer", exist_ok=True)
locus_prom={}
for _,r in reg.iterrows():
    s=locus_seq(r.chrom,int(r.anchor))
    assert len(s)==SEQ_LEN, len(s)
    t0=time.time(); p=predict(s)
    np.save(f"{CACHE}/enformer/{r.region}_full.npy", p.astype(np.float32))
    locus_prom[r.region]=promoter_value(p)
    print(f"{r.region}: {time.time()-t0:.1f}s max={p.max():.2f}", flush=True)

# ---------------- background promoters --------------------------------------
import pyfaidx
fa=pyfaidx.Fasta(f"{CACHE}/genome/hg38.fa", as_raw=True, sequence_always_upper=True)
rs=pd.read_csv(f"{CACHE}/genome/ncbiRefSeqSelect.txt.gz", sep="\t", header=None,
               names=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd",
                      "exonCount","exonStarts","exonEnds","score","name2","cdsStartStat",
                      "cdsEndStat","exonFrames"])
rs=rs[rs.name.str.startswith("NM_")]
MAINCHR=[f"chr{i}" for i in list(range(1,23))+["X"]]
rs=rs[rs.chrom.isin(MAINCHR)].drop_duplicates("name2")
rs["tss"]=np.where(rs.strand=="-", rs.txEnd, rs.txStart)
rng=np.random.default_rng(20260909)
ok=rs[(rs.tss>SEQ_LEN) & (rs.tss < rs.chrom.map({c:len(fa[c]) for c in MAINCHR})-SEQ_LEN)]
bg=ok.sample(n=min(NBG,len(ok)), random_state=20260909)
print("background promoters:", len(bg), flush=True)

bgv=[]; bgnames=[]
t0=time.time()
for i,(_,r) in enumerate(bg.iterrows()):
    s=str(fa[r.chrom][int(r.tss)-SEQ_LEN//2:int(r.tss)+SEQ_LEN//2])
    if len(s)!=SEQ_LEN or s.count('N')>SEQ_LEN*0.1: continue
    bgv.append(promoter_value(predict(s))); bgnames.append(r.name2)
    if (i+1)%50==0: print(f"  bg {i+1}/{len(bg)} {time.time()-t0:.0f}s", flush=True)
BG=np.vstack(bgv)
np.save(f"{CACHE}/enformer/background_promoter_values.npy", BG)
pd.Series(bgnames).to_csv(f"{CACHE}/enformer/background_promoter_genes.csv", index=False, header=["gene"])
print("BG matrix", BG.shape, flush=True)

# ---------------- assemble table --------------------------------------------
rows=[]
for regname,v in locus_prom.items():
    pct=(BG < v[None,:]).mean(axis=0)*100
    d=targets[['index','description','assay','tf','cell']].copy()
    d['region']=regname; d['value']=v; d['bg_percentile']=pct
    d['bg_mean']=BG.mean(axis=0); d['bg_sd']=BG.std(axis=0)
    rows.append(d)
out=pd.concat(rows, ignore_index=True)
out.to_csv(f"{RES}/seqreg_enformer_promoter_tracks.csv.gz", index=False, compression="gzip")
print("wrote", len(out), "rows")

# a compact, readable summary for the focus factors and assays
foc=out[(out.tf.isin(['SP1','ETS1','RELB','JUND','JUN','FOS','FOSL1','FOSL2','SMAD1','SMAD2',
                      'SMAD5','CTCF','EP300','POLR2A','TEAD1','TEAD4','RUNX1','RUNX2','SRF',
                      'MYC','STAT3','TP53','NFIC','NFIB','KLF4','EGR1'])) |
        (out.assay.isin(['DNASE','ATAC','CAGE']))]
foc.to_csv(f"{RES}/seqreg_enformer_focus_tracks.csv", index=False)
print("focus rows", len(foc))
