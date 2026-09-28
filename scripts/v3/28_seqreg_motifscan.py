#!/usr/bin/env python3
"""28_seqreg_motifscan.py -- FIMO-equivalent PWM scan of the COL1A1 / COL3A1 (and miRNA-locus)
promoters with JASPAR 2024 CORE vertebrates non-redundant and HOCOMOCO v11 core human,
calibrated against a GC-matched genome-wide promoter background.

Method (matches MEME-FIMO): PFM -> log-odds PWM against a 0-order background estimated from
1,000 random human promoters; exact p-value per site by dynamic programming on the discretised
score distribution; site threshold p < 1e-4; Benjamini-Hochberg q within each promoter/TF.
Because SP1-type GC-boxes are ubiquitous in GC-rich promoters, every hit count is additionally
reported as an empirical percentile against GC-matched random promoters.
"""
import os, re, gzip, time
import numpy as np, pandas as pd

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
PROM_UP, PROM_DN = 1000, 500          # promoter = TSS-1000 .. TSS+500 (strand aware)
NBG = int(os.environ.get("NBG_MOTIF","1000"))
PTHRESH = 1e-4
rng=np.random.default_rng(20260909)

# ---------- motif parsing ----------------------------------------------------
def read_jaspar(path):
    ms={}; name=None; rows=[]
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name: ms[name]=np.array(rows,dtype=float)
            parts=line[1:].split()
            name=f"{parts[0]}|{parts[1]}"; rows=[]
        elif line.strip():
            v=[float(x) for x in re.findall(r"[-\d.]+", line.split("[")[-1])]
            rows.append(v)
    if name: ms[name]=np.array(rows,dtype=float)
    return {k:v.T for k,v in ms.items()}      # -> (L,4) counts, order ACGT

def read_hocomoco(path):
    ms={}; name=None; rows=[]
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name: ms[name]=np.array(rows,dtype=float)
            name=line[1:].strip(); rows=[]
        elif line.strip():
            rows.append([float(x) for x in line.split()])
    if name: ms[name]=np.array(rows,dtype=float)
    return ms                                  # already log-odds (L,4) vs uniform

JAS=read_jaspar(f"{CACHE}/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
HOC=read_hocomoco(f"{CACHE}/HOCOMOCOv11_core_pwms_HUMAN_mono.txt")
print("JASPAR", len(JAS), "HOCOMOCO", len(HOC))

# ---------- sequences --------------------------------------------------------
MAPI={'A':0,'C':1,'G':2,'T':3}
def encode(s):
    a=np.full(len(s),4,dtype=np.int8)
    b=np.frombuffer(s.encode(),dtype=np.uint8)
    for c,i in MAPI.items(): a[b==ord(c)]=i
    return a
def revcomp(s):
    return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

reg=pd.read_csv(f"{RES}/seqreg_regions.csv"); HALF=300_000
def promoter(chrom, anchor, strand):
    s=open(f"{CACHE}/seq/{chrom}_{anchor-HALF}_{anchor+HALF}.txt").read().strip()
    c=HALF
    if strand=='+': sub=s[c-PROM_UP : c+PROM_DN]
    else:           sub=revcomp(s[c-PROM_DN : c+PROM_UP])
    return sub

targets={}
for _,r in reg.iterrows():
    targets[r.region]=promoter(r.chrom,int(r.anchor),r.strand)

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
bgsel=rs.sample(n=NBG*2, random_state=1)
bgseqs=[]; bgnames=[]
for _,r in bgsel.iterrows():
    if len(bgseqs)>=NBG: break
    t=int(r.tss)
    if t-PROM_UP<0 or t+PROM_UP>len(fa[r.chrom]): continue
    s=str(fa[r.chrom][t-PROM_UP:t+PROM_DN]) if r.strand=="+" else revcomp(str(fa[r.chrom][t-PROM_DN:t+PROM_UP]))
    if s.count('N')>0 or len(s)!=PROM_UP+PROM_DN: continue
    bgseqs.append(s); bgnames.append(r.name2)
print("background promoters:", len(bgseqs))

# 0-order background from the background promoter set
allbg="".join(bgseqs)
cnt=np.array([allbg.count(c) for c in "ACGT"],dtype=float); BG0=cnt/cnt.sum()
print("background base composition ACGT:", np.round(BG0,4))

# ---------- PWM construction and exact p-values ------------------------------
def pfm_to_pwm(pfm, bg, pseudo=0.1):
    n=pfm.sum(axis=1, keepdims=True); n[n==0]=1
    p=(pfm+pseudo*bg[None,:])/(n+pseudo)
    return np.log2(p/bg[None,:])

def hoc_to_pwm(lw, bg):
    # HOCOMOCO PWMs are log-odds vs uniform 0.25; convert probabilities then re-score vs bg
    p=np.exp(lw)*0.25
    p=p/p.sum(axis=1,keepdims=True)
    p=np.clip(p,1e-6,None); p=p/p.sum(axis=1,keepdims=True)
    return np.log2(p/bg[None,:])

SCALE=100
def score_threshold(pwm, bg, pval):
    """Exact FIMO-style: discretise, DP the null distribution, return score cut-off."""
    q=np.round(pwm*SCALE).astype(int)
    off=q.min(axis=1); q=q-off[:,None]
    maxs=int(q.max(axis=1).sum())
    dist=np.zeros(maxs+1); dist[0]=1.0
    for l in range(q.shape[0]):
        nd=np.zeros(maxs+1)
        for j in range(4):
            nd[q[l,j]:q[l,j]+len(dist)-q[l,j]] += bg[j]*dist[:len(dist)-q[l,j]]
        dist=nd
    tail=np.cumsum(dist[::-1])[::-1]          # tail[i] = P(shifted score >= i)
    ok=np.nonzero(tail<=pval)[0]
    cut_shift=int(ok[0]) if len(ok) else int(maxs)
    return (cut_shift+off.sum())/SCALE, tail, int(off.sum())

def scan_counts(seqs_idx, pwm, cutoff):
    """count sites >= cutoff on both strands over a concatenated int-encoded array,
       with separators (value 4) preventing cross-sequence hits."""
    L=pwm.shape[0]; M=len(seqs_idx)
    if M<L: return 0, np.array([])
    lut=np.concatenate([pwm, np.full((L,1), -1e6)], axis=1)   # N -> -inf
    sc=np.zeros(M-L+1)
    for l in range(L):
        sc += lut[l][seqs_idx[l:M-L+1+l]]
    pwm_rc=pwm[::-1,::-1]
    lut2=np.concatenate([pwm_rc, np.full((L,1), -1e6)], axis=1)
    sc2=np.zeros(M-L+1)
    for l in range(L):
        sc2 += lut2[l][seqs_idx[l:M-L+1+l]]
    hits=np.where((sc>=cutoff)|(sc2>=cutoff))[0]
    best=np.maximum(sc,sc2)
    return len(hits), best

SEP=np.full(30,4,dtype=np.int8)
bg_concat=np.concatenate([np.concatenate([encode(s),SEP]) for s in bgseqs])
tgt_idx={k:encode(v) for k,v in targets.items()}
bg_gc=np.array([(s.count('G')+s.count('C'))/len(s) for s in bgseqs])

motifs={}
for k,v in JAS.items(): motifs[("JASPAR",k)]=pfm_to_pwm(v,BG0)
for k,v in HOC.items(): motifs[("HOCOMOCO",k)]=hoc_to_pwm(v,BG0)
print("total motifs:", len(motifs))

rows=[]; sites=[]
t0=time.time()
for i,((db,name),pwm) in enumerate(motifs.items()):
    if pwm.shape[0]>30: pwm=pwm[:30]
    cut,tail,offsum=score_threshold(pwm,BG0,PTHRESH)
    nbg,_=scan_counts(bg_concat,pwm,cut)
    bg_rate=nbg/ (len(bg_concat)-pwm.shape[0]+1)     # sites per scanned position (2 strands)
    for regname,idx in tgt_idx.items():
        n,best=scan_counts(idx,pwm,cut)
        exp = bg_rate*(len(idx)-pwm.shape[0]+1)
        rows.append(dict(db=db, motif=name, region=regname, pwm_len=pwm.shape[0],
                         score_cutoff=round(cut,4), n_sites=n, best_score=round(float(best.max()),4) if len(best) else np.nan,
                         expected_sites_bg=round(exp,4),
                         enrich=round(n/exp,3) if exp>0 else np.nan,
                         bg_sites_total=nbg))
        if n>0:
            L=pwm.shape[0]
            lut=np.concatenate([pwm, np.full((L,1),-1e6)],axis=1)
            sc=np.zeros(len(idx)-L+1)
            for l in range(L): sc+=lut[l][idx[l:len(idx)-L+1+l]]
            pr=pwm[::-1,::-1]; lut2=np.concatenate([pr,np.full((L,1),-1e6)],axis=1)
            sc2=np.zeros(len(idx)-L+1)
            for l in range(L): sc2+=lut2[l][idx[l:len(idx)-L+1+l]]
            for pos in np.where((sc>=cut)|(sc2>=cut))[0]:
                strand='+' if sc[pos]>=sc2[pos] else '-'
                s=float(max(sc[pos],sc2[pos]))
                pv=float(tail[int(np.clip(round(s*SCALE)-offsum,0,len(tail)-1))])
                sites.append(dict(db=db, motif=name, region=regname, offset_in_promoter=int(pos),
                                  rel_to_TSS=int(pos)-PROM_UP, strand=strand, score=round(s,4),
                                  pvalue=pv, seq=targets[regname][pos:pos+L]))
    if (i+1)%100==0: print(f"  {i+1}/{len(motifs)} {time.time()-t0:.0f}s", flush=True)

df=pd.DataFrame(rows); df.to_csv(f"{RES}/seqreg_motif_scan_counts.csv", index=False)
sd=pd.DataFrame(sites); sd.to_csv(f"{RES}/seqreg_motif_sites.csv", index=False)
print("counts", len(df), "sites", len(sd))

# ---------- per-promoter empirical percentile for the focus TFs --------------
FOCUS_PAT = r"(?i)(^|\|)(NFKB1|NFKB2|RELA|RELB|REL|SP1|SP2|SP3|ETS1|ETS2|ELK1|ELK4|GABPA|SMAD3|SMAD4|TWIST1|RUNX1|RUNX2|TEAD1|SRF|MKL|JUN|JUND|FOS|FOSL1|FOSL2|EGR1|KLF4|CTCF|STAT3|TP53|MYC|NFIC|NFIB|SOX9|ZEB1|SNAI2)($|_)"
foc={k:v for k,v in motifs.items() if re.search(FOCUS_PAT, k[1].split("|")[-1])}
print("focus motifs:", len(foc), sorted(set(k[1] for k in foc))[:40])

perc=[]
for (db,name),pwm in foc.items():
    if pwm.shape[0]>30: pwm=pwm[:30]
    cut,tail,offsum=score_threshold(pwm,BG0,PTHRESH)
    # per-background-promoter counts
    counts=[]
    for s in bgseqs:
        n,_=scan_counts(encode(s),pwm,cut); counts.append(n)
    counts=np.array(counts)
    for regname,idx in tgt_idx.items():
        n,_=scan_counts(idx,pwm,cut)
        gc=(targets[regname].count('G')+targets[regname].count('C'))/len(targets[regname])
        sel=np.abs(bg_gc-gc)<=0.05
        if sel.sum()<50: sel=np.abs(bg_gc-gc)<=0.10
        perc.append(dict(db=db, motif=name, region=regname, n_sites=n, promoter_gc=round(gc,4),
                         bg_all_mean=round(counts.mean(),3),
                         pct_all=round(float((counts<n).mean()*100),2),
                         pct_all_ge=round(float((counts>=n).mean()*100),3),
                         n_gc_matched=int(sel.sum()),
                         bg_gcmatched_mean=round(float(counts[sel].mean()),3),
                         pct_gcmatched=round(float((counts[sel]<n).mean()*100),2),
                         emp_p_gcmatched=round(float((counts[sel]>=n).mean()),4)))
    print("  focus", name, flush=True)
pd.DataFrame(perc).to_csv(f"{RES}/seqreg_motif_focus_calibrated.csv", index=False)
print("done")
