#!/usr/bin/env python3
"""22b_seqreg_enformer_ism.py -- sequence attribution at the COL1A1 and COL3A1 promoters
with Enformer.

Two complementary experiments, both purely in silico and both reported as changes relative
to the reference sequence:
 (A) single-base in-silico mutagenesis (ISM) across TSS-400..TSS+100, all 3 alternative
     bases at every position, scored on (i) mean predicted CAGE over the 25 fibroblast CAGE
     tracks and (ii) mean predicted DNase over the 31 fibroblast DNase tracks, at the
     promoter bins;
 (B) motif occlusion: every JASPAR/HOCOMOCO site called at p<1e-4 in the promoter by
     28_seqreg_motifscan.py is replaced, in place, by 20 independent dinucleotide-preserving
     shuffles of that site; the reported effect is the mean change over the 20 shuffles.
     This asks whether the model treats a given NFKB/SP1/ETS1 site as load-bearing.
"""
import os, time, itertools, random
import numpy as np, pandas as pd, torch
from enformer_pytorch import from_pretrained

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
SEQ_LEN=196_608; N_BINS=896; BIN=128
CROP=(SEQ_LEN-N_BINS*BIN)//2; CENTRE_BIN=(SEQ_LEN//2-CROP)//BIN
HALF=300_000
ISM_LO, ISM_HI = -400, 100         # relative to TSS, in transcription orientation

targets=pd.read_csv(f"{CACHE}/enformer_targets_human.txt", sep="\t")
targets['assay']=targets.description.str.split(':').str[0]
FIB_CAGE=targets.index[(targets.assay=='CAGE')&targets.description.str.contains('ibroblast')].to_numpy()
FIB_DNASE=targets.index[(targets.assay=='DNASE')&targets.description.str.contains('ibroblast')].to_numpy()
ALL_CAGE=targets.index[targets.assay=='CAGE'].to_numpy()
print("fibroblast CAGE tracks:",len(FIB_CAGE)," fibroblast DNase tracks:",len(FIB_DNASE))

MAP={'A':0,'C':1,'G':2,'T':3}
def onehot(s):
    a=np.zeros((len(s),4),dtype=np.float32)
    idx=np.frombuffer(s.encode(),dtype=np.uint8)
    for b,i in MAP.items(): a[idx==ord(b),i]=1.0
    return a
def revcomp(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

model=from_pretrained('EleutherAI/enformer-official-rough').eval().cuda()
BATCH=2
@torch.no_grad()
def predict_batch(seqs):
    outs=[]
    for k in range(0,len(seqs),BATCH):
        x=torch.from_numpy(np.stack([onehot(s) for s in seqs[k:k+BATCH]])).cuda()
        with torch.autocast('cuda', dtype=torch.bfloat16):
            o=model(x)['human']
        outs.append(o.float().cpu().numpy())
        del x,o
    torch.cuda.empty_cache()
    return np.concatenate(outs,axis=0)
def score(pred):
    prom=pred[:,CENTRE_BIN-1:CENTRE_BIN+2,:].mean(axis=1)
    return dict(fib_cage=prom[:,FIB_CAGE].mean(axis=1),
                fib_dnase=prom[:,FIB_DNASE].mean(axis=1),
                all_cage=prom[:,ALL_CAGE].mean(axis=1))

reg=pd.read_csv(f"{RES}/seqreg_regions.csv").set_index("region")
def genomic_window(regname):
    r=reg.loc[regname]; a=int(r.anchor)
    s=open(f"{CACHE}/seq/{r.chrom}_{a-HALF}_{a+HALF}.txt").read().strip()
    return s[HALF-SEQ_LEN//2:HALF+SEQ_LEN//2], r.strand, r.chrom, a

def rel_to_index(rel, strand):
    """rel is relative to the TSS in transcription orientation; return index in the
       plus-strand 196,608 bp window (whose centre index SEQ_LEN//2 is the TSS)."""
    return SEQ_LEN//2 + (rel if strand=='+' else -rel)

BASES="ACGT"
rng=random.Random(20260909)
def dinuc_shuffle(s):
    # Altschul-Erikson dinucleotide shuffle (simple edge-ordering implementation)
    if len(s)<3: return s
    for _ in range(50):
        edges={}
        for a,b in zip(s,s[1:]): edges.setdefault(a,[]).append(b)
        for k in edges: rng.shuffle(edges[k])
        out=[s[0]]; cur=s[0]; ok=True
        for _ in range(len(s)-1):
            if not edges.get(cur): ok=False; break
            nxt=edges[cur].pop(); out.append(nxt); cur=nxt
        if ok and len(out)==len(s): return "".join(out)
    return "".join(rng.sample(s,len(s)))

ism_rows=[]; occ_rows=[]
sites=pd.read_csv(f"{RES}/seqreg_motif_sites.csv") if os.path.exists(f"{RES}/seqreg_motif_sites.csv") else pd.DataFrame()
# occlusion is restricted to the NF-kB / SP / ETS focus families: occluding all 618-1008
# called sites per promoter would need >12,000 forward passes per locus.
if len(sites):
    sites["tf_name"]=sites.motif.str.split("|").str[-1].str.split("_HUMAN").str[0]
    FOCUS_TFS={"NFKB1","NFKB2","RELA","RELB","REL","SP1","SP2","SP3",
               "ETS1","ETS2","ELK1","ELK4","GABPA"}
    sites=sites[sites.tf_name.isin(FOCUS_TFS)].copy()
    print("focus sites for occlusion:", len(sites), sites.groupby("region").size().to_dict(), flush=True)

for regname in ["COL1A1","COL3A1"]:
    ref, strand, chrom, anchor = genomic_window(regname)
    base=predict_batch([ref]); b=score(base)
    print(f"{regname} reference: fib_cage={b['fib_cage'][0]:.3f} fib_dnase={b['fib_dnase'][0]:.3f}", flush=True)

    # ---- (A) ISM
    if regname in ("COL1A1","COL3A1"):
        muts=[]; meta=[]
        for rel in range(ISM_LO, ISM_HI):
            i=rel_to_index(rel, strand)
            wt=ref[i]
            for alt in BASES:
                if alt==wt: continue
                muts.append(ref[:i]+alt+ref[i+1:])
                meta.append((rel, i, wt, alt))
        t0=time.time(); out=[]
        for k in range(0,len(muts),16):
            out.append(score(predict_batch(muts[k:k+16])))
            if (k//16)%10==0: print(f"  ISM {k}/{len(muts)} {time.time()-t0:.0f}s", flush=True)
        S={k:np.concatenate([o[k] for o in out]) for k in out[0]}
        for j,(rel,i,wt,alt) in enumerate(meta):
            ism_rows.append(dict(region=regname, rel_to_TSS=rel, genome_pos=anchor+(rel if strand=='+' else -rel),
                                 ref_base=wt, alt_base=alt,
                                 d_fib_cage=float(S['fib_cage'][j]-b['fib_cage'][0]),
                                 d_fib_dnase=float(S['fib_dnase'][j]-b['fib_dnase'][0]),
                                 d_all_cage=float(S['all_cage'][j]-b['all_cage'][0])))
        pd.DataFrame(ism_rows).to_csv(f"{RES}/seqreg_enformer_ism.csv", index=False)

    # ---- (B) motif occlusion
    if len(sites):
        ss=sites[sites.region==regname]
        for _,r in ss.iterrows():
            rel=int(r.rel_to_TSS); L=len(str(r.seq))
            i=rel_to_index(rel, strand)
            if strand=='+': lo,hi=i,i+L
            else:           lo,hi=i-L+1,i+1
            wt=ref[lo:hi]
            shuffles=[]
            for _ in range(20):
                sh=dinuc_shuffle(wt if strand=='+' else revcomp(wt))
                if strand=='-': sh=revcomp(sh)
                shuffles.append(ref[:lo]+sh+ref[hi:])
            o=[]
            for k in range(0,len(shuffles),10): o.append(score(predict_batch(shuffles[k:k+10])))
            S={k:np.concatenate([x[k] for x in o]) for k in o[0]}
            occ_rows.append(dict(region=regname, db=r.db, motif=r.motif, rel_to_TSS=rel,
                                 motif_strand=r.strand, motif_score=r.score, motif_p=r.pvalue,
                                 site_seq=r.seq,
                                 ref_fib_cage=float(b['fib_cage'][0]),
                                 mean_shuffled_fib_cage=float(S['fib_cage'].mean()),
                                 d_fib_cage=float(S['fib_cage'].mean()-b['fib_cage'][0]),
                                 sd_fib_cage=float(S['fib_cage'].std()),
                                 d_fib_dnase=float(S['fib_dnase'].mean()-b['fib_dnase'][0]),
                                 d_all_cage=float(S['all_cage'].mean()-b['all_cage'][0])))
        pd.DataFrame(occ_rows).to_csv(f"{RES}/seqreg_enformer_motif_occlusion.csv", index=False)
        print(f"  occlusion done for {regname}: {len(ss)} sites", flush=True)
print("ISM/occlusion complete")
