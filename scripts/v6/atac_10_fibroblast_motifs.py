#!/usr/bin/env python3
"""(b) In the compartment where the collagens are actually transcribed:
do ETS1/NFKB1/RELA/SP1 motifs occur in the fibroblast-accessible COL1A1 and
COL3A1 promoters?  ENCODE DNase-seq, fibroblast of mammary gland (ENCSR000ENW,
file ENCFF032AZY, hg38).  Peaks recentred to 501 bp.  Null = GC-decile matched
fibroblast peaks, 10,000 resamples; BH across gene x motif tests."""
import os, subprocess, tempfile, numpy as np, pandas as pd
import importlib.util
ROOT="/path/to/revision"; C=f"{ROOT}/cache/v6/atac"; R=f"{ROOT}/results/v6"
spec=importlib.util.spec_from_file_location("ms",f"{ROOT}/scripts/v6/atac_03_motif_scan.py")
ms=importlib.util.module_from_spec(spec); spec.loader.exec_module(ms)
rng=np.random.default_rng(7)
FILES={"fibroblast_of_mammary_gland_DNase":"ENCFF032AZY","MCF7_ATAC":"ENCFF821OEF"}
PANEL={"ETS1":"MA0098.4","NFKB1":"MA0105.4","RELA":"MA0107.1","SP1":"MA0079.5",
       "FOS::JUN":"MA0099.4","TEAD1":"MA0090.3","RUNX2":"MA0511.2","CTCF":"MA0139.2",
       "SMAD3":"MA0795.1","TWIST1":"MA1123.2","PRRX1":"MA0716.2","NFYA":"MA0060.4"}
mot=ms.read_jaspar(f"{C}/JASPAR2024_CORE_vert_nr_pfms.txt")
PANEL={k:v for k,v in PANEL.items() if v in mot}
tss=pd.read_csv(f"{C}/refseq_select_tss.tsv",sep="\t").drop_duplicates("symbol").set_index("symbol")
GENES=["COL1A1","COL3A1","COL1A2","DCN","LUM","POSTN","FAP","EPCAM","ESR1","GAPDH"]
tmp=tempfile.mkdtemp(); out_rows=[]
for label,fid in FILES.items():
    pk=pd.read_csv(f"{C}/encode/{fid}.bed.gz",sep="\t",header=None,usecols=[0,1,2],names=["chrom","start","end"])
    pk=pk[pk.chrom.str.match(r"^chr[0-9XY]+$")].drop_duplicates().copy()   # ENCODE beds repeat rows
    ctr=((pk.start+pk.end)//2).values
    pk["s"]=np.maximum(0,ctr-250); pk["e"]=ctr+251
    bed=f"{tmp}/{fid}.bed"; pk[["chrom","s","e"]].assign(n=[f"p{i}" for i in range(len(pk))]).to_csv(bed,sep="\t",header=False,index=False)
    fa=f"{tmp}/{fid}.tab"
    subprocess.run(["bedtools","getfasta","-fi",f"{ROOT}/cache/seqreg/genome/hg38.fa","-bed",bed,"-name","-tab","-fo",fa],check=True)
    names=[];seqs=[]
    for line in open(fa):
        a,b=line.rstrip("\n").split("\t"); names.append(a.split("::")[0]); seqs.append(b.upper())
    L=501; keep=[i for i,s in enumerate(seqs) if len(s)==L]
    pk=pk.iloc[keep].reset_index(drop=True); seqs=[seqs[i] for i in keep]
    code=np.full(256,4,dtype=np.int8)
    for i,c in enumerate("ACGT"): code[ord(c)]=i
    S=code[np.frombuffer("".join(seqs).encode(),dtype=np.uint8)].reshape(len(seqs),L)
    cnt=np.bincount(S.ravel(),minlength=5)[:4].astype(float); bg=cnt/cnt.sum()
    gc=((S==1)|(S==2)).mean(1)
    print(f"{label}: {len(seqs)} peaks, bg ACGT {np.round(bg,3)}",flush=True)
    H=np.zeros((len(seqs),len(PANEL)),dtype=np.int16); thrs={}
    for mi,(nm,mid) in enumerate(PANEL.items()):
        _,counts=mot[mid]; lo=ms.logodds(counts,bg); thr=ms.score_threshold(lo,bg); thrs[nm]=thr
        w=lo.shape[0]; lo_rc=lo[::-1,::-1]; npos=L-w+1
        for st in range(0,len(seqs),20000):
            en=min(st+20000,len(seqs)); sub=S[st:en]
            f=np.zeros((en-st,npos),dtype=np.float32); r_=np.zeros_like(f)
            for j in range(w):
                col=sub[:,j:j+npos]
                f+=lo[j][np.minimum(col,3)]*(col<4)+(col==4)*(-1e4)
                r_+=lo_rc[j][np.minimum(col,3)]*(col<4)+(col==4)*(-1e4)
            H[st:en,mi]=(f>=thr).sum(1)+(r_>=thr).sum(1)
    np.save(f"{C}/encode_motifhits_{fid}.npy",H)
    gcq=pd.qcut(pd.Series(gc),10,labels=False,duplicates="drop").values
    for g in GENES:
        if g not in tss.index: continue
        t=tss.loc[g]
        # any accessible element whose ORIGINAL interval overlaps TSS +/- 1 kb
        sel=np.where((pk.chrom.values==t.chrom)&(pk.end.values>t.tss-1000)&(pk.start.values<t.tss+1000))[0]
        for mi,(nm,mid) in enumerate(PANEL.items()):
            if len(sel)==0:
                out_rows.append(dict(dataset=label,gene=g,motif=nm,n_promoter_peaks=0,obs_with_hit=np.nan,
                                     frac_exp=np.nan,p_emp=np.nan)); continue
            obs=int((H[sel,mi]>0).sum()); draws=np.zeros(10000,dtype=np.int32)
            for s in np.unique(gcq[sel]):
                pool=np.where(gcq==s)[0]; k=int((gcq[sel]==s).sum())
                v=(H[pool,mi]>0).astype(np.int8)
                draws+=v[rng.integers(0,len(v),size=(10000,k))].sum(1)
            exp=draws.mean(); phi=(draws>=obs).mean(); plo=(draws<=obs).mean()
            out_rows.append(dict(dataset=label,gene=g,motif=nm,n_promoter_peaks=len(sel),
                                 obs_with_hit=obs,frac_obs=obs/len(sel),exp_with_hit=exp,frac_exp=exp/len(sel),
                                 p_emp=min(1.0,2*min(phi,plo)),threshold_bits=thrs[nm]))
D=pd.DataFrame(out_rows)
m=D.p_emp.notna(); p=D.loc[m,"p_emp"].values; o=np.argsort(p); q=np.empty(len(p))
q[o]=np.minimum.accumulate((p[o]*len(p)/(np.arange(len(p))+1))[::-1])[::-1]
D.loc[m,"fdr_BH"]=np.clip(q,0,1)
D.to_csv(f"{R}/atac_fibroblast_promoter_motifs.csv",index=False)
pd.set_option("display.width",240)
print(D[D.gene.isin(["COL1A1","COL3A1"])].to_string(index=False))
print("DONE")
