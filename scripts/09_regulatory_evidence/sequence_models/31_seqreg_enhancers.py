#!/usr/bin/env python3
"""31_seqreg_enhancers.py -- candidate enhancers of COL1A1 and COL3A1: ENCODE cCREs within
+/-100 kb, intersected with ReMap 2022 peaks, reported by TF and by cell-type class, with
fibroblast/mesenchymal peaks called out. Also reports which cCREs carry the strongest
Enformer-predicted fibroblast DNase signal within the 114,688 bp Enformer output window.
"""
import numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
COLS=["chrom","start","end","name","score","strand","tstart","tend","rgb","tf","biotype"]
FIBRO={"BJ","BJ1-hTERT","GM23248","HDF","HFF","IMR-90","MRC-5","WI-38","WI-38VA13",
       "dermal-fibroblast","fibroblast","hMSC","hMSC-TERT","hMSC-TERT4","mesenchymal",
       "myofibroblast","primary-dermal-fibroblasts","primary-lung-fibroblast",
       "proliferating-human-fibroblast"}
peaks=pd.read_csv(f"{CACHE}/remap2022_hg38_loci.bed", sep="\t", header=None, names=COLS)
peaks["fibro"]=peaks.biotype.isin(FIBRO)
ccre=pd.read_csv(f"{RES}/seqreg_encode_ccre.csv")
reg=pd.read_csv(f"{RES}/seqreg_regions.csv").set_index("region")

# Enformer fibroblast DNase profile over the 114,688 bp output window
targets=pd.read_csv(f"{CACHE}/enformer_targets_human.txt", sep="\t")
targets['assay']=targets.description.str.split(':').str[0]
FIB_DNASE=targets.index[(targets.assay=='DNASE')&targets.description.str.contains('ibroblast')].to_numpy()
SEQ_LEN=196_608; N_BINS=896; BIN=128; CROP=(SEQ_LEN-N_BINS*BIN)//2

rows=[]
for g in ["COL1A1","COL3A1"]:
    r=reg.loc[g]; anchor=int(r.anchor)
    pred=np.load(f"{CACHE}/enformer/{g}_full.npy")
    fib=pred[:,FIB_DNASE].mean(axis=1)
    win_start=anchor-SEQ_LEN//2+CROP
    cc=ccre[(ccre.region==g)].drop_duplicates("accession")
    for _,c in cc.iterrows():
        b0=int((c.start-win_start)//BIN); b1=int((c.end-win_start)//BIN)+1
        if b1<=0 or b0>=N_BINS: enf=np.nan
        else: enf=float(fib[max(0,b0):min(N_BINS,b1)].max())
        sub=peaks[(peaks.chrom==c.chrom)&(peaks.start<c.end)&(peaks.end>c.start)]
        fib_sub=sub[sub.fibro]
        rows.append(dict(gene=g, accession=c.accession, ccre_class=c.ccre_class,
                         chrom=c.chrom, start=int(c.start), end=int(c.end),
                         dist_to_TSS=int(c.dist_to_anchor), zscore=c.zscore,
                         enformer_fibroblast_dnase_max=round(enf,4) if enf==enf else np.nan,
                         n_remap_peaks=len(sub), n_remap_tfs=sub.tf.nunique(),
                         n_fibro_peaks=len(fib_sub), n_fibro_tfs=fib_sub.tf.nunique(),
                         fibro_tfs=";".join(sorted(fib_sub.tf.unique()))[:400],
                         has_NFKB_family=int(sub.tf.isin(["NFKB1","NFKB2","RELA","RELB","REL"]).sum()),
                         has_SP=int(sub.tf.isin(["SP1","SP2","SP3"]).sum()),
                         has_ETS=int(sub.tf.isin(["ETS1","ETS2","ELK1","ELK4","GABPA"]).sum())))
E=pd.DataFrame(rows).sort_values(["gene","dist_to_TSS"])
E.to_csv(f"{RES}/seqreg_collagen_enhancers.csv", index=False)
pd.set_option("display.width",250)
for g in ["COL1A1","COL3A1"]:
    s=E[E.gene==g]
    print(f"\n===== {g}: {len(s)} cCREs within +/-100 kb; {int(s.n_fibro_peaks.sum())} fibroblast ReMap peaks")
    print("-- top 12 cCREs by Enformer-predicted fibroblast DNase:")
    print(s.nlargest(12,"enformer_fibroblast_dnase_max")[["accession","ccre_class","dist_to_TSS",
          "enformer_fibroblast_dnase_max","n_remap_tfs","n_fibro_tfs","has_NFKB_family","has_SP","has_ETS","fibro_tfs"]].to_string(index=False))
    print("-- top 8 cCREs by number of fibroblast TFs bound:")
    print(s.nlargest(8,"n_fibro_tfs")[["accession","ccre_class","dist_to_TSS","n_fibro_tfs","fibro_tfs"]].to_string(index=False))
