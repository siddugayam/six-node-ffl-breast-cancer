#!/usr/bin/env python3
"""35_seqreg_nfkb_context.py -- every NF-kB-family ReMap 2022 peak within +/-100 kb of COL1A1
and COL3A1, annotated by cell type and by the treatment string encoded in the ReMap dataset
name (TNF, LPS, IL1, etc.). The NF-kB -> collagen edges of the network concern resting tumour
tissue; this asks in what cellular context NF-kB has been observed at these loci.
"""
import re, numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
COLS=["chrom","start","end","name","score","strand","tstart","tend","rgb","tf","biotype"]
p=pd.read_csv(f"{CACHE}/remap2022_hg38_loci.bed", sep="\t", header=None, names=COLS)
reg=pd.read_csv(f"{RES}/seqreg_regions.csv").set_index("region")
FIBRO={"BJ","BJ1-hTERT","GM23248","HDF","HFF","IMR-90","MRC-5","WI-38","WI-38VA13",
       "dermal-fibroblast","fibroblast","hMSC","hMSC-TERT","hMSC-TERT4","mesenchymal",
       "myofibroblast","primary-dermal-fibroblasts","primary-lung-fibroblast",
       "proliferating-human-fibroblast"}
MESEN_LINEAGE=FIBRO|{"SGBS","myoblast","osteoblast","primary-chondrocyte","preadipocyte"}
STIM=r"(?i)(TNF|LPS|IL1|IL-1|PMA|TPA|CD40|ETOPOSIDE|IR|Gy|NUTLIN|SIFOXA1|POLYIC|IFN|LTA|SENDAI|VIRUS)"
rows=[]
for g in ["COL1A1","COL3A1"]:
    r=reg.loc[g]; a=int(r.anchor)
    sub=p[(p.chrom==r.chrom)&(p.start<a+100000)&(p.end>a-100000)&
          (p.tf.isin(["NFKB1","NFKB2","RELA","RELB","REL"]))].copy()
    sub["mid"]=(sub.start+sub.end)//2
    sub["dist_to_TSS"]=(sub.mid-a)*(1 if r.strand=='+' else -1)
    sub["dataset"]=sub["name"].str.split(".").str[0]
    sub["condition"]=sub["name"].str.split(".").str[2]
    sub["stimulated"]=sub.condition.str.contains(STIM, regex=True, na=False)
    sub["fibroblast"]=sub.biotype.isin(FIBRO)
    sub["mesenchymal_lineage"]=sub.biotype.isin(MESEN_LINEAGE)
    sub["gene"]=g
    sub["window"]=np.where(sub.dist_to_TSS.abs()<=1000,"promoter_1kb",
                    np.where(sub.dist_to_TSS.abs()<=10000,"proximal_10kb","distal_100kb"))
    rows.append(sub[["gene","tf","biotype","dataset","condition","stimulated","fibroblast",
                     "mesenchymal_lineage","chrom","start","end","dist_to_TSS","window"]])
D=pd.concat(rows, ignore_index=True)
D.to_csv(f"{RES}/seqreg_nfkb_peaks_context.csv", index=False)
pd.set_option("display.width",260)
print("Total NF-kB-family peaks within +/-100 kb:", len(D))
print(D.groupby(["gene","tf"]).size().to_string())
print("\nby window and stimulation:")
print(pd.crosstab([D.gene,D.window],[D.stimulated]).to_string())
print("\nmesenchymal-lineage peaks (fibroblast whitelist + SGBS preadipocyte, myoblast, osteoblast, chondrocyte):")
m=D[D.mesenchymal_lineage]
print(m[["gene","tf","biotype","dataset","condition","stimulated","dist_to_TSS","window"]].sort_values(["gene","dist_to_TSS"]).to_string(index=False))
print("\nsummary counts:")
S=[]
for g in ["COL1A1","COL3A1"]:
    d=D[D.gene==g]
    S.append(dict(gene=g, total_peaks=len(d), n_biotypes=d.biotype.nunique(),
                  promoter_peaks=int((d.window=="promoter_1kb").sum()),
                  promoter_fibroblast=int(((d.window=="promoter_1kb")&d.fibroblast).sum()),
                  promoter_mesenchymal=int(((d.window=="promoter_1kb")&d.mesenchymal_lineage).sum()),
                  fibroblast_peaks=int(d.fibroblast.sum()),
                  mesenchymal_peaks=int(d.mesenchymal_lineage.sum()),
                  stimulated_fraction=round(float(d.stimulated.mean()),3),
                  fibroblast_all_stimulated=bool(d[d.fibroblast].stimulated.all()) if d.fibroblast.any() else None))
Sd=pd.DataFrame(S); Sd.to_csv(f"{RES}/seqreg_nfkb_context_summary.csv", index=False)
print(Sd.to_string(index=False))
