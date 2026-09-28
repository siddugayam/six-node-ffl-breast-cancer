#!/usr/bin/env python3
"""25b -- recompute the ReMap 2022 region tables with a CURATED cell-type (biotype) whitelist.
The regex used in 25_ matched BJAB (B-cell lymphoma), Detroit-562 (pharyngeal carcinoma),
HS578T (breast carcinoma) and hTERT-HME1 (mammary epithelial) as 'fibroblast-like'; those are
excluded here. Whitelists were curated by eye from the 721 biotypes present at these loci.
"""
import pandas as pd, numpy as np, os
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
COLS=["chrom","start","end","name","score","strand","tstart","tend","rgb","tf","biotype"]

FIBRO={"BJ","BJ1-hTERT","GM23248","HDF","HFF","IMR-90","MRC-5","WI-38","WI-38VA13",
       "dermal-fibroblast","fibroblast","hMSC","hMSC-TERT","hMSC-TERT4","mesenchymal",
       "myofibroblast","primary-dermal-fibroblasts","primary-lung-fibroblast",
       "proliferating-human-fibroblast"}
MESEN_OTHER={"myoblast","osteoblast","primary-chondrocyte"}

p=pd.read_csv(f"{CACHE}/remap2022_hg38_loci.bed", sep="\t", header=None, names=COLS)
bt=set(p.biotype.unique())
BREAST={b for b in bt if any(k in b.upper() for k in
        ["MCF","T-47D","T47D","MDA-MB","SK-BR","SKBR","ZR-75","ZR75","BT-474","BT474",
         "HCC1","HS578T","HME","MAMMARY","BREAST","CAMA","UACC","LM2","SUM1"])}
p["celltype_class"]=np.where(p.biotype.isin(FIBRO),"fibroblast_mesenchymal",
                     np.where(p.biotype.isin(BREAST),"breast",
                      np.where(p.biotype.isin(MESEN_OTHER),"mesenchymal_other","other")))
pd.Series(sorted(FIBRO)).to_csv(f"{RES}/seqreg_remap_biotype_whitelist_fibroblast.csv",index=False,header=["biotype"])
pd.Series(sorted(BREAST)).to_csv(f"{RES}/seqreg_remap_biotype_whitelist_breast.csv",index=False,header=["biotype"])
print("fibroblast biotypes present:",sorted(FIBRO & bt))
print("breast biotypes present:",sorted(BREAST))

reg=pd.read_csv(f"{RES}/seqreg_regions.csv")
def window(r,up,dn):
    a=int(r.anchor); return (a-up,a+dn) if r.strand=='+' else (a-dn,a+up)

rows=[]
for _,r in reg.iterrows():
    for label,(up,dn) in {"prom_1kb":(1000,500),"prom_5kb":(5000,2000),
                          "prom_10kb":(10000,5000),"locus_100kb":(100000,100000)}.items():
        s,e=window(r,up,dn)
        sub=p[(p.chrom==r.chrom)&(p.start<e)&(p.end>s)]
        for tf,g in sub.groupby("tf"):
            cc=g.celltype_class.value_counts()
            rows.append(dict(region=r.region, window=label, tf=tf, n_peaks=len(g),
                n_biotypes=g.biotype.nunique(),
                n_peaks_fibroblast=int(cc.get("fibroblast_mesenchymal",0)),
                n_peaks_breast=int(cc.get("breast",0)),
                fibroblast_biotypes=";".join(sorted(set(g.loc[g.celltype_class=="fibroblast_mesenchymal","biotype"]))),
                breast_biotypes=";".join(sorted(set(g.loc[g.celltype_class=="breast","biotype"])))[:300],
                all_biotypes=";".join(sorted(set(g.biotype))[:30])))
out=pd.DataFrame(rows)
out.to_csv(f"{RES}/seqreg_remap_tf_by_region_curated.csv", index=False)

# experiment-level table restricted to fibroblast/mesenchymal and breast biotypes
sel=p[p.celltype_class.isin(["fibroblast_mesenchymal","breast"])].copy()
recs=[]
for _,r in reg.iterrows():
    for label,(up,dn) in {"prom_1kb":(1000,500),"prom_5kb":(5000,2000),"prom_10kb":(10000,5000)}.items():
        s,e=window(r,up,dn)
        sub=sel[(sel.chrom==r.chrom)&(sel.start<e)&(sel.end>s)]
        for _,q in sub.iterrows():
            recs.append(dict(region=r.region, window=label, tf=q.tf, biotype=q.biotype,
                             celltype_class=q.celltype_class, dataset=q["name"].split(".")[0],
                             peak_start=q.start, peak_end=q.end))
pd.DataFrame(recs).to_csv(f"{RES}/seqreg_remap_fibroblast_breast_peaks.csv", index=False)
print("curated rows:",len(out)," fibro/breast peak rows:",len(recs))

FOC=["NFKB1","NFKB2","RELA","RELB","REL","SP1","SP2","SP3","ETS1","ETS2","ELK1","ELK4","GABPA",
     "SMAD3","SMAD4","TWIST1","RUNX1","RUNX2","TEAD1","SRF","JUN","JUND","FOS","FOSL1","FOSL2",
     "STAT3","EP300","MED1","CTCF","YY1","KLF4","EGR1","ESR1","MYC","TP53","HIF1A","SOX9","ZEB1"]
f=out[(out.region.isin(["COL1A1","COL3A1"]))&(out.tf.isin(FOC))]
f.to_csv(f"{RES}/seqreg_remap_focus_tfs_collagen.csv", index=False)
pd.set_option("display.width",250)
print(f[["region","window","tf","n_peaks","n_biotypes","n_peaks_fibroblast","n_peaks_breast","fibroblast_biotypes"]].to_string(index=False))

# --- appended note (see 35_seqreg_nfkb_context.py): SGBS is a Simpson-Golabi-Behmel
# preadipocyte line, i.e. mesenchymal lineage but not a fibroblast. It is NOT in the
# fibroblast whitelist above; 35_ treats it separately as "mesenchymal_lineage".
