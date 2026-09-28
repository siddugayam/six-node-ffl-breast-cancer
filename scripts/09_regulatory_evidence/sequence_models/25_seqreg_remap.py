#!/usr/bin/env python3
"""25_seqreg_remap.py -- ReMap 2022 (Hammal et al., NAR 2022) experimentally mapped TF peaks
at the COL1A1 / COL3A1 / miR-29 / miR-130a loci (hg38), with cell-type (biotype) resolution
and a random-promoter background so that "TF X binds here" is calibrated.

Peaks pulled by range query from the UCSC-hosted ReMap 2022 bigBed
(/gbdb/hg38/reMap/reMap2022.bb) with bigBedToBed. Columns 10/11 are TF and biotype.
"""
import os, subprocess, tempfile, time
import numpy as np, pandas as pd

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
BB="https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMap2022.bb"
BBTOOL="/path/to/home/bin/bigBedToBed"
COLS=["chrom","start","end","name","score","strand","tstart","tend","rgb","tf","biotype"]

peaks=pd.read_csv(f"{CACHE}/remap2022_hg38_loci.bed", sep="\t", header=None, names=COLS)
print("peaks loaded:", len(peaks), "TFs:", peaks.tf.nunique(), "biotypes:", peaks.biotype.nunique())

FIB_PAT=r"(?i)(fibroblast|BJ|IMR-?90|MRC-?5|WI-?38|HFF|NHLF|NHDF|HDF|CCD-|myofibroblast|stellate|mesenchymal|MSC|HS27|HS5|Detroit|GM23248|AG0|hTERT[- ]?HM)"
peaks["is_fibro_like"]=peaks.biotype.str.contains(FIB_PAT, regex=True, na=False)
print("fibroblast-like biotypes:", sorted(peaks.loc[peaks.is_fibro_like,'biotype'].unique())[:60])

reg=pd.read_csv(f"{RES}/seqreg_regions.csv")
def window(r, up, dn):
    a=int(r.anchor)
    return (a-up, a+dn) if r.strand=='+' else (a-dn, a+up)

rows=[]; percell=[]
for _,r in reg.iterrows():
    for label,(up,dn) in {"prom_1kb":(1000,500),"prom_5kb":(5000,2000),"locus_100kb":(100000,100000)}.items():
        s,e=window(r,up,dn)
        sub=peaks[(peaks.chrom==r.chrom)&(peaks.start<e)&(peaks.end>s)]
        for tf,g in sub.groupby("tf"):
            rows.append(dict(region=r.region, window=label, tf=tf, n_peaks=len(g),
                             n_biotypes=g.biotype.nunique(),
                             n_fibro_peaks=int(g.is_fibro_like.sum()),
                             fibro_biotypes=";".join(sorted(set(g.loc[g.is_fibro_like,'biotype'])))[:400],
                             biotypes=";".join(sorted(set(g.biotype))[:25])))
        if label=="prom_1kb":
            for (tf,bt),g in sub.groupby(["tf","biotype"]):
                percell.append(dict(region=r.region, tf=tf, biotype=bt, n_peaks=len(g),
                                    is_fibro_like=bool(g.is_fibro_like.iloc[0])))
        print(r.region,label,"TFs:",sub.tf.nunique(), flush=True)
pd.DataFrame(rows).to_csv(f"{RES}/seqreg_remap_tf_by_region.csv", index=False)
pd.DataFrame(percell).to_csv(f"{RES}/seqreg_remap_promoter_tf_celltype.csv", index=False)

# ------- random-promoter background: how often is each TF seen at a promoter? ------
import pyfaidx
rs=pd.read_csv(f"{CACHE}/genome/ncbiRefSeqSelect.txt.gz", sep="\t", header=None,
               names=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd",
                      "exonCount","exonStarts","exonEnds","score","name2","cdsStartStat",
                      "cdsEndStat","exonFrames"])
rs=rs[rs.name.str.startswith("NM_")]
MAIN=[f"chr{i}" for i in list(range(1,23))+["X"]]
rs=rs[rs.chrom.isin(MAIN)].drop_duplicates("name2")
rs["tss"]=np.where(rs.strand=="-", rs.txEnd, rs.txStart)
NBG=int(os.environ.get("NBG_REMAP","300"))
bg=rs.sample(n=NBG, random_state=20260909)

recs=[]; t0=time.time()
tmp=tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
for i,(_,r) in enumerate(bg.iterrows()):
    s=int(r.tss)-1000 if r.strand=='+' else int(r.tss)-500
    e=int(r.tss)+500  if r.strand=='+' else int(r.tss)+1000
    try:
        subprocess.run([BBTOOL,f"-chrom={r.chrom}",f"-start={max(0,s)}",f"-end={e}",BB,tmp],
                       check=True, capture_output=True, timeout=300)
        d=pd.read_csv(tmp, sep="\t", header=None, names=COLS)
    except Exception as ex:
        print("bg fail", r.name2, ex); continue
    for tf in d.tf.unique():
        recs.append((r.name2, tf))
    if (i+1)%25==0: print(f"  bg {i+1}/{NBG} {time.time()-t0:.0f}s", flush=True)
bgdf=pd.DataFrame(recs, columns=["gene","tf"])
bgdf.to_csv(f"{CACHE}/remap_background_promoter_tfs.csv", index=False)
ngenes=bgdf.gene.nunique()
freq=bgdf.groupby("tf").gene.nunique()/ngenes
freq.name="frac_bg_promoters_bound"
freq.to_csv(f"{RES}/seqreg_remap_background_tf_frequency.csv")
print("background promoters used:", ngenes)

tab=pd.DataFrame(rows)
tab=tab[tab.window=="prom_1kb"].merge(freq.reset_index(), on="tf", how="left")
tab["frac_bg_promoters_bound"]=tab.frac_bg_promoters_bound.fillna(0)
tab.to_csv(f"{RES}/seqreg_remap_promoter_tfs_calibrated.csv", index=False)
print("done")
