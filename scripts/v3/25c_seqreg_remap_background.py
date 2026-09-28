#!/usr/bin/env python3
"""25c -- ReMap 2022 random-promoter background: for 300 random RefSeq-Select promoters,
how often is each TF observed with a peak in the same TSS-1000..+500 window? This turns
"TF X has a ReMap peak at the COL1A1 promoter" into a calibrated statement.
"""
import subprocess, tempfile, time, os
import numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
BB="https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMap2022.bb"; BBTOOL="/path/to/home/bin/bigBedToBed"
COLS=["chrom","start","end","name","score","strand","tstart","tend","rgb","tf","biotype"]
NBG=int(os.environ.get("NBG_REMAP","300"))

rs=pd.read_csv(f"{CACHE}/genome/ncbiRefSeqSelect.txt.gz", sep="\t", header=None,
    names=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount",
           "exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"])
rs=rs[rs.name.str.startswith("NM_")]
MAIN=[f"chr{i}" for i in list(range(1,23))+["X"]]
rs=rs[rs.chrom.isin(MAIN)].drop_duplicates("name2")
rs["tss"]=np.where(rs.strand=="-",rs.txEnd,rs.txStart)
bg=rs.sample(n=NBG, random_state=20260909)

from concurrent.futures import ThreadPoolExecutor
def one(args):
    gene,chrom,tss,strand=args
    s=tss-1000 if strand=='+' else tss-500
    e=tss+500  if strand=='+' else tss+1000
    tmp=tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    try:
        subprocess.run([BBTOOL,f"-chrom={chrom}",f"-start={max(0,s)}",f"-end={e}",BB,tmp],
                       check=True, capture_output=True, timeout=300)
        d=pd.read_csv(tmp, sep="\t", header=None, names=COLS)
        out=(gene, sorted(d.tf.unique()))
    except Exception as ex:
        out=(gene, None)
    finally:
        try: os.unlink(tmp)
        except OSError: pass
    return out

jobs=[(r.name2,r.chrom,int(r.tss),r.strand) for _,r in bg.iterrows()]
recs=[]; fails=0; t0=time.time()
with ThreadPoolExecutor(max_workers=6) as ex:
    for i,(gene,tfs) in enumerate(ex.map(one, jobs)):
        if tfs is None: fails+=1; continue
        for t in tfs: recs.append((gene,t))
        if (i+1)%50==0: print(f"  {i+1}/{len(jobs)} fails={fails} {time.time()-t0:.0f}s", flush=True)
bgdf=pd.DataFrame(recs, columns=["gene","tf"])
bgdf.to_csv(f"{CACHE}/remap_background_promoter_tfs.csv", index=False)
n=bgdf.gene.nunique()
freq=(bgdf.groupby("tf").gene.nunique()/n).rename("frac_bg_promoters_bound").reset_index()
freq["n_bg_promoters"]=n
freq.to_csv(f"{RES}/seqreg_remap_background_tf_frequency.csv", index=False)
print("background promoters used:", n, "failures:", fails)

cur=pd.read_csv(f"{RES}/seqreg_remap_tf_by_region_curated.csv")
cur=cur[cur.window=="prom_1kb"].merge(freq, on="tf", how="left")
cur["frac_bg_promoters_bound"]=cur.frac_bg_promoters_bound.fillna(0.0)
cur["n_bg_promoters"]=cur.n_bg_promoters.fillna(n)
cur.to_csv(f"{RES}/seqreg_remap_promoter_tfs_calibrated.csv", index=False)
pd.set_option("display.width",240)
FOC=["NFKB1","NFKB2","RELA","RELB","REL","SP1","SP2","SP3","ETS1","ETS2","MED1","EP300","CTCF","ESR1","TWIST1","STAT3","SMAD3"]
print(cur[(cur.region.isin(["COL1A1","COL3A1","MIR29B2CHG","MIR130AHG","LINC_PINT_prox"]))&(cur.tf.isin(FOC))]
      [["region","tf","n_peaks","n_biotypes","n_peaks_fibroblast","frac_bg_promoters_bound"]].to_string(index=False))
