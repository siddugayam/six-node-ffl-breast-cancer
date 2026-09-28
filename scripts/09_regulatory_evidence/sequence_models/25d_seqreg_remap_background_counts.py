#!/usr/bin/env python3
"""25d -- as 25c but storing the NUMBER of ReMap peaks per TF per background promoter, so
that a TF's peak count at COL1A1/COL3A1 can be given as a percentile of its peak count at
300 random promoters (presence/absence alone is uninformative: most TFs have a peak at
~80-90% of promoters)."""
import subprocess, tempfile, time, os
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
BB="https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMap2022.bb"; BBTOOL="/path/to/home/bin/bigBedToBed"
COLS=["chrom","start","end","name","score","strand","tstart","tend","rgb","tf","biotype"]
NBG=int(os.environ.get("NBG_REMAP2","300"))
rs=pd.read_csv(f"{CACHE}/genome/ncbiRefSeqSelect.txt.gz", sep="\t", header=None,
    names=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount",
           "exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"])
rs=rs[rs.name.str.startswith("NM_")]
MAIN=[f"chr{i}" for i in list(range(1,23))+["X"]]
rs=rs[rs.chrom.isin(MAIN)].drop_duplicates("name2")
rs["tss"]=np.where(rs.strand=="-",rs.txEnd,rs.txStart)
bg=rs.sample(n=NBG, random_state=424242)
def one(a):
    gene,chrom,tss,strand=a
    s=tss-1000 if strand=='+' else tss-500
    e=tss+500  if strand=='+' else tss+1000
    tmp=tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    try:
        subprocess.run([BBTOOL,f"-chrom={chrom}",f"-start={max(0,s)}",f"-end={e}",BB,tmp],
                       check=True, capture_output=True, timeout=300)
        d=pd.read_csv(tmp, sep="\t", header=None, names=COLS)
        g=d.groupby("tf").agg(n_peaks=("tf","size"), n_biotypes=("biotype","nunique")).reset_index()
        g["gene"]=gene
        return g
    except Exception:
        return None
    finally:
        try: os.unlink(tmp)
        except OSError: pass
jobs=[(r.name2,r.chrom,int(r.tss),r.strand) for _,r in bg.iterrows()]
res=[];fails=0;t0=time.time()
with ThreadPoolExecutor(max_workers=6) as ex:
    for i,g in enumerate(ex.map(one,jobs)):
        if g is None: fails+=1
        else: res.append(g)
        if (i+1)%50==0: print(f"  {i+1}/{len(jobs)} fails={fails} {time.time()-t0:.0f}s", flush=True)
B=pd.concat(res, ignore_index=True)
B.to_csv(f"{CACHE}/remap_background_promoter_tf_counts.csv.gz", index=False, compression="gzip")
ngene=B.gene.nunique(); print("background promoters:",ngene,"failures:",fails)

cur=pd.read_csv(f"{RES}/seqreg_remap_tf_by_region_curated.csv")
cur=cur[cur.window=="prom_1kb"].copy()
piv=B.pivot_table(index="gene", columns="tf", values="n_peaks", fill_value=0)
rows=[]
for _,r in cur.iterrows():
    v=piv[r.tf].values if r.tf in piv.columns else np.zeros(ngene)
    rows.append(dict(region=r.region, tf=r.tf, n_peaks=r.n_peaks, n_biotypes=r.n_biotypes,
                     n_peaks_fibroblast=r.n_peaks_fibroblast, n_peaks_breast=r.n_peaks_breast,
                     bg_mean_peaks=round(float(v.mean()),3), bg_median_peaks=float(np.median(v)),
                     bg_max_peaks=float(v.max()),
                     pct_bg_bound=round(float((v>0).mean()*100),2),
                     percentile_vs_bg=round(float((v<r.n_peaks).mean()*100),2),
                     emp_p_ge=round(float((v>=r.n_peaks).mean()),4), n_bg=ngene))
out=pd.DataFrame(rows)
out.to_csv(f"{RES}/seqreg_remap_promoter_tfs_calibrated.csv", index=False)
pd.set_option("display.width",240)
FOC=["NFKB1","NFKB2","RELA","RELB","SP1","SP2","SP3","ETS1","ETS2","MED1","EP300","CTCF","ESR1","TWIST1","STAT3","SMAD3"]
for reg in ["COL1A1","COL3A1","MIR29B2CHG","LINC_PINT_prox","MIR130AHG"]:
    print(f"\n===== {reg}")
    print(out[(out.region==reg)&(out.tf.isin(FOC))][["tf","n_peaks","n_peaks_fibroblast","bg_mean_peaks",
          "bg_median_peaks","pct_bg_bound","percentile_vs_bg","emp_p_ge"]].to_string(index=False))
