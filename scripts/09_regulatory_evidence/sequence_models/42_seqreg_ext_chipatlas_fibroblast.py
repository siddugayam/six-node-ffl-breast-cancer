#!/usr/bin/env python3
"""42_seqreg_ext_chipatlas_fibroblast.py

EXTENSION of 24_seqreg_chipatlas.py.

The first pass answered "is TF X bound at COL1A1?" from ChIP-Atlas's PRE-COMPUTED
target-gene tables, and 203 of 3,702 of those requests failed (HTTP 404), including
MKL1/MRTFA. This script instead downloads the ACTUAL PEAK FILES of every
fibroblast/mesenchymal TF ChIP-seq experiment in ChIP-Atlas hg38 and intersects them with
the loci directly, so the fibroblast result is a direct peak-level observation rather than
an inherited annotation.

Threshold: bed05 = MACS2 q < 1e-05, ChIP-Atlas's standard "significant peak" set.

Cell-type selection is by explicit regex on the ChIP-Atlas celltype field, and the matched
celltype list is written out so the reader can audit exactly what was called a fibroblast.
"""
import os, re, io, sys, time
import concurrent.futures as cf
import numpy as np, pandas as pd, requests

ROOT = "/path/to/revision"
RES, CACHE = f"{ROOT}/results/v3", f"{ROOT}/cache/seqreg"
BEDDIR = f"{CACHE}/chipatlas_beds"
os.makedirs(BEDDIR, exist_ok=True)

FIB_PAT = re.compile(
    r"fibroblast|IMR-?90|WI-?38|BJ\b|MRC-?5|HFF|NHDF|HDF|MSC|mesenchymal|"
    r"myofibroblast|stellate|HS27|CCD-?1|TIG-?|AG0|GM23248|LX-2", re.I)

cols = ["srx", "genome", "ag_class", "antigen", "ct_class", "celltype", "desc"]
d = pd.read_csv(f"{CACHE}/chipatlas_hg38.tsv", sep="\t", names=cols, low_memory=False)
tf = d[d.ag_class == "TFs and others"].copy()
fb = tf[tf.celltype.str.contains(FIB_PAT, na=False)].copy()
print(f"hg38 TF experiments {len(tf)}; fibroblast/mesenchymal {len(fb)} "
      f"({fb.celltype.nunique()} celltypes, {fb.antigen.nunique()} antigens)", flush=True)
fb[["srx", "antigen", "ct_class", "celltype", "desc"]].to_csv(
    f"{RES}/seqreg_ext_chipatlas_fibroblast_experiment_list.csv", index=False)

reg = pd.read_csv(f"{RES}/seqreg_regions.csv")
WINDOWS = {"promoter_1kb_500": (-1000, 500), "pm5kb": (-5000, 5000), "pm100kb": (-100000, 100000)}
regions = []
for _, r in reg.iterrows():
    a = int(r.anchor); strand = r.strand
    for wname, (lo, hi) in WINDOWS.items():
        # convert transcription-oriented offsets to plus-strand genomic coordinates
        if strand == "+":
            gs, ge = a + lo, a + hi
        else:
            gs, ge = a - hi, a - lo
        regions.append(dict(region=r.region, window=wname, chrom=r.chrom,
                            start=gs, end=ge, anchor=a, strand=strand))
regions = pd.DataFrame(regions)
CHROMS = set(regions.chrom)
regions.to_csv(f"{RES}/seqreg_ext_chipatlas_windows.csv", index=False)

URL = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bed05/{srx}.05.bed"
sess = requests.Session()

def fetch(srx):
    """Download a peak file once, keeping only the chromosomes we care about."""
    p = f"{BEDDIR}/{srx}.bed"
    if os.path.exists(p):
        return p, None
    for attempt in range(3):
        try:
            r = sess.get(URL.format(srx=srx), timeout=180)
            if r.status_code != 200:
                if attempt == 2:
                    return None, f"HTTP {r.status_code}"
                time.sleep(2); continue
            keep = []
            for line in r.text.splitlines():
                if not line or line.startswith("track"):
                    continue
                f0 = line.split("\t", 1)[0]
                if f0 in CHROMS:
                    keep.append(line)
            with open(p, "w") as fh:
                fh.write("\n".join(keep))
            return p, None
        except Exception as e:
            if attempt == 2:
                return None, str(e)[:100]
            time.sleep(2)
    return None, "unreachable"

srxs = fb.srx.tolist()
paths, errors = {}, []
t0 = time.time()
with cf.ThreadPoolExecutor(max_workers=8) as ex:
    for n, (srx, (p, err)) in enumerate(zip(srxs, ex.map(fetch, srxs))):
        if p:
            paths[srx] = p
        else:
            errors.append(dict(srx=srx, error=err))
        if n % 50 == 0:
            print(f"  downloaded {n}/{len(srxs)}  {time.time()-t0:.0f}s  errors={len(errors)}", flush=True)
print(f"downloaded {len(paths)}/{len(srxs)}; {len(errors)} errors  {time.time()-t0:.0f}s", flush=True)
pd.DataFrame(errors).to_csv(f"{RES}/seqreg_ext_chipatlas_download_errors.csv", index=False)

meta = fb.set_index("srx")
rows = []
for n, (srx, p) in enumerate(paths.items()):
    try:
        b = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2, 4],
                        names=["chrom", "start", "end", "score"], low_memory=False)
    except Exception:
        continue
    if not len(b):
        continue
    m = meta.loc[srx]
    for _, w in regions.iterrows():
        sel = b[(b.chrom == w.chrom) & (b.end > w.start) & (b.start < w.end)]
        if not len(sel):
            continue
        centres = ((sel.start + sel.end) / 2).to_numpy()
        dist = centres - w.anchor
        if w.strand == "-":
            dist = -dist
        j = int(np.argmin(np.abs(dist)))
        rows.append(dict(srx=srx, antigen=m.antigen, celltype=m.celltype,
                         ct_class=m.ct_class, region=w.region, window=w.window,
                         n_peaks=len(sel), max_score=float(sel.score.max()),
                         nearest_dist_to_anchor=int(dist[j])))
    if n % 100 == 0:
        print(f"  intersected {n}/{len(paths)}", flush=True)

hits = pd.DataFrame(rows)
hits.to_csv(f"{RES}/seqreg_ext_chipatlas_fibroblast_peak_hits.csv", index=False)
print(f"\nhit rows: {len(hits)}", flush=True)

# per TF x region x window: how many of that TF's fibroblast experiments have a peak
n_exp = fb.groupby("antigen").srx.nunique().rename("n_fibroblast_experiments")
summ = []
for (ag, rg, w), g in hits.groupby(["antigen", "region", "window"]):
    summ.append(dict(antigen=ag, region=rg, window=w,
                     n_experiments_with_peak=g.srx.nunique(),
                     n_fibroblast_experiments=int(n_exp.get(ag, 0)),
                     celltypes_with_peak=";".join(sorted(g.celltype.unique())),
                     max_score=float(g.max_score.max()),
                     nearest_dist_to_anchor=int(g.nearest_dist_to_anchor.abs().min())))
summ = pd.DataFrame(summ).sort_values(["region", "window", "n_experiments_with_peak"],
                                      ascending=[True, True, False])
summ.to_csv(f"{RES}/seqreg_ext_chipatlas_fibroblast_tf_summary.csv", index=False)

# explicit focus-TF record, including the zero rows that matter most
FOCUS = ["NFKB1", "RELA", "RELB", "NFKB2", "REL", "SP1", "SP2", "SP3",
         "ETS1", "ETS2", "MKL1", "MRTFA", "SRF", "JUN", "FOS", "SMAD3", "CTCF"]
frows = []
for ag in FOCUS:
    ne = int(n_exp.get(ag, 0))
    for rg in ["COL1A1", "COL3A1", "MIR29A", "MIR29B1", "MIR29C", "MIR29B2",
               "MIR29B2CHG", "LINC_PINT_prox"]:
        for w in WINDOWS:
            g = hits[(hits.antigen == ag) & (hits.region == rg) & (hits.window == w)]
            frows.append(dict(antigen=ag, region=rg, window=w,
                              n_fibroblast_experiments=ne,
                              n_experiments_with_peak=int(g.srx.nunique()),
                              celltypes_with_peak=";".join(sorted(g.celltype.unique())) if len(g) else "",
                              nearest_dist_to_anchor=(int(g.nearest_dist_to_anchor.abs().min())
                                                      if len(g) else np.nan)))
pd.DataFrame(frows).to_csv(f"{RES}/seqreg_ext_chipatlas_focus_tfs.csv", index=False)

print("\n--- COL1A1 / COL3A1 promoter (TSS-1000..+500), fibroblast experiments ---", flush=True)
for rg in ["COL1A1", "COL3A1"]:
    s = summ[(summ.region == rg) & (summ.window == "promoter_1kb_500")]
    print(f"\n{rg}: {len(s)} TFs with a fibroblast peak in the promoter")
    print(s.head(25)[["antigen", "n_experiments_with_peak", "n_fibroblast_experiments",
                      "celltypes_with_peak", "nearest_dist_to_anchor"]].to_string(index=False))
print("\nDONE", flush=True)
