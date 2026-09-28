#!/usr/bin/env python3
"""
PART B. Pan-cancer generalisation of the miR-29 -> collagen and TF -> collagen axes.
Independent implementation, written fresh 2026-09-09.
Data: UCSC Xena PANCAN EB++Adjust gene expression + pancanMiRs miRNA + phenotype.
"""
import gzip, os, sys, json, math
import numpy as np
from scipy import stats

BASE  = "/path/to/revision"
ML    = "/path/to/home/Desktop/DD/R_GPR/ML"
OUT   = os.path.join(BASE, "results/v2")
GEXP  = os.path.join(ML, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
MEXP  = os.path.join(ML, "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz")
PHENO = os.path.join(ML, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz")

MIN_N = 30          # minimum paired primary tumours per cohort

# ---------------------------------------------------------------- gene sets --
cafA = [l.strip() for l in open(os.path.join(OUT,"cafA_gene_list_v3.txt")) if l.strip()]
imm  = [l.strip() for l in open(os.path.join(OUT,"immune_gene_list_v3.txt")) if l.strip()]
cafB = ["DCN","LUM","FAP","THY1"]
FOCUS = ["COL1A1","COL3A1","COL1A2","COL5A1","FN1","SPARC","LOX",
         "ETS1","NFKB1","SP1","RELA","MYC","TP53","EZH2",
         "EPCAM","KRT8","KRT18","CDH1","PTPRC"]
WANT = sorted(set(cafA) | set(cafB) | set(imm) | set(FOCUS))
print(f"CAF-A genes requested: {len(cafA)}; immune: {len(imm)}; total wanted rows: {len(WANT)}", flush=True)

# ---------------------------------------------------------------- phenotype --
ph = {}
with gzip.open(PHENO, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        ph[p[0]] = (p[2], p[3])          # sample_type, primary_disease
print(f"phenotype rows: {len(ph)}", flush=True)

# ------------------------------------------------------------ miRNA matrix ---
with gzip.open(MEXP, "rt") as fh:
    mhdr = fh.readline().rstrip("\n").split("\t")
    msamp = mhdr[1:]
    mnames, mrows = [], []
    for line in fh:
        p = line.rstrip("\n").split("\t")
        mnames.append(p[0])
        mrows.append(np.array([np.nan if v in ("","NA") else float(v) for v in p[1:]], dtype=np.float32))
Mmat = np.vstack(mrows)
midx = {n:i for i,n in enumerate(mnames)}
print(f"miRNA matrix: {Mmat.shape[0]} x {Mmat.shape[1]}", flush=True)
assert len(set(mnames)) == len(mnames), "duplicate miRNA rownames"
_dupm = sorted({s for s in msamp if msamp.count(s) > 1})
print(f"duplicated miRNA sample columns: {len(_dupm)} -> {_dupm}; keeping FIRST occurrence", flush=True)

# ------------------------------------------------- stream gene matrix rows ---
want = set(WANT)
grows, gnames = [], []
with gzip.open(GEXP, "rt") as fh:
    ghdr = fh.readline().rstrip("\n").split("\t")
    gsamp = ghdr[1:]
    ncol = len(gsamp)
    colsum = np.zeros(ncol, dtype=np.float64)
    colcnt = np.zeros(ncol, dtype=np.int64)
    nrow_total = 0
    for line in fh:
        p = line.rstrip("\n").split("\t")
        nrow_total += 1
        vals = np.array([np.nan if v in ("","NA") else float(v) for v in p[1:]], dtype=np.float32)
        ok = np.isfinite(vals)
        colsum[ok] += vals[ok]; colcnt[ok] += 1
        if p[0] in want:
            gnames.append(p[0]); grows.append(vals)
Gmat = np.vstack(grows)
gidx = {}
for i,n in enumerate(gnames):
    gidx.setdefault(n, []).append(i)
dupg = [n for n,v in gidx.items() if len(v)>1]
print(f"gene matrix: {nrow_total} rows total x {ncol} samples; extracted {Gmat.shape[0]} rows "
      f"({len(gidx)} unique symbols); duplicated symbols among extracted: {dupg}", flush=True)
gmean_all = colsum/np.maximum(colcnt,1)     # per-sample mean over ALL genes
gidx1 = {n: v[0] for n,v in gidx.items()}   # first occurrence
missing = sorted(want - set(gidx1))
print(f"wanted genes not found in matrix: {len(missing)} -> {missing[:20]}", flush=True)

gpos = {}
for i,s in enumerate(gsamp): gpos.setdefault(s, i)     # keep FIRST occurrence
mpos = {}
for i,s in enumerate(msamp): mpos.setdefault(s, i)
print(f"unique gene-matrix samples {len(gpos)} of {len(gsamp)}; "
      f"unique miRNA-matrix samples {len(mpos)} of {len(msamp)}", flush=True)

# ------------------------------------------------------------ cohort split ---
common = [s for s in gsamp if s in mpos and s in ph and ph[s][0] == "Primary Tumor"]
print(f"primary tumours with BOTH gene and miRNA data: {len(common)}", flush=True)
coh = {}
for s in common:
    coh.setdefault(ph[s][1], []).append(s)
coh = {k:v for k,v in sorted(coh.items(), key=lambda kv: -len(kv[1]))}
for k,v in coh.items():
    print(f"  {k:45s} {len(v)}", flush=True)

# ------------------------------------------------------------- statistics ----
def spearman(x, y):
    k = np.isfinite(x) & np.isfinite(y)
    n = int(k.sum())
    if n < 12: return (np.nan, np.nan, n)
    if np.std(x[k]) == 0 or np.std(y[k]) == 0: return (np.nan, np.nan, n)
    r, p = stats.spearmanr(x[k], y[k])
    return (float(r), float(p), n)

def partial_spearman(x, y, z):
    k = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(k.sum())
    if n < 15: return (np.nan, np.nan, n)
    rx = stats.rankdata(x[k]); ry = stats.rankdata(y[k]); rz = stats.rankdata(z[k])
    if np.std(rx)==0 or np.std(ry)==0 or np.std(rz)==0: return (np.nan, np.nan, n)
    rxy = np.corrcoef(rx,ry)[0,1]; rxz = np.corrcoef(rx,rz)[0,1]; ryz = np.corrcoef(ry,rz)[0,1]
    den = math.sqrt(max((1-rxz**2)*(1-ryz**2), 1e-12))
    r = (rxy - rxz*ryz)/den
    r = max(min(r, 0.999999), -0.999999)
    t = r*math.sqrt((n-3)/(1-r**2))
    p = 2*stats.t.sf(abs(t), df=n-3)
    return (float(r), float(p), n)

MIRS = {"hsa-miR-29a":["hsa-miR-29a-3p","hsa-miR-29a-5p"],
        "hsa-miR-29b":["hsa-miR-29b-1-5p","hsa-miR-29b-2-5p","hsa-miR-29b-3p"],
        "hsa-miR-29c":["hsa-miR-29c-3p","hsa-miR-29c-5p"],
        "hsa-miR-29a-3p":["hsa-miR-29a-3p"],
        "hsa-miR-29b-3p":["hsa-miR-29b-3p"],
        "hsa-miR-29c-3p":["hsa-miR-29c-3p"],
        "hsa-miR-101":["hsa-miR-101-3p","hsa-miR-101-5p"],
        "hsa-let-7b":["hsa-let-7b-3p","hsa-let-7b-5p"]}
for lab, arms in MIRS.items():
    have = [a for a in arms if a in midx]
    print(f"  miRNA {lab}: arms present {have}", flush=True)
    MIRS[lab] = have

TFS  = ["ETS1","NFKB1","SP1","RELA","MYC","TP53"]
COLS = ["COL1A1","COL3A1"]

rows, csum, eco = [], [], []
for cohort, ss in coh.items():
    n = len(ss)
    if n < MIN_N:
        print(f"skipping {cohort}: n={n} < {MIN_N}", flush=True); continue
    gi = np.array([gpos[s] for s in ss]); mi = np.array([mpos[s] for s in ss])
    Gc = Gmat[:, gi]; Mc = Mmat[:, mi]
    gall = gmean_all[gi]

    def gvec(sym):
        i = gidx1.get(sym)
        return None if i is None else Gc[i].astype(np.float64)
    def mvec(lab):
        arms = MIRS[lab]
        if not arms: return None
        A = np.vstack([Mc[midx[a]] for a in arms]).astype(np.float64)
        return np.nanmean(A, axis=0)

    # ---- CAF scores: z within cohort, mean over genes
    def score(genes):
        rs = []
        for g in genes:
            v = gvec(g)
            if v is None: continue
            sd = np.nanstd(v)
            if not np.isfinite(sd) or sd == 0: continue
            rs.append((v - np.nanmean(v))/sd)
        return (np.mean(np.vstack(rs), axis=0), len(rs)) if rs else (None, 0)
    CAFA, nA = score(cafA); CAFB, nB = score(cafB); IMM, nI = score(imm)

    # ---- absolute (cross-cohort comparable) stromal index
    absA = []
    for g in cafA:
        v = gvec(g)
        if v is not None: absA.append(v)
    absA = np.vstack(absA)
    strom_abs = np.nanmean(absA, axis=0) - gall      # stromal mean minus all-gene mean
    col1 = gvec("COL1A1"); col3 = gvec("COL3A1")
    epi = [gvec(g) for g in ["EPCAM","KRT8","KRT18","CDH1"]]
    epi = np.nanmean(np.vstack([e for e in epi if e is not None]), axis=0) - gall

    csum.append(dict(cohort=cohort, n=n, cafA_genes_used=nA, cafB_genes_used=nB,
        immune_genes_used=nI,
        strom_abs_median=float(np.median(strom_abs)), strom_abs_mean=float(np.mean(strom_abs)),
        strom_abs_sd=float(np.std(strom_abs, ddof=1)),
        strom_abs_iqr=float(np.percentile(strom_abs,75)-np.percentile(strom_abs,25)),
        epi_abs_median=float(np.median(epi)),
        COL1A1_median=float(np.nanmedian(col1)), COL1A1_sd=float(np.nanstd(col1, ddof=1)),
        COL3A1_median=float(np.nanmedian(col3)), COL3A1_sd=float(np.nanstd(col3, ddof=1)),
        allgene_mean=float(np.mean(gall)),
        rho_CAFA_COL1A1=spearman(CAFA, col1)[0], rho_CAFA_COL3A1=spearman(CAFA, col3)[0],
        rho_COL1A1_COL3A1=spearman(col1, col3)[0]))

    pairs = []
    for m in ["hsa-miR-29a","hsa-miR-29b","hsa-miR-29c",
              "hsa-miR-29a-3p","hsa-miR-29b-3p","hsa-miR-29c-3p","hsa-miR-101","hsa-let-7b"]:
        for c in COLS: pairs.append((m,"miRNA",c))
    for t in TFS:
        for c in COLS: pairs.append((t,"gene",c))
    pairs += [("COL1A1","gene","COL3A1"), ("CAF_A","score","COL1A1"), ("CAF_A","score","COL3A1"),
              ("hsa-miR-101","miRNA","EZH2")]

    for src, sty, tgt in pairs:
        y = gvec(tgt)
        if y is None: continue
        x = mvec(src) if sty=="miRNA" else (CAFA if sty=="score" else gvec(src))
        if x is None: continue
        r,p,nn = spearman(x,y)
        pa = partial_spearman(x,y,CAFA) if sty!="score" else (np.nan,np.nan,nn)
        pb = partial_spearman(x,y,CAFB) if sty!="score" else (np.nan,np.nan,nn)
        rows.append(dict(cohort=cohort, n=nn, source=src, source_type=sty, target=tgt,
            axis=f"{src} -> {tgt}", rho=r, p=p,
            rho_partial_CAFA=pa[0], p_partial_CAFA=pa[1],
            rho_partial_CAFB=pb[0], p_partial_CAFB=pb[1]))
    print(f"done {cohort} (n={n})", flush=True)

import csv as _csv
def wr(path, dicts):
    if not dicts: return
    ks = list(dicts[0].keys())
    with open(path,"w",newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=ks); w.writeheader()
        for d in dicts: w.writerow(d)
    print("wrote", path, len(dicts), "rows", flush=True)

# BH across all reported primary correlations
ps = np.array([r["p"] if r["p"] is not None and np.isfinite(r["p"]) else 1.0 for r in rows])
order = np.argsort(ps); ranked = ps[order]; m = len(ps)
q = ranked*m/np.arange(1,m+1); q = np.minimum.accumulate(q[::-1])[::-1]
qq = np.empty(m); qq[order] = np.minimum(q,1.0)
for i,r in enumerate(rows): r["q_BH"] = float(qq[i])

wr(os.path.join(OUT,"pancancer_axis.csv"), rows)
wr(os.path.join(OUT,"pancancer_cohort_summary_v3.csv"), csum)

# ------------------------------------------------------- ecological test -----
import collections
byax = collections.defaultdict(dict)
for r in rows: byax[r["axis"]][r["cohort"]] = r
csum_by = {c["cohort"]: c for c in csum}
cohorts = [c["cohort"] for c in csum]
eco_rows = []
for ax, d in byax.items():
    cs = [c for c in cohorts if c in d and np.isfinite(d[c]["rho"])]
    if len(cs) < 8: continue
    rho = np.array([d[c]["rho"] for c in cs])
    nn  = np.array([d[c]["n"] for c in cs], dtype=float)
    z   = np.arctanh(np.clip(rho,-0.999999,0.999999))
    for stat in ["strom_abs_median","strom_abs_sd","strom_abs_iqr","epi_abs_median",
                 "COL1A1_sd","rho_CAFA_COL1A1","n"]:
        v = np.array([csum_by[c][stat] for c in cs], dtype=float)
        r1,p1 = stats.spearmanr(v, rho)
        r2,p2 = stats.pearsonr(v, z)
        eco_rows.append(dict(axis=ax, cohort_statistic=stat, k_cohorts=len(cs),
            spearman_rho=float(r1), spearman_p=float(p1),
            pearson_on_fisherz=float(r2), pearson_p=float(p2)))
    # heterogeneity across cohorts
    w = nn-3; zbar = float(np.sum(w*z)/np.sum(w))
    Q = float(np.sum(w*(z-zbar)**2)); df = len(cs)-1
    eco_rows.append(dict(axis=ax, cohort_statistic="HETEROGENEITY_Q", k_cohorts=len(cs),
        spearman_rho=Q, spearman_p=float(stats.chi2.sf(Q,df)),
        pearson_on_fisherz=float(np.tanh(zbar)), pearson_p=float(max(0,(Q-df)/Q)*100)))
wr(os.path.join(OUT,"pancancer_ecological_v3.csv"), eco_rows)

json.dump(dict(min_n=MIN_N, n_cohorts=len(csum), n_paired_samples=len(common),
               gene_rows_total=nrow_total, gene_rows_extracted=int(Gmat.shape[0]),
               mirna_rows=int(Mmat.shape[0]), missing_genes=missing),
          open(os.path.join(OUT,"pancancer_meta_v3.json"),"w"), indent=1)
print("ALL DONE", flush=True)
