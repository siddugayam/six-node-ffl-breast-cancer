#!/usr/bin/env python3
"""
TASK B: pan-cancer generalisation of the miR-29 -> collagen axis and of the
        TF -> collagen (stromal-confounded) axis, across all TCGA cohorts.

Also: an ECOLOGICAL test -- does the strength of the TF->collagen correlation
      track cohort-level stromal (CAF) content ACROSS cohorts?

All numbers computed in this run. Iteration counts are printed explicitly.
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats

SCR  = "/path/to/scratch"
ML   = "/path/to/home/Desktop/DD/R_GPR/ML"
BASE = "/path/to/revision"
OUT  = os.path.join(BASE, "results/v2")
os.makedirs(OUT, exist_ok=True)
MIN_N = 30                       # minimum paired primary tumours per cohort

print("="*78); print("PART B : PAN-CANCER"); print("="*78)

# ------------------------------------------------------------------ inputs --
G = pd.read_csv(os.path.join(SCR, "pancan_genes_subset.tsv"), sep="\t", index_col=0)
M = pd.read_csv(os.path.join(SCR, "pancan_mir29.tsv"),        sep="\t", index_col=0)
print(f"gene subset matrix : {G.shape[0]} genes x {G.shape[1]} samples")
print(f"miR-29 matrix      : {M.shape[0]} arms  x {M.shape[1]} samples")
assert G.index.duplicated().sum() == 0, "duplicated gene rownames"
assert M.index.duplicated().sum() == 0, "duplicated miRNA rownames"
assert "COL1A1" in G.index and "COL3A1" in G.index and "ETS1" in G.index and "NFKB1" in G.index

PH = pd.read_csv(os.path.join(ML, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"),
                 sep="\t", index_col=0)
print(f"phenotype rows     : {PH.shape[0]}  cols: {list(PH.columns)}")

# primary solid tumours only (sample_type_id 01)
prim = PH[PH["sample_type_id"] == 1].copy()
print(f"primary tumour samples in phenotype: {len(prim)}")

sigA = [l.strip() for l in open(os.path.join(SCR, "sigA.txt")) if l.strip()]
sigB = [l.strip() for l in open(os.path.join(SCR, "sigB.txt")) if l.strip()]
assert len(sigA) == 123 and len(sigB) == 4
assert not any(g.startswith("COL") for g in sigA + sigB), "collagen leaked into CAF signature"
sigA_have = [g for g in sigA if g in G.index]
sigB_have = [g for g in sigB if g in G.index]
print(f"CAF-A genes present in pan-cancer matrix: {len(sigA_have)}/123")
print(f"CAF-B genes present: {len(sigB_have)}/4 -> {sigB_have}")

# --------------------------------------------------------------- alignment --
gs = set(G.columns); ms = set(M.columns); ps = set(prim.index)
paired = sorted(gs & ms & ps)
print(f"\nprimary tumours with BOTH gene and miRNA assays: n = {len(paired)}")
gene_only = sorted(gs & ps)
print(f"primary tumours with gene assay (for TF axes)  : n = {len(gene_only)}")

cohort = prim["_primary_disease"]

# ---------------------------------------------------------- stat helpers ---
def sr(x, y):
    """Spearman rho, p, n -- NaN tolerant (pairwise complete)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
        return np.nan, np.nan, int(ok.sum())
    r = stats.spearmanr(x[ok], y[ok])
    return float(r.statistic), float(r.pvalue), int(ok.sum())

def prho(x, y, z):
    """partial Spearman of x,y given z -- NaN tolerant."""
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 10: return np.nan
    rx, ry, rz = (stats.rankdata(v[ok]) for v in (x, y, z))
    rxy = np.corrcoef(rx, ry)[0, 1]; rxz = np.corrcoef(rx, rz)[0, 1]; ryz = np.corrcoef(ry, rz)[0, 1]
    d = np.sqrt((1 - rxz**2) * (1 - ryz**2))
    return float((rxy - rxz * ryz) / d) if d > 0 else np.nan

# ------------------------------------------------------------- CAF scoring --
# z-score each signature gene ACROSS ALL primary tumours pooled, so that the
# per-cohort MEAN of the score is interpretable as that cohort's stromal content.
Gp = G[gene_only]
nan_per_gene = Gp.isna().sum(axis=1)
nan_samples = int((Gp.isna().sum(axis=0) > 0).sum())
print(f"missing values: {int(Gp.isna().sum().sum())} cells; "
      f"{int((nan_per_gene>0).sum())} genes affected; {nan_samples} samples affected")
if (nan_per_gene > 0).any():
    print("  genes with missing values:", list(nan_per_gene[nan_per_gene > 0].index))

def pooled_score(sig):
    """mean of per-gene z-scores, z computed across ALL pooled primary tumours;
    NaN-tolerant so the 211 samples missing 9 genes still get a score."""
    X = Gp.loc[sig].astype(float).values
    mu = np.nanmean(X, axis=1)[:, None]
    sd = np.nanstd(X, axis=1, ddof=1)[:, None]
    keep = (sd[:, 0] > 0) & np.isfinite(sd[:, 0])
    Z = (X[keep] - mu[keep]) / sd[keep]
    return pd.Series(np.nanmean(Z, axis=0), index=Gp.columns), int(keep.sum())

cafA, nA = pooled_score(sigA_have)
cafB, nB = pooled_score(sigB_have)
assert np.isfinite(cafA).all() and np.isfinite(cafB).all(), "CAF score has non-finite values"
print(f"CAF-A score built from {nA} genes; CAF-B from {nB} genes (pooled z across "
      f"{len(gene_only)} primary tumours)")
print(f"CAF-A vs CAF-B Spearman across all primary tumours: {sr(cafA, cafB)[0]:.3f}")
print(f"CAF-A vs COL1A1 rho = {sr(cafA, Gp.loc['COL1A1'])[0]:.3f}; "
      f"CAF-A vs COL3A1 rho = {sr(cafA, Gp.loc['COL3A1'])[0]:.3f}")

# ------------------------------------------------------------ miRNA arms ----
# two mappings, exactly as in the BRCA/CPTAC work:
#   arm-averaged : mean over every mature arm carrying the label
#   dominant arm : the single most abundant arm (-3p for all three miR-29s)
arms = {
    "hsa-miR-29a": ["hsa-miR-29a-3p", "hsa-miR-29a-5p"],
    "hsa-miR-29b": ["hsa-miR-29b-1-5p", "hsa-miR-29b-2-5p", "hsa-miR-29b-3p"],
    "hsa-miR-29c": ["hsa-miR-29c-3p", "hsa-miR-29c-5p"],
}
print("\nmean abundance per arm (all samples):")
for a in M.index:
    print(f"   {a:20s} {M.loc[a].mean():.3f}")
dom = {lab: max(v, key=lambda a: M.loc[a].mean()) for lab, v in arms.items()}
print("dominant arm per label:", dom)

Mavg = pd.DataFrame({lab: M.loc[v].mean(axis=0) for lab, v in arms.items()}).T
Mdom = pd.DataFrame({lab: M.loc[a] for lab, a in dom.items()}).T

MIRS = ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"]
TFS  = ["ETS1", "NFKB1", "SP1", "RELA"]
COLS = ["COL1A1", "COL3A1"]

# ------------------------------------------------------------- per cohort ---
cohorts = sorted(cohort.loc[gene_only].dropna().unique())
rows, coh_rows = [], []
print(f"\ncandidate cohorts: {len(cohorts)}")

for ct in cohorts:
    s_gene = [s for s in gene_only if cohort.get(s) == ct]
    s_pair = [s for s in paired    if cohort.get(s) == ct]
    if len(s_gene) < MIN_N:
        print(f"  SKIP {ct}: only {len(s_gene)} gene-assay primary tumours")
        continue

    cA = cafA[s_gene]; cB = cafB[s_gene]
    c1 = Gp.loc["COL1A1", s_gene]; c3 = Gp.loc["COL3A1", s_gene]

    # rho of each collagen with the cohort's CAF score (needed for the
    # confounding decomposition rho_confound = rho(src,CAF) * rho(tgt,CAF))
    rc1 = sr(c1, cA)[0]; rc3 = sr(c3, cA)[0]

    # --- TF -> collagen (gene-assay samples) ---
    for tf in TFS:
        tv = Gp.loc[tf, s_gene]
        r_src_caf = sr(tv, cA)[0]
        for col, cv, r_tgt_caf in (("COL1A1", c1, rc1), ("COL3A1", c3, rc3)):
            r, p, n = sr(tv, cv)
            rows.append(dict(cohort=ct, axis=f"{tf} -> {col}", source=tf, target=col,
                             layer="mRNA", n=n, rho=r, p=p,
                             rho_partial_cafA=prho(tv, cv, cA),
                             rho_partial_cafB=prho(tv, cv, cB),
                             rho_source_cafA=r_src_caf, rho_target_cafA=r_tgt_caf,
                             rho_confound=r_src_caf * r_tgt_caf,
                             mir_mapping=""))
    # --- COL1A1 ~ COL3A1 ---
    r, p, n = sr(c1, c3)
    rows.append(dict(cohort=ct, axis="COL1A1 ~ COL3A1", source="COL1A1", target="COL3A1",
                     layer="mRNA", n=n, rho=r, p=p,
                     rho_partial_cafA=prho(c1, c3, cA), rho_partial_cafB=prho(c1, c3, cB),
                     rho_source_cafA=rc1, rho_target_cafA=rc3, rho_confound=rc1 * rc3,
                     mir_mapping=""))
    # --- miR-29 -> collagen (paired samples) ---
    n_pair = len(s_pair)
    if n_pair >= MIN_N:
        cAp = cafA[s_pair]; cBp = cafB[s_pair]
        pc1 = sr(Gp.loc["COL1A1", s_pair], cAp)[0]
        pc3 = sr(Gp.loc["COL3A1", s_pair], cAp)[0]
        for mp_name, MM in (("arm_averaged", Mavg), ("dominant_arm", Mdom)):
            for mir in MIRS:
                mv = MM.loc[mir, s_pair]
                r_src_caf = sr(mv, cAp)[0]
                for col, r_tgt_caf in (("COL1A1", pc1), ("COL3A1", pc3)):
                    cv = Gp.loc[col, s_pair]
                    r, p, n = sr(mv, cv)
                    rows.append(dict(cohort=ct, axis=f"{mir} -> {col}", source=mir, target=col,
                                     layer="mRNA", n=n, rho=r, p=p,
                                     rho_partial_cafA=prho(mv, cv, cAp),
                                     rho_partial_cafB=prho(mv, cv, cBp),
                                     rho_source_cafA=r_src_caf, rho_target_cafA=r_tgt_caf,
                                     rho_confound=r_src_caf * r_tgt_caf,
                                     mir_mapping=mp_name))
    # --- cohort-level summary (stromal content etc.) ---
    coh_rows.append(dict(
        cohort=ct, n_gene=len(s_gene), n_paired=n_pair,
        cafA_mean=float(cA.mean()), cafA_sd=float(cA.std(ddof=1)), cafA_median=float(cA.median()),
        cafB_mean=float(cB.mean()), cafB_sd=float(cB.std(ddof=1)),
        COL1A1_mean=float(c1.mean()), COL3A1_mean=float(c3.mean()),
        COL1A1_sd=float(c1.std(ddof=1)),
        EPCAM_mean=float(Gp.loc["EPCAM", s_gene].mean()) if "EPCAM" in Gp.index else np.nan,
        PTPRC_mean=float(Gp.loc["PTPRC", s_gene].mean()) if "PTPRC" in Gp.index else np.nan,
    ))

R = pd.DataFrame(rows)
C = pd.DataFrame(coh_rows)
print(f"\ncohorts analysed (gene assay, n>={MIN_N}): {len(C)}")
print(f"cohorts with paired miRNA (n>={MIN_N})    : {int((C.n_paired>=MIN_N).sum())}")

# ---------------------------------------------------------- headline table --
def show(axis, mapping=""):
    d = R[(R.axis == axis) & (R.mir_mapping == mapping)].sort_values("rho")
    print(f"\n--- {axis} ({mapping or 'mRNA'}) ---")
    for _, x in d.iterrows():
        star = "***" if x.p < 1e-4 else ("**" if x.p < 1e-2 else ("*" if x.p < .05 else "  "))
        print(f"  {x.cohort[:44]:44s} n={x.n:4d} rho={x.rho:+.3f} p={x.p:.2e} {star}  "
              f"partial(CAF-A)={x.rho_partial_cafA:+.3f}")

for mir in MIRS:
    for col in COLS:
        show(f"{mir} -> {col}", "arm_averaged")
for tf in ["ETS1", "NFKB1"]:
    show(f"{tf} -> COL1A1")

# ------------------------------------------- how consistent is miR-29 axis? --
print("\n" + "="*78)
print("CONSISTENCY OF THE miR-29 -> COLLAGEN AXIS ACROSS COHORTS")
print("="*78)
summ = []
for mp in ("arm_averaged", "dominant_arm"):
    for mir in MIRS:
        for col in COLS:
            d = R[(R.axis == f"{mir} -> {col}") & (R.mir_mapping == mp)]
            if not len(d): continue
            neg = int((d.rho < 0).sum()); k = len(d)
            signif_neg = int(((d.rho < 0) & (d.p < 0.05)).sum())
            signif_pos = int(((d.rho > 0) & (d.p < 0.05)).sum())
            # sign test that the axis is negative more often than chance
            bt = stats.binomtest(neg, k, 0.5, alternative="greater").pvalue
            # one-sample t on Fisher-z (are cohort correlations negative on average?)
            z = np.arctanh(d.rho.values)
            tt = stats.ttest_1samp(z, 0.0)
            summ.append(dict(mapping=mp, axis=f"{mir} -> {col}", k_cohorts=k,
                             n_negative=neg, n_signif_negative=signif_neg,
                             n_signif_positive=signif_pos,
                             median_rho=float(d.rho.median()),
                             min_rho=float(d.rho.min()), max_rho=float(d.rho.max()),
                             rho_BRCA=float(d[d.cohort.str.contains("breast")].rho.iloc[0])
                                       if (d.cohort.str.contains("breast")).any() else np.nan,
                             sign_test_p=float(bt), t_on_fisherz=float(tt.statistic),
                             p_fisherz=float(tt.pvalue)))
S = pd.DataFrame(summ)
print(S.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# ---------------------------------------------------- ECOLOGICAL CONFOUND ---
print("\n" + "="*78)
print("ECOLOGICAL TEST: does TF->collagen correlation strength track cohort stromal content?")
print("="*78)
eco = []
ECO_VARS = ("rho_confound",   # rho(source,CAF) * rho(target,CAF) -- the direct confound term
            "cafA_sd",        # within-cohort spread of stromal content
            "cafA_mean",      # cohort mean stromal level
            "cafB_mean", "COL1A1_sd")
for axis in [f"{tf} -> {col}" for tf in TFS for col in COLS] + \
            [f"{m} -> {c}" for m in MIRS for c in COLS] + ["COL1A1 ~ COL3A1"]:
    mp = "arm_averaged" if axis.startswith("hsa-") else ""
    d = R[(R.axis == axis) & (R.mir_mapping == mp)].merge(C, on="cohort")
    if len(d) < 8: continue
    for xvar in ECO_VARS:
        v = d[xvar].values.astype(float)
        y = d.rho.values.astype(float)
        ok = np.isfinite(v) & np.isfinite(y)
        rr = stats.spearmanr(v[ok], y[ok])
        pp = stats.pearsonr(v[ok], np.arctanh(y[ok]))
        eco.append(dict(axis=axis, cohort_variable=xvar, k_cohorts=int(ok.sum()),
                        spearman_rho=float(rr.statistic), spearman_p=float(rr.pvalue),
                        pearson_r_fisherz=float(pp.statistic), pearson_p=float(pp.pvalue)))
E = pd.DataFrame(eco)
for xv, lab in (("rho_confound", "rho(source,CAF) x rho(target,CAF)  [the confound term]"),
                ("cafA_sd",      "within-cohort SD of the CAF score  [stromal heterogeneity]"),
                ("cafA_mean",    "cohort MEAN CAF score              [stromal level]")):
    print(f"\n>>> across-cohort predictor: {lab}")
    print(E[E.cohort_variable == xv].drop(columns="cohort_variable")
           .to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# how much of each cohort's TF->collagen rho is explained by the confound term?
print("\n--- per-cohort confound decomposition for ETS1 -> COL1A1 ---")
d = R[R.axis == "ETS1 -> COL1A1"].merge(C, on="cohort").sort_values("cafA_sd", ascending=False)
print(f"{'cohort':46s} {'n':>5s} {'CAFsd':>6s} {'r(ETS1,CAF)':>11s} {'r(COL,CAF)':>10s} "
      f"{'product':>8s} {'rho':>7s} {'partial':>8s}")
for _, x in d.iterrows():
    print(f"  {x.cohort[:44]:44s} {int(x.n):5d} {x.cafA_sd:6.3f} {x.rho_source_cafA:+11.3f} "
          f"{x.rho_target_cafA:+10.3f} {x.rho_confound:+8.3f} {x.rho:+7.3f} {x.rho_partial_cafA:+8.3f}")

print("\n--- cohort stromal content (CAF-A mean, pooled z-scale), sorted ---")
for _, x in C.sort_values("cafA_mean", ascending=False).iterrows():
    print(f"  {x.cohort[:44]:44s} n={int(x.n_gene):4d} CAF-A={x.cafA_mean:+.3f} "
          f"CAF-B={x.cafB_mean:+.3f} COL1A1mean={x.COL1A1_mean:.2f}")

# ------------------------------- DOES THE AXIS SURVIVE CAF ADJUSTMENT? ------
# This is the non-tautological comparison. (Note: rho_xy is algebraically tied to
# rho_confound via rho_xy = rho_xz*rho_yz + rho_partial*sqrt((1-rxz^2)(1-ryz^2)),
# so the rho_confound predictor above is partly definitional; cafA_sd is not, and
# the retention statistics below are not.)
print("\n" + "="*78)
print("SURVIVAL OF EACH AXIS AFTER WITHIN-COHORT CAF ADJUSTMENT (across cohorts)")
print("="*78)
surv = []
for axis in [f"{tf} -> {col}" for tf in TFS for col in COLS] + \
            [f"{m} -> {c}" for m in MIRS for c in COLS]:
    mp = "arm_averaged" if axis.startswith("hsa-") else ""
    d = R[(R.axis == axis) & (R.mir_mapping == mp)]
    if len(d) < 8: continue
    exp_sign = -1 if axis.startswith("hsa-") else +1     # miRNA repress, TF activate
    raw = d.rho.values.astype(float); par = d.rho_partial_cafA.values.astype(float)
    ok = np.isfinite(raw) & np.isfinite(par); raw, par = raw[ok], par[ok]
    k = len(raw)
    wil = stats.wilcoxon(np.arctanh(raw), np.arctanh(par))
    surv.append(dict(axis=axis, k_cohorts=k,
                     median_rho=float(np.median(raw)),
                     median_rho_partial_cafA=float(np.median(par)),
                     n_expected_sign_raw=int((np.sign(raw) == exp_sign).sum()),
                     n_expected_sign_partial=int((np.sign(par) == exp_sign).sum()),
                     pct_retained_median=float(100 * np.median(par) / np.median(raw))
                                          if np.median(raw) != 0 else np.nan,
                     wilcoxon_raw_vs_partial_p=float(wil.pvalue)))
SV = pd.DataFrame(surv)
print(SV.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# --------------------------------------------------------------- write out --
R.to_csv(os.path.join(OUT, "pancancer_axis.csv"), index=False)
C.to_csv(os.path.join(OUT, "pancancer_cohort_summary.csv"), index=False)
S.to_csv(os.path.join(OUT, "pancancer_axis_consistency.csv"), index=False)
E.to_csv(os.path.join(OUT, "pancancer_ecological.csv"), index=False)
SV.to_csv(os.path.join(OUT, "pancancer_caf_survival.csv"), index=False)
print(f"\nWROTE pancancer_axis.csv            ({len(R)} rows)")
print(f"WROTE pancancer_cohort_summary.csv  ({len(C)} rows)")
print(f"WROTE pancancer_axis_consistency.csv({len(S)} rows)")
print(f"WROTE pancancer_ecological.csv      ({len(E)} rows)")
print("DONE PART B")
