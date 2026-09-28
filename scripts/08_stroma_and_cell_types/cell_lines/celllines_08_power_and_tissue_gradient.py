#!/usr/bin/env python3
"""
celllines_08_power_and_tissue_gradient.py

Two things needed to read the cell-line null correctly.

(1) POWER.  What effect size can n = 50 (miRNA) / n = 69-70 (mRNA) breast lines actually
    exclude?  Reports the minimum detectable |rho| at 80% power and the exact 95% CI of
    each observed cell-line correlation, so "no effect" is never claimed beyond the data.

(2) TISSUE GRADIENT.  In TCGA-BRCA, recompute the same correlations inside tertiles of a
    4-marker CAF score (DCN, LUM, FAP, THY1).  If a relationship is tumour-cell-intrinsic
    it should be present in the CAF-low tertile; if it is compositional it should scale
    with stromal content.  This is the in-tissue counterpart of the cell-line test.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")

# ------------------------------------------------------------------ 1. power
def mde_rho(n, power=0.80, alpha=0.05):
    """minimum |rho| detectable at given power via Fisher z."""
    za = stats.norm.ppf(1 - alpha / 2)
    zb = stats.norm.ppf(power)
    z = (za + zb) / np.sqrt(n - 3)
    return float(np.tanh(z))


def ci_rho(r, n, alpha=0.05):
    if n < 5:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    za = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - za * se)), float(np.tanh(z + za * se))


pw = [dict(n=n, MDE_rho_80pct_power=mde_rho(n), MDE_rho_90pct_power=mde_rho(n, 0.90),
           note=note)
      for n, note in [(50, "breast lines with miRNA + mRNA"),
                      (69, "breast cancer lines with mRNA"),
                      (70, "all breast lineage lines with mRNA"),
                      (47, "breast lines in GDSC2"),
                      (22, "breast lines in PRISM secondary"),
                      (1066, "TCGA-BRCA paired tumours"),
                      (1097, "TCGA-BRCA primary tumours")]]
pd.DataFrame(pw).to_csv(os.path.join(RES, "celllines_F_power.csv"), index=False)
print(pd.DataFrame(pw).round(4).to_string(index=False))

# attach CIs to the main correlation table
d = pd.read_csv(os.path.join(RES, "celllines_B_correlations.csv"))
lo, hi = zip(*[ci_rho(r, n) if np.isfinite(r) and np.isfinite(n) else (np.nan, np.nan)
               for r, n in zip(d["rho"], d["n"])])
d["rho_CI95_lo"], d["rho_CI95_hi"] = lo, hi
d.to_csv(os.path.join(RES, "celllines_B_correlations.csv"), index=False)
key = d[(d.subset.isin(["breast_cancer_lines_only", "miRNA_matched_cancer_only"]))
        & (d.regulator.isin(["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "ETS1", "NFKB1",
                             "RELA", "SP1", "COL1A1"]))]
print("\ncell-line estimates with 95% CI:")
print(key[["regulator", "target", "n", "rho", "rho_CI95_lo", "rho_CI95_hi", "p"]]
      .round(4).to_string(index=False))

# ------------------------------------------------------------------ 2. tissue gradient
tg = pd.read_csv(os.path.join(RES, "_tcga_genes_tumour.tsv"), sep="\t", index_col=0).T
tgp = pd.read_csv(os.path.join(RES, "_tcga_genes_paired.tsv"), sep="\t", index_col=0).T
tm = pd.read_csv(os.path.join(RES, "_tcga_mirs_tumour.tsv"), sep="\t", index_col=0).T
zc = lambda df: (df - df.mean()) / df.std(ddof=1)
CAFG = ["DCN", "LUM", "FAP", "THY1"]
caf_g = zc(tg[CAFG]).mean(axis=1)
caf_p = zc(tgp[CAFG]).mean(axis=1)
print(f"\nTCGA CAF score: gene panel {CAFG}; n_gene={len(caf_g)} n_paired={len(caf_p)}")
print("CAF score vs COL1A1 rho = %.3f ; vs COL3A1 rho = %.3f" %
      (stats.spearmanr(caf_g, tg["COL1A1"])[0], stats.spearmanr(caf_g, tg["COL3A1"])[0]))

rows = []
for label, X, Y, caf in [("TF_target", tg, tg, caf_g), ("miRNA_target", tm, tgp, caf_p)]:
    q = np.quantile(caf.values, [1 / 3, 2 / 3])
    grp = np.digitize(caf.values, q)          # 0 = CAF-low, 1 = mid, 2 = CAF-high
    regs = (["ETS1", "NFKB1", "RELA", "SP1", "MKL1", "STAT6", "TFAP2A", "MYB"]
            if label == "TF_target" else
            ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-miR-101", "hsa-miR-21"])
    for r_ in regs:
        if r_ not in X.columns:
            continue
        for c in ["COL1A1", "COL3A1"]:
            all_r, all_p = stats.spearmanr(X[r_].values, Y[c].values)
            rec = dict(pair_type=label, regulator=r_, target=c, n_all=len(caf),
                       rho_all=float(all_r), p_all=float(all_p))
            for g, name in [(0, "CAFlow"), (1, "CAFmid"), (2, "CAFhigh")]:
                m = grp == g
                rr, pp = stats.spearmanr(X[r_].values[m], Y[c].values[m])
                rec[f"n_{name}"] = int(m.sum())
                rec[f"rho_{name}"] = float(rr)
                rec[f"p_{name}"] = float(pp)
            rows.append(rec)
t = pd.DataFrame(rows)
t.to_csv(os.path.join(RES, "celllines_F_tcga_CAF_tertiles.csv"), index=False)
print("\nTCGA correlations within CAF tertiles:")
print(t[["pair_type", "regulator", "target", "rho_all", "rho_CAFlow", "rho_CAFmid",
         "rho_CAFhigh", "p_CAFlow"]].round(4).to_string(index=False))

# ------------------------------------------------------------------ 3. collagen dynamic range
print("\nCOL1A1/COL3A1 dispersion, TCGA bulk vs cell lines "
      "(sd of log2 values on each platform's own scale):")
E = pd.read_csv(os.path.join(RES, "celllines_expr_breast.tsv.gz"), sep="\t", index_col=0)
disp = []
for c in ["COL1A1", "COL3A1"]:
    disp.append(dict(gene=c, dataset="TCGA-BRCA bulk (log2 norm_count+1)",
                     n=len(tg), sd=float(tg[c].std(ddof=1)),
                     IQR=float(np.subtract(*np.percentile(tg[c], [75, 25]))),
                     median=float(tg[c].median())))
    disp.append(dict(gene=c, dataset="DepMap breast lines (log2 TPM+1)",
                     n=len(E), sd=float(E[c].std(ddof=1)),
                     IQR=float(np.subtract(*np.percentile(E[c], [75, 25]))),
                     median=float(E[c].median())))
disp = pd.DataFrame(disp)
disp.to_csv(os.path.join(RES, "celllines_F_dispersion.csv"), index=False)
print(disp.round(3).to_string(index=False))
print("DONE")
