#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04c_mir29_modelfree_prediction.py
=================================
The circuit fit in 04b gives "a 2.02-fold rise in miR-29 halves COL1A1 mRNA".
That number depends on a Hill form whose lambda and K are not separately
identifiable from cross-sectional data.  This script therefore states the SAME
prediction in a form that depends on no model at all -- the local logarithmic
slope of the measured relationship -- for every miRNA / target / layer / cohort
combination, with bootstrap confidence intervals, so the prediction can be
quoted with an honest interval and checked against an experiment directly.

  fold-change of miR-29 needed to halve the target = 2^(1 / |slope|)
  where slope = d log2(target) / d log2(miR-29)

Both the unadjusted slope and the slope adjusted for the TF drive (NFKB1) are
reported, because the circuit contains a TF arm onto the same target.

Out: results/v3/dynamics_mir29_modelfree_predictions.csv
"""
import numpy as np, pandas as pd

RES = "/path/to/revision/results/v3"
rng = np.random.default_rng(20260909)
B = 4000

T = pd.read_csv(f"{RES}/dynamics_mir29_tcga_measurements.csv")
C = pd.read_csv(f"{RES}/dynamics_mir29_cptac_measurements.csv")


def slope_ci(x, y, z=None):
    """OLS slope of y on x (optionally adjusting for z), with a percentile bootstrap CI."""
    ok = np.isfinite(x) & np.isfinite(y)
    if z is not None:
        ok &= np.isfinite(z)
    x, y = x[ok], y[ok]
    z = z[ok] if z is not None else None
    n = len(x)

    def fit(idx):
        X = np.column_stack([np.ones(len(idx)), x[idx]] + ([z[idx]] if z is not None else []))
        b, *_ = np.linalg.lstsq(X, y[idx], rcond=None)
        return b[1]
    est = fit(np.arange(n))
    bs = np.array([fit(rng.integers(0, n, n)) for _ in range(B)])
    return est, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), n, bs


def fold_to_halve(s):
    return np.where(s < 0, 2.0 ** (1.0 / np.abs(np.where(s < 0, s, np.nan))), np.nan)


rows = []
COMBOS = []
for m in ["hsa_miR_29a", "hsa_miR_29b", "hsa_miR_29c", "miR29_family"]:
    for g in ["COL1A1", "COL3A1"]:
        COMBOS.append(("TCGA-BRCA", "mRNA", m, g, T, m, g, "NFKB1"))
for m in ["hsa_miR_29a_3p", "hsa_miR_29b_3p", "hsa_miR_29c_3p", "miR29_family_3p"]:
    for g in ["COL1A1", "COL3A1"]:
        for layer in ["mRNA", "protein"]:
            COMBOS.append(("CPTAC-BRCA", layer, m, g, C, m, f"{g}_{layer}", "NFKB1_mRNA"))

for cohort, layer, mlab, glab, DF, mcol, gcol, tfcol in COMBOS:
    if mcol not in DF.columns or gcol not in DF.columns:
        continue
    x = DF[mcol].to_numpy(float)
    y = DF[gcol].to_numpy(float)
    z = DF[tfcol].to_numpy(float) if tfcol in DF.columns else None
    for adj, zz in (("unadjusted", None), (f"adjusted for {tfcol}", z)):
        est, lo, hi, n, bs = slope_ci(x, y, zz)
        f_est = fold_to_halve(np.array([est]))[0]
        fbs = fold_to_halve(bs)
        rows.append(dict(cohort=cohort, layer=layer, miRNA=mlab, target=glab,
                         adjustment=adj, n=n,
                         slope_log2_target_per_log2_miR=est,
                         slope_ci_lo=lo, slope_ci_hi=hi,
                         fold_miR29_to_halve_target=f_est,
                         fold_ci_lo=float(np.nanpercentile(fbs, 2.5)),
                         fold_ci_hi=float(np.nanpercentile(fbs, 97.5)),
                         pct_bootstraps_with_positive_slope=100 * float(np.mean(bs >= 0))))

R = pd.DataFrame(rows)
R.to_csv(f"{RES}/dynamics_mir29_modelfree_predictions.csv", index=False)
pd.set_option("display.width", 260)
print(R.round(3).to_string(index=False))
print()
key = R[(R.miRNA.str.startswith("miR29_family")) & (R.target == "COL1A1")]
print("HEADLINE (miR-29 family -> COL1A1):")
print(key.round(3).to_string(index=False))
