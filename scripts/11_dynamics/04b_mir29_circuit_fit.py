#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04b_mir29_circuit_fit.py  --  Part D: parameterise one concrete circuit from the data.

CIRCUIT (all four edges are in the manuscript's own network, data/canonical_edges.tsv):
    NFKB1  -| hsa-miR-29a/b/c     TransmiR, Repression, PMID 18977326 / 20564213
    miR-29 -| COL1A1              miRNA_target, sign -1
    miR-29 -| NFKB1               miRNA_target, sign -1   <- reciprocal arm => COMPOSITE
    NFKB1  ->  COL1A1             TRRUST mode "Unknown" (PMID 11277268); the sign used
                                  here is taken from the data (Spearman rho = +0.129,
                                  p = 1.7e-5, n = 1097), not assumed
This is a composite FFL of type C2 with a DOUBLE-NEGATIVE (mutual repression)
TF<->miRNA arm, i.e. a toggle switch embedded in a feed-forward loop.

WHAT IS FITTED AND WHAT IS NOT.  Cross-sectional tumour data identify the circuit's
STEADY-STATE input-output relation, and nothing else.  Time constants cannot be
identified from cross-sectional data and are therefore taken from published
measurements and swept:
    COL1A1 mRNA half-life  1.5 h (quiescent) / 24 h (activated myofibroblast)
                           Stefanovic et al. 1997 Mol Cell Biol doi:10.1128/mcb.17.9.5201
    mature miRNA half-life 119 h  Gantier et al. 2011 NAR doi:10.1093/nar/gkr148
    protein half-life      46 h (median mammalian)
                           Schwanhausser et al. 2011 Nature doi:10.1038/nature10098
This separation is stated in every output table.
"""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D
from scipy.optimize import least_squares
from scipy import stats

RES = "/path/to/revision/results/v3"
rng = np.random.default_rng(20260908)

T = pd.read_csv(f"{RES}/dynamics_mir29_tcga_measurements.csv")
C = pd.read_csv(f"{RES}/dynamics_mir29_cptac_measurements.csv")
print(f"TCGA n={len(T)}  CPTAC n={len(C)}")

med_m = float(np.median(T.miR29_family))
med_x = float(np.median(T.NFKB1))
T["m"] = 2.0 ** (T.miR29_family - med_m)          # miR-29 in units of the cohort median
T["dx"] = T.NFKB1 - med_x

# ---------------------------------------------------------------- the model --
def suppress(m, lam, K, n):
    """Total post-transcriptional removal factor 1 + lam*q(m); q = Hill occupancy."""
    q = m ** n / (K ** n + m ** n)
    return 1.0 + lam * q

def predict(p, m, dx):
    c0, a, lam, K, n = p
    return c0 + a * dx - np.log2(suppress(m, lam, K, n))

def resid(p, m, dx, y):
    return predict(p, m, dx) - y

def fit(m, dx, y, p0=(15.5, 0.5, 3.0, 1.0, 1.5)):
    lo = [0.0, -5.0, 0.0, 0.05, 0.5]
    hi = [30.0, 5.0, 500.0, 20.0, 4.0]
    r = least_squares(resid, p0, bounds=(lo, hi), args=(m, dx, y), max_nfev=20000)
    return r

y = T.COL1A1.values; m = T.m.values; dx = T.dx.values
r = fit(m, dx, y)
c0, a, lam, K, n = r.x
res = resid(r.x, m, dx, y)
sse = float(np.sum(res ** 2)); sst = float(np.sum((y - y.mean()) ** 2))
R2 = 1 - sse / sst
N = len(y); k = 5
AIC = N * np.log(sse / N) + 2 * k
print(f"\nFITTED (COL1A1 mRNA, TCGA n={N}):")
print(f"  c0={c0:.4f}  a(NFKB1 drive)={a:.4f}  lam={lam:.4f}  K={K:.4f}  n={n:.4f}  R2={R2:.4f}")

# ---- comparison models ------------------------------------------------------
X1 = np.column_stack([np.ones(N), dx])
b1, *_ = np.linalg.lstsq(X1, y, rcond=None)
sse1 = float(np.sum((y - X1 @ b1) ** 2)); AIC1 = N * np.log(sse1 / N) + 2 * 2
X2 = np.column_stack([np.ones(N), dx, T.miR29_family.values - med_m])
b2, *_ = np.linalg.lstsq(X2, y, rcond=None)
sse2 = float(np.sum((y - X2 @ b2) ** 2)); AIC2 = N * np.log(sse2 / N) + 2 * 3
print(f"  null (TF only)          R2={1-sse1/sst:.4f}  AIC={AIC1:.1f}")
print(f"  log-linear (TF + miR)   R2={1-sse2/sst:.4f}  AIC={AIC2:.1f}  miR slope={b2[2]:.4f}")
print(f"  Hill circuit model      R2={R2:.4f}  AIC={AIC:.1f}")

# ---- bootstrap CIs ----------------------------------------------------------
B = 2000
bs = np.empty((B, 5))
for i in range(B):
    idx = rng.integers(0, N, N)
    try:
        bs[i] = fit(m[idx], dx[idx], y[idx], p0=r.x).x
    except Exception:
        bs[i] = np.nan
bs = bs[np.isfinite(bs).all(axis=1)]
ci = np.nanpercentile(bs, [2.5, 97.5], axis=0)
names = ["c0", "a_NFKB1", "lambda", "K", "n"]
fit_rows = [dict(parameter=nm, estimate=float(v), ci_lo=float(ci[0, j]), ci_hi=float(ci[1, j]),
                 units=("log2 units" if nm == "c0" else
                        "log2 COL1A1 per log2 NFKB1" if nm == "a_NFKB1" else
                        "max fold-acceleration of COL1A1 mRNA decay" if nm == "lambda" else
                        "miR-29 fold of cohort median" if nm == "K" else "dimensionless"))
            for j, (nm, v) in enumerate(zip(names, r.x))]

# ---------------------------------------------- identifiability diagnostics --
# K rails at its upper bound whenever the Hill function stays in its linear regime
# over the observed miR-29 range.  We therefore ALSO fit a version with K confined
# to the observed range, and we report the local logarithmic sensitivity, which is
# the combination that cross-sectional data actually identify.
def fit_bounded(m_, dx_, y_, Klo, Khi, p0):
    lo = [0.0, -5.0, 0.0, Klo, 0.5]; hi = [30.0, 5.0, 500.0, Khi, 4.0]
    return least_squares(resid, p0, bounds=(lo, hi), args=(m_, dx_, y_), max_nfev=20000)

obs_lo, obs_hi = float(np.percentile(m, 1)), float(np.percentile(m, 99))
rb = fit_bounded(m, dx, y, max(obs_lo, 0.1), min(obs_hi, 5.0), (17.4, 0.31, 3.0, 1.0, 1.5))
c0b, ab, lamb, Kb, nb = rb.x
sseb = float(np.sum(resid(rb.x, m, dx, y) ** 2))
print(f"\n  CONSTRAINED fit (K forced into the observed miR-29 range "
      f"[{max(obs_lo,0.1):.2f},{min(obs_hi,5.0):.2f}]): "
      f"lam={lamb:.3f} K={Kb:.3f} n={nb:.3f} R2={1-sseb/sst:.4f}")

def local_sens(lam_, K_, n_, m0=1.0, h=0.01):
    return -(np.log2(suppress(m0 * (1 + h), lam_, K_, n_)) -
             np.log2(suppress(m0, lam_, K_, n_))) / np.log2(1 + h)
S_local = local_sens(lam, K, n)
S_local_b = local_sens(lamb, Kb, nb)
bs_S = np.array([local_sens(pp[2], pp[3], pp[4]) for pp in bs])
print(f"  local logarithmic sensitivity at the cohort median "
      f"d log2 COL1A1 / d log2 miR-29 = {S_local:.4f} "
      f"[95% CI {np.nanpercentile(bs_S,2.5):.3f} to {np.nanpercentile(bs_S,97.5):.3f}]; "
      f"constrained fit gives {S_local_b:.4f}; log-linear regression gives {b2[2]:.4f}")
print(f"  K estimate hit its upper bound: {abs(K-20.0)<1e-3}  -> lambda and K are NOT "
      f"separately identifiable from cross-sectional data; the identified quantity is the "
      f"local slope above.")

# ---- is the reciprocal composite arm (miR-29 -| NFKB1) supported by the data? --
from scipy.stats import spearmanr
rec = spearmanr(T.miR29_family, T.NFKB1)
print(f"\n  RECIPROCAL ARM CHECK  miR-29 vs NFKB1: Spearman rho={rec.statistic:+.4f} "
      f"p={rec.pvalue:.3g} (n={len(T)}). The annotated edge is REPRESSIVE; the observed "
      f"association has the OPPOSITE sign.")

# ---------------------------------------------------- what the circuit predicts
q1 = 1.0 ** n / (K ** n + 1.0 ** n)               # occupancy at the cohort median
F1 = 1.0 + lam * q1                                # removal factor at the median
Fmax = 1.0 + lam                                   # removal factor at saturation
max_fold_suppression = Fmax / F1                   # best achievable from the median

def fold_to_reduce(factor, lam_, K_, n_):
    """miR-29 fold-change (relative to cohort median) that divides COL1A1 by `factor`."""
    q1_ = 1.0 / (K_ ** n_ + 1.0)
    F1_ = 1.0 + lam_ * q1_
    Fneed = factor * F1_
    qneed = (Fneed - 1.0) / lam_
    if not (0 < qneed < 1):
        return np.nan
    return K_ * (qneed / (1.0 - qneed)) ** (1.0 / n_)

f_half = fold_to_reduce(2.0, lam, K, n)
bs_fhalf = np.array([fold_to_reduce(2.0, p[2], p[3], p[4]) for p in bs])
bs_maxsupp = np.array([(1 + p[2]) / (1 + p[2] / (p[3] ** p[4] + 1)) for p in bs])
print(f"\n  removal factor at cohort median  F(1)   = {F1:.3f}")
print(f"  max achievable suppression from median  = {max_fold_suppression:.3f}-fold "
      f"[95% CI {np.nanpercentile(bs_maxsupp,2.5):.2f}-{np.nanpercentile(bs_maxsupp,97.5):.2f}]")
print(f"  miR-29 fold needed to HALVE COL1A1 mRNA = {f_half:.3f} "
      f"[95% CI {np.nanpercentile(bs_fhalf,2.5):.2f}-{np.nanpercentile(bs_fhalf,97.5):.2f}] "
      f"({np.isnan(bs_fhalf).mean()*100:.1f}% of bootstraps: unattainable)")

# ---- model-free check -------------------------------------------------------
sl = b2[2]
f_half_loglin = 2.0 ** (1.0 / abs(sl)) if sl < 0 else np.nan
print(f"  model-free (log-linear slope {sl:.3f}): fold needed = {f_half_loglin:.2f}")

# ---- indirect (composite) path: miR-29 -| NFKB1 -> COL1A1 -------------------
lm_nf = stats.linregress(T.miR29_family - med_m, T.NFKB1)
indirect = a * lm_nf.slope
print(f"\n  miR-29 -| NFKB1 slope e = {lm_nf.slope:.4f} (p={lm_nf.pvalue:.3g})")
print(f"  indirect path a*e = {indirect:.4f} log2 COL1A1 per log2 miR-29")
print(f"  direct post-transcriptional slope at the median = "
      f"{-(np.log2(suppress(1.02**1, lam, K, n)) - np.log2(suppress(1.0, lam, K, n)))/np.log2(1.02):.4f}")

# ---------------------------------------------------------- protein-level fit
Cc = C.dropna(subset=["COL1A1_protein", "miR29_family_3p", "NFKB1_mRNA"]).copy()
med_mp = float(np.median(Cc.miR29_family_3p)); med_xp = float(np.median(Cc.NFKB1_mRNA))
mp = 2.0 ** (Cc.miR29_family_3p - med_mp); dxp = (Cc.NFKB1_mRNA - med_xp).values
yp = Cc.COL1A1_protein.values
rp = fit(mp.values, dxp, yp, p0=(np.median(yp), 0.2, 3.0, 1.0, 1.5))
c0p, ap, lamp, Kp, np_ = rp.x
ssep = float(np.sum(resid(rp.x, mp.values, dxp, yp) ** 2))
sstp = float(np.sum((yp - yp.mean()) ** 2))
f_half_prot = fold_to_reduce(2.0, lamp, Kp, np_)
print(f"\nFITTED (COL1A1 protein, CPTAC n={len(Cc)}): lam={lamp:.3f} K={Kp:.3f} n={np_:.3f} "
      f"R2={1-ssep/sstp:.3f}; fold to halve = {f_half_prot:.2f}")

# ============================================================ DYNAMICS ========
# Dimensional scenarios.  tau = COL1A1 mRNA lifetime = t_half / ln2.
SCEN = [
    dict(state="activated myofibroblast", mRNA_half_life_h=24.0,
         src="Stefanovic 1997 doi:10.1128/mcb.17.9.5201"),
    dict(state="quiescent stellate/fibroblast", mRNA_half_life_h=1.5,
         src="Stefanovic 1997 doi:10.1128/mcb.17.9.5201"),
    dict(state="median mammalian mRNA", mRNA_half_life_h=9.0,
         src="Schwanhausser 2011 doi:10.1038/nature10098"),
]
PROT_HL = [12.0, 46.0, 168.0]     # h: fast turnover / median mammalian / long-lived ECM pool
MIR_HL = 119.0                    # h, Gantier 2011

dyn_rows = []
for sc in SCEN:
    tau_h = sc["mRNA_half_life_h"] / np.log(2)          # hours per dimensionless time unit
    gm = sc["mRNA_half_life_h"] / MIR_HL                # alpha_miR / alpha_mRNA
    for phl in PROT_HL:
        gp = sc["mRNA_half_life_h"] / phl
        # --- isolated post-transcriptional arm (analytic + numerical) --------
        for fold in (1.0, 1.5, 2.0, f_half if np.isfinite(f_half) else np.nan, 3.0, 5.0, 10.0):
            if not np.isfinite(fold):
                continue
            F0, F1n = suppress(1.0, lam, K, n), suppress(fold, lam, K, n)
            # --- exact solution of the linear 2-compartment relaxation ----------
            #   y' = 1 - F1n*y            y(0) = 1/F0     y_inf = 1/F1n
            #   p' = gp*(y - p)           p(0) = 1/F0     p_inf = 1/F1n
            y0 = 1.0 / F0; yinf = 1.0 / F1n
            if abs(F1n - gp) < 1e-9:
                F1n = F1n * (1 + 1e-6)
            A = (y0 - yinf) * gp / (gp - F1n)
            Bc = (y0 - yinf) - A
            pfun = lambda tt: yinf + A * np.exp(-F1n * tt) + Bc * np.exp(-gp * tt)
            yfun = lambda tt: yinf + (y0 - yinf) * np.exp(-F1n * tt)
            if abs(yinf - y0) < 1e-12:
                t_half_mRNA_h = np.nan; t_half_prot_h = np.nan; t90_prot_h = np.nan
            else:
                t_half_mRNA_h = np.log(2) / F1n / np.log(2) * 0 + (np.log(2) / F1n) * tau_h
                from scipy.optimize import brentq
                def solve_frac(f_, fun):
                    tgt = y0 + f_ * (yinf - y0)
                    g = lambda tt: fun(tt) - tgt
                    hi_t = 1.0
                    while g(hi_t) * g(1e-9) > 0 and hi_t < 1e5:
                        hi_t *= 2
                    return brentq(g, 1e-9, hi_t) if hi_t < 1e5 else np.nan
                t_half_prot_h = solve_frac(0.5, pfun) * tau_h
                t90_prot_h = solve_frac(0.9, pfun) * tau_h
            dyn_rows.append(dict(
                cell_state=sc["state"], mRNA_half_life_h=sc["mRNA_half_life_h"],
                protein_half_life_h=phl, miRNA_half_life_h=MIR_HL,
                tau_hours_per_time_unit=tau_h,
                miR29_fold_over_cohort_median=fold,
                steady_state_COL1A1_mRNA_rel=F0 / F1n,
                steady_state_COL1A1_mRNA_pct_drop=100 * (1 - F0 / F1n),
                t_half_COL1A1_mRNA_hours=t_half_mRNA_h,
                t_half_COL1A1_protein_hours=t_half_prot_h,
                t90_COL1A1_protein_hours=t90_prot_h,
                source_mRNA_half_life=sc["src"]))
DY = pd.DataFrame(dyn_rows)
DY.to_csv(f"{RES}/dynamics_mir29_predictions.csv", index=False)

pd.DataFrame(fit_rows + [
    dict(parameter="R2_circuit_model", estimate=R2, ci_lo=np.nan, ci_hi=np.nan, units="fraction"),
    dict(parameter="R2_TF_only", estimate=1 - sse1 / sst, ci_lo=np.nan, ci_hi=np.nan, units="fraction"),
    dict(parameter="R2_loglinear", estimate=1 - sse2 / sst, ci_lo=np.nan, ci_hi=np.nan, units="fraction"),
    dict(parameter="AIC_circuit_model", estimate=AIC, ci_lo=np.nan, ci_hi=np.nan, units="AIC"),
    dict(parameter="AIC_TF_only", estimate=AIC1, ci_lo=np.nan, ci_hi=np.nan, units="AIC"),
    dict(parameter="AIC_loglinear", estimate=AIC2, ci_lo=np.nan, ci_hi=np.nan, units="AIC"),
    dict(parameter="removal_factor_at_cohort_median", estimate=F1, ci_lo=np.nan, ci_hi=np.nan, units="fold"),
    dict(parameter="max_achievable_suppression_from_median", estimate=max_fold_suppression,
         ci_lo=float(np.nanpercentile(bs_maxsupp, 2.5)), ci_hi=float(np.nanpercentile(bs_maxsupp, 97.5)), units="fold"),
    dict(parameter="miR29_fold_to_halve_COL1A1_mRNA", estimate=f_half,
         ci_lo=float(np.nanpercentile(bs_fhalf, 2.5)), ci_hi=float(np.nanpercentile(bs_fhalf, 97.5)), units="fold of cohort median"),
    dict(parameter="miR29_fold_to_halve_COL1A1_mRNA_loglinear_check", estimate=f_half_loglin,
         ci_lo=np.nan, ci_hi=np.nan, units="fold of cohort median"),
    dict(parameter="miR29_fold_to_halve_COL1A1_protein_CPTAC", estimate=f_half_prot,
         ci_lo=np.nan, ci_hi=np.nan, units="fold of cohort median"),
    dict(parameter="lambda_protein_CPTAC", estimate=lamp, ci_lo=np.nan, ci_hi=np.nan, units="fold"),
    dict(parameter="K_protein_CPTAC", estimate=Kp, ci_lo=np.nan, ci_hi=np.nan, units="fold of median"),
    dict(parameter="n_protein_CPTAC", estimate=np_, ci_lo=np.nan, ci_hi=np.nan, units="dimensionless"),
    dict(parameter="indirect_path_slope_miR29_NFKB1_COL1A1", estimate=indirect, ci_lo=np.nan, ci_hi=np.nan,
         units="log2 COL1A1 per log2 miR-29"),
    dict(parameter="slope_miR29_to_NFKB1", estimate=lm_nf.slope, ci_lo=np.nan, ci_hi=np.nan,
         units="log2 NFKB1 per log2 miR-29"),
    dict(parameter="loglinear_total_slope_miR29_COL1A1", estimate=b2[2], ci_lo=np.nan, ci_hi=np.nan,
         units="log2 COL1A1 per log2 miR-29"),
    dict(parameter="local_log_sensitivity_at_median", estimate=S_local,
         ci_lo=float(np.nanpercentile(bs_S, 2.5)), ci_hi=float(np.nanpercentile(bs_S, 97.5)),
         units="d log2 COL1A1 / d log2 miR-29"),
    dict(parameter="local_log_sensitivity_constrained_fit", estimate=S_local_b,
         ci_lo=np.nan, ci_hi=np.nan, units="d log2 COL1A1 / d log2 miR-29"),
    dict(parameter="lambda_constrained_K_in_observed_range", estimate=lamb,
         ci_lo=np.nan, ci_hi=np.nan, units="max fold-acceleration of decay"),
    dict(parameter="K_constrained", estimate=Kb, ci_lo=np.nan, ci_hi=np.nan,
         units="miR-29 fold of cohort median"),
    dict(parameter="n_constrained", estimate=nb, ci_lo=np.nan, ci_hi=np.nan, units="dimensionless"),
    dict(parameter="R2_constrained", estimate=1 - sseb / sst, ci_lo=np.nan, ci_hi=np.nan, units="fraction"),
    dict(parameter="K_hit_upper_bound_nonidentifiable", estimate=float(abs(K - 20.0) < 1e-3),
         ci_lo=np.nan, ci_hi=np.nan, units="1 = yes"),
    dict(parameter="reciprocal_arm_miR29_vs_NFKB1_spearman_rho", estimate=float(rec.statistic),
         ci_lo=np.nan, ci_hi=np.nan, units="rho (annotated edge is repressive; observed sign is opposite)"),
    dict(parameter="reciprocal_arm_miR29_vs_NFKB1_p", estimate=float(rec.pvalue),
         ci_lo=np.nan, ci_hi=np.nan, units="p"),
    dict(parameter="observed_miR29_fold_range_p1_p99_lo", estimate=obs_lo, ci_lo=np.nan, ci_hi=np.nan,
         units="fold of cohort median"),
    dict(parameter="observed_miR29_fold_range_p1_p99_hi", estimate=obs_hi, ci_lo=np.nan, ci_hi=np.nan,
         units="fold of cohort median"),
]).to_csv(f"{RES}/dynamics_mir29_fit.csv", index=False)

# steady-state dose-response curve for the figure
folds = np.logspace(-2, 2, 200)
pd.DataFrame(dict(miR29_fold_over_median=folds,
                  COL1A1_mRNA_rel=suppress(1.0, lam, K, n) / suppress(folds, lam, K, n),
                  COL1A1_protein_rel=suppress(1.0, lamp, Kp, np_) / suppress(folds, lamp, Kp, np_))) \
  .to_csv(f"{RES}/dynamics_mir29_doseresponse.csv", index=False)

print("\nPredictions written. Head:")
print(DY.head(12).to_string(index=False))
