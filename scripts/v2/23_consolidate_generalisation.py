#!/usr/bin/env python3
"""
Consolidation for the subtype / pan-cancer / CPTAC-pan-cancer task.

(1) Is TCGA-BRCA an OUTLIER among the 32 TCGA cohorts for each axis, or is it
    an ordinary member?  -> answers "breast-specific or general".
(2) Random- and fixed-effects meta-analysis of the CPTAC protein-level
    miR-29 -> collagen correlations across the 10 CPTAC cohorts, with Cochran's
    Q / I^2 heterogeneity.
(3) One tidy verdict table.
"""
import os
import numpy as np, pandas as pd
from scipy import stats

BASE = "/path/to/revision"
OUT  = os.path.join(BASE, "results/v2")

R  = pd.read_csv(os.path.join(OUT, "pancancer_axis.csv"))
CP = pd.read_csv(os.path.join(OUT, "cptac_pancancer_mir29.csv"))
ST = pd.read_csv(os.path.join(OUT, "subtype_stratified.csv"))

print("="*78); print("(1) IS BREAST AN OUTLIER AMONG TCGA COHORTS?"); print("="*78)
BR = "breast invasive carcinoma"
rows = []
for (axis, mp), d in R.groupby(["axis", "mir_mapping"], dropna=False):
    mp = "" if pd.isna(mp) else mp
    if mp not in ("", "arm_averaged"): continue
    d = d.dropna(subset=["rho"])
    if BR not in set(d.cohort) or len(d) < 8: continue
    rb = float(d.loc[d.cohort == BR, "rho"].iloc[0])
    others = d.loc[d.cohort != BR, "rho"].values.astype(float)
    pct = float(100 * (others < rb).mean())
    # is BRCA's rho consistent with the other cohorts' distribution?
    z_vs_dist = (np.arctanh(rb) - np.arctanh(others).mean()) / np.arctanh(others).std(ddof=1)
    rows.append(dict(axis=axis, k_other_cohorts=len(others), rho_BRCA=rb,
                     median_other=float(np.median(others)),
                     min_other=float(others.min()), max_other=float(others.max()),
                     BRCA_percentile=pct, z_BRCA_vs_other_cohorts=float(z_vs_dist),
                     p_two_sided=float(2 * stats.norm.sf(abs(z_vs_dist)))))
OUTL = pd.DataFrame(rows).sort_values("axis")
print(OUTL.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

print("\n" + "="*78)
print("(2) META-ANALYSIS OF THE CPTAC PROTEIN-LEVEL miR-29 -> COLLAGEN RESULT")
print("="*78)

def meta(r, n):
    """DerSimonian-Laird random effects on Fisher-z, plus fixed effect."""
    ok = np.isfinite(r) & (n > 3)
    r, n = np.asarray(r)[ok], np.asarray(n)[ok]
    k = len(r)
    z = np.arctanh(r); w = (n - 3).astype(float)
    zf = (w * z).sum() / w.sum(); se_f = np.sqrt(1 / w.sum())
    Q = (w * (z - zf) ** 2).sum(); df = k - 1
    pQ = stats.chi2.sf(Q, df) if df > 0 else np.nan
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wr = 1 / (1 / w + tau2)
    zr = (wr * z).sum() / wr.sum(); se_r = np.sqrt(1 / wr.sum())
    return dict(k=k, n_total=int(n.sum()),
                rho_fixed=float(np.tanh(zf)),
                p_fixed=float(2 * stats.norm.sf(abs(zf / se_f))),
                rho_random=float(np.tanh(zr)),
                ci_lo=float(np.tanh(zr - 1.96 * se_r)), ci_hi=float(np.tanh(zr + 1.96 * se_r)),
                p_random=float(2 * stats.norm.sf(abs(zr / se_r))),
                Q=float(Q), df=int(df), p_heterogeneity=float(pQ), I2_pct=float(I2),
                tau2=float(tau2))

mrows = []
for mp in ("arm_averaged", "dominant_arm"):
    for mir in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"]:
        for g in ["COL1A1", "COL3A1"]:
            d = CP[(CP.mir_mapping == mp) & (CP.mir == mir) & (CP.gene == g)]
            d = d.dropna(subset=["rho_protein"])
            if len(d) < 3: continue
            m = meta(d.rho_protein.values, d.n_protein.values)
            m.update(layer="protein", mapping=mp, axis=f"{mir} -> {g}",
                     n_neg=int((d.rho_protein < 0).sum()),
                     n_sig_neg=int(((d.rho_protein < 0) & (d.p_protein < .05)).sum()),
                     n_sig_pos=int(((d.rho_protein > 0) & (d.p_protein < .05)).sum()))
            mrows.append(m)
            dm = d.dropna(subset=["rho_mrna"])
            if len(dm) >= 3:
                m2 = meta(dm.rho_mrna.values, dm.n_mrna.values)
                m2.update(layer="mRNA", mapping=mp, axis=f"{mir} -> {g}",
                          n_neg=int((dm.rho_mrna < 0).sum()),
                          n_sig_neg=int(((dm.rho_mrna < 0) & (dm.p_mrna < .05)).sum()),
                          n_sig_pos=int(((dm.rho_mrna > 0) & (dm.p_mrna < .05)).sum()))
                mrows.append(m2)
META = pd.DataFrame(mrows)[["axis", "mapping", "layer", "k", "n_total", "n_neg",
                            "n_sig_neg", "n_sig_pos", "rho_fixed", "p_fixed",
                            "rho_random", "ci_lo", "ci_hi", "p_random",
                            "Q", "df", "p_heterogeneity", "I2_pct", "tau2"]]
print("\n-- arm-averaged, protein layer --")
print(META[(META.mapping == "arm_averaged") & (META.layer == "protein")]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print("\n-- arm-averaged, mRNA layer (same CPTAC samples) --")
print(META[(META.mapping == "arm_averaged") & (META.layer == "mRNA")]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print("\n-- dominant-arm (-3p), protein layer --")
print(META[(META.mapping == "dominant_arm") & (META.layer == "protein")]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# pooled over all 6 miR x collagen combinations, arm-averaged, protein
d = CP[(CP.mir_mapping == "arm_averaged") & (CP.gene.isin(["COL1A1", "COL3A1"]))].dropna(subset=["rho_protein"])
per_cohort = d.groupby("cohort").rho_protein.mean()
print("\nmean protein-level rho over the 6 miR-29 x (COL1A1,COL3A1) pairs, per cohort:")
for c, v in per_cohort.sort_values().items():
    print(f"   {c:8s} {v:+.3f}")
neg = int((per_cohort < 0).sum())
print(f"   => {neg}/{len(per_cohort)} CPTAC cohorts have a mean NEGATIVE miR-29->collagen "
      f"protein correlation; sign-test p={stats.binomtest(neg, len(per_cohort), 0.5, alternative='greater').pvalue:.4g}")

print("\n" + "="*78)
print("(3) PROTEIN-vs-mRNA AMPLIFICATION: BREAST-SPECIFIC OR GENERAL?")
print("="*78)
d = CP[(CP.mir_mapping == "arm_averaged") & CP.gene.isin(["COL1A1", "COL3A1"])].dropna(
        subset=["rho_protein", "rho_mrna"])
amp = []
for c, dd in d.groupby("cohort"):
    bigger = int((dd.rho_protein.abs() > dd.rho_mrna.abs()).sum())
    amp.append(dict(cohort=c, k_pairs=len(dd), n_stronger_at_protein=bigger,
                    mean_absrho_mrna=float(dd.rho_mrna.abs().mean()),
                    mean_absrho_protein=float(dd.rho_protein.abs().mean()),
                    wilcoxon_p=float(stats.wilcoxon(dd.rho_protein.abs(), dd.rho_mrna.abs()).pvalue)
                                if len(dd) >= 6 else np.nan))
AMP = pd.DataFrame(amp)
print(AMP.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
tot_b = int(AMP.n_stronger_at_protein.sum()); tot_k = int(AMP.k_pairs.sum())
print(f"\npooled over all 10 cohorts: |rho| stronger at protein in {tot_b}/{tot_k} pairs, "
      f"binomial p={stats.binomtest(tot_b, tot_k, 0.5).pvalue:.4g}")

OUTL.to_csv(os.path.join(OUT, "pancancer_brca_outlier.csv"), index=False)
META.to_csv(os.path.join(OUT, "cptac_pancancer_meta.csv"), index=False)
AMP.to_csv(os.path.join(OUT, "cptac_protein_vs_mrna_by_cohort.csv"), index=False)
print(f"\nWROTE pancancer_brca_outlier.csv ({len(OUTL)}), cptac_pancancer_meta.csv ({len(META)}), "
      f"cptac_protein_vs_mrna_by_cohort.csv ({len(AMP)})")
