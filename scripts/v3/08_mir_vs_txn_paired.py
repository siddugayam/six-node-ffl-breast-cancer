#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_mir_vs_txn_paired.py
(a) PAIRED comparison, at identical parameter vectors, of each miRNA-mediated FFL
    against the transcriptional FFL of the same sign structure.  This is the direct
    answer to "what distinguishes these from conventional FFLs".
(b) The TITRATION THRESHOLD: with stoichiometric miRNA-target pairing the target's
    input-output curve acquires a threshold that no transcriptional repressor arm
    produces (Levine et al. 2007 doi:10.1371/journal.pbio.0050229;
    Mukherji et al. 2011 doi:10.1038/ng.905).  Solved exactly from the steady state.
Out: dynamics_mir_vs_txn_paired.csv, dynamics_titration_threshold.csv
"""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

RES = "/path/to/revision/results/v3"

# ---------------------------------------------------------------- (a) paired --
E = pd.read_csv(f"{RES}/dynamics_3node_metrics_ensemble.csv")
E = E[E.param_set >= 0]
METRICS = ["rel_T_ON", "rel_T_OFF", "overshoot_ratio", "adaptation", "is_pulse",
           "filter_index", "noise_transmission", "n_eff_protein", "n_eff_mRNA",
           "dyn_range_protein", "FCD_error", "is_FCD", "sign_sensitive_delay_index"]
PAIRS = [
    ("I1_MIR_AND", "I1_TXN_AND", "incoherent type 1, AND gate: miRNA arm vs protein-repressor arm"),
    ("I1_MIR_OR",  "I1_TXN_OR",  "incoherent type 1, OR gate: miRNA arm vs protein-repressor arm"),
    ("C2_MIR_AND", "C1_TXN_AND", "coherent, AND gate: realisable miRNA coherent (C2) vs Alon C1"),
    ("C2_MIR_OR",  "C1_TXN_OR",  "coherent, OR gate: realisable miRNA coherent (C2) vs Alon C1"),
    ("COMP_I1_AND", "I1_MIR_AND", "composite I1 (adds miRNA -| TF) vs its non-composite core"),
    ("COMP_C2_AND", "C2_MIR_AND", "composite C2 (adds miRNA -| TF) vs its non-composite core"),
    ("COMP_I1_OR",  "I1_MIR_OR",  "composite I1 (adds miRNA -| TF) vs its non-composite core"),
    ("COMP_C2_OR",  "C2_MIR_OR",  "composite C2 (adds miRNA -| TF) vs its non-composite core"),
]
rows = []
for a, b, why in PAIRS:
    A = E[E.topology == a].set_index("param_set")
    B = E[E.topology == b].set_index("param_set")
    idx = A.index.intersection(B.index)
    for m in METRICS:
        x = A.loc[idx, m].astype(float).values
        y = B.loc[idx, m].astype(float).values
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 30:
            continue
        d = x[ok] - y[ok]
        try:
            p = wilcoxon(x[ok], y[ok], zero_method="zsplit").pvalue if np.any(d != 0) else 1.0
        except Exception:
            p = np.nan
        rows.append(dict(comparison=f"{a} vs {b}", meaning=why, metric=m, n_pairs=int(ok.sum()),
                         median_A=float(np.median(x[ok])), median_B=float(np.median(y[ok])),
                         median_paired_difference=float(np.median(d)),
                         pct_A_greater=100 * float(np.mean(d > 0)),
                         wilcoxon_p=p))
P = pd.DataFrame(rows)
P.to_csv(f"{RES}/dynamics_mir_vs_txn_paired.csv", index=False)
print(f"paired comparisons: {len(P)} rows")
print(P[P.metric.isin(["rel_T_ON", "overshoot_ratio", "filter_index", "n_eff_protein", "is_pulse"])]
      .to_string(index=False))

# --------------------------------------------------- (b) titration threshold --
# Steady state of  dy/dt = by - y - theta*m*y ,  dm/dt = bm - gm*m - theta*m*y
#   -> theta*y^2 + y*(gm + theta*bm - theta*by) - by*gm = 0
def y_ss(by, bm, gm, theta):
    by = np.asarray(by, float)
    if theta == 0:
        return by
    A = theta
    B = gm + theta * bm - theta * by
    Cc = -by * gm
    disc = B ** 2 - 4 * A * Cc
    return (-B + np.sqrt(np.maximum(disc, 0.0))) / (2 * A)

def y_txn(by, r, K, n):
    """Same target under a transcriptional repressor at level r: production is scaled."""
    return by * (K ** n / (K ** n + r ** n))

bys = np.logspace(-2, 1.5, 400)
rows = []
gm = 0.3
for theta in (0.0, 1.0, 5.0, 20.0, 100.0):
    for bm in (1.0,):
        y = y_ss(bys, bm, gm, theta)
        ls = np.gradient(np.log(np.maximum(y, 1e-300)), np.log(bys))   # local log-log slope
        for byi, yi, si in zip(bys, y, ls):
            rows.append(dict(arm="miRNA_titration", theta=theta, beta_miRNA=bm, gamma_miRNA=gm,
                             beta_target=byi, target_ss=yi, local_loglog_slope=si))
# transcriptional repressor comparator at matched repression strength
for rlev in (0.5, 1.0, 2.0):
    y = y_txn(bys, rlev, 1.0, 2.0)
    ls = np.gradient(np.log(np.maximum(y, 1e-300)), np.log(bys))
    for byi, yi, si in zip(bys, y, ls):
        rows.append(dict(arm="transcriptional_repressor", theta=np.nan, beta_miRNA=np.nan,
                         gamma_miRNA=np.nan, beta_target=byi, target_ss=yi,
                         local_loglog_slope=si, repressor_level=rlev))
TT = pd.DataFrame(rows)
TT.to_csv(f"{RES}/dynamics_titration_threshold.csv", index=False)
summ = (TT[TT.arm == "miRNA_titration"].groupby("theta")
        .local_loglog_slope.max().rename("max_local_loglog_slope").reset_index())
summ2 = (TT[TT.arm == "transcriptional_repressor"].groupby("repressor_level")
         .local_loglog_slope.max().rename("max_local_loglog_slope").reset_index())
print("\nmaximum local log-log slope of target vs its own transcription rate")
print("  miRNA titration arm:"); print(summ.to_string(index=False))
print("  transcriptional repressor arm:"); print(summ2.to_string(index=False))
print("  (slope 1 = linear; slope > 1 = a threshold. A transcriptional repressor "
      "rescales the curve but cannot create a threshold.)")
