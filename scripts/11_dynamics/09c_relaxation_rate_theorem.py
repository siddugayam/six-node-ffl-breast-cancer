#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09c_relaxation_rate_theorem.py
==============================
THE EXACT STATEMENT THAT SEPARATES A miRNA FFL FROM A TRANSCRIPTIONAL FFL.

In Alon's transcriptional FFL both arms enter the PRODUCTION term of the target.
The target's relaxation rate is then its own decay rate, whatever the arms do:

    dy/dt = f(X, Y) - alpha_y y      =>   d(dy/dt)/dy = -alpha_y   (always)

so the response half-time is ln2 / alpha_y = ln2 in scaled time, for every
parameter set, every input function and every loop sign.

A miRNA arm acts on the transcript, so it enters the DEGRADATION term:

    dy/dt = f(X) - (1 + lambda_y q(M_eff)) y - theta M_eff y

    =>   d(dy/dt)/dy = -(1 + lambda_y q(M_eff) + theta M_eff)

and the half-time is ln2 / (1 + lambda_y q + theta M_eff) <= ln2, with equality
only when the miRNA is absent.  The miRNA arm therefore ACCELERATES the target
transcript's relaxation for ANY loop sign - coherent or incoherent - which no
promoter-bound repressor can do.  Acceleration is not a property of incoherence
here, as it is in the transcriptional I1-FFL; it is a property of the layer the
repressor acts on.

This script verifies the identity numerically to machine precision on a random
parameter sample, for miRNA-arm and transcriptional-arm circuits alike.

Out: results/v3/dynamics_relaxation_rate_theorem.csv
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
N = 1024
DT = 0.005
S_OP = 2.0
rng = np.random.default_rng(20260909)

P = {k: np.asarray(v, float) for k, v in SMP.sample(N, seed=4242).items()}
# gm is bounded above so the explicit RK4 stays stable at dt = 0.005
P["gm"] = np.clip(P["gm"], 0.1, 1.5)
reg = D.build_registry()

rows = []
for nm, topo in reg.items():
    act = topo.active_states()
    rhs = D.make_rhs(topo, P)
    Z = D.integrate(rhs, np.zeros((D.NSTATE, N)), lambda t: S_OP, 80.0, DT, active=act)
    m = Z[D.IM1]
    Meff = m / (1.0 + m / P["Kc"])
    h = 1e-7
    Zp = Z.copy(); Zp[D.IY1] = Zp[D.IY1] + h
    rate_obs = -(rhs(0.0, Zp, S_OP)[D.IY1] - rhs(0.0, Z, S_OP)[D.IY1]) / h
    if topo.mir_arm:
        q = 1.0 - D.hill_rep(Meff, P["Kmy"], P["nmy"])
        rate_pred = 1.0 + P["lam_y"] * q + P["theta"] * Meff
    else:
        rate_pred = np.ones(N)
    d = np.abs(rate_obs - rate_pred)
    th_obs = np.log(2.0) / np.maximum(rate_obs, 1e-12)
    rows.append(dict(
        topology=nm, miRNA_arm_post_transcriptional=topo.mir_arm, gate=topo.gate,
        census_class=topo.core_class,
        n_param_sets=N,
        predicted_relaxation_rate_median=float(np.median(rate_pred)),
        observed_relaxation_rate_median=float(np.median(rate_obs)),
        max_abs_deviation_from_identity=float(d.max()),
        median_abs_deviation_from_identity=float(np.median(d)),
        mRNA_half_time_median_tau=float(np.median(th_obs)),
        mRNA_half_time_p05_tau=float(np.percentile(th_obs, 5)),
        mRNA_half_time_p95_tau=float(np.percentile(th_obs, 95)),
        # strict "<" against ln2 would count last-bit rounding for the transcriptional
        # circuits, whose rate is exactly 1; a relative tolerance of 1e-9 is used
        pct_sets_faster_than_ln2=100 * float(np.mean(th_obs < np.log(2) * (1 - 1e-9))),
    ))
T = pd.DataFrame(rows)
T.to_csv(f"{RES}/dynamics_relaxation_rate_theorem.csv", index=False)
pd.set_option("display.width", 250)
print(T.round(6).to_string(index=False))
print()
print("largest deviation from the identity anywhere:",
      f"{T.max_abs_deviation_from_identity.max():.3e}")
mir = T[T.miRNA_arm_post_transcriptional]
txn = T[~T.miRNA_arm_post_transcriptional]
print(f"miRNA-arm circuits: mRNA half-time median {mir.mRNA_half_time_median_tau.median():.4f} tau; "
      f"{mir.pct_sets_faster_than_ln2.min():.1f}-{mir.pct_sets_faster_than_ln2.max():.1f}% of "
      f"parameter sets faster than ln2")
print(f"transcriptional-arm circuits: mRNA half-time exactly ln2 = {np.log(2):.4f} tau in "
      f"{100 - txn.pct_sets_faster_than_ln2.max():.1f}% of sets")
