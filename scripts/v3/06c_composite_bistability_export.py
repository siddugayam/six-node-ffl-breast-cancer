#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06c_composite_bistability_export.py
===================================
CORRECTION TO THE BISTABILITY FIGURE.

06b_composite_showcase.py certified bistability with a two-initial-condition test
(integrate from Z = 0 and from Z = 3 at the same input, both to a certified steady
state) but then exported CONTINUATION sweeps for the figure.  Those are not the
same protocol: continuation in S only reveals the two branches if the input can
carry the system across a saddle-node inside the scanned range.  For the composite
C2 toggle the two attractors are separated in the STATE variables, not along S, so
the continuation sweeps collapse onto one curve and the exported figure showed a
maximum branch gap of 0.0009 - i.e. no visible hysteresis - even though 1,219 of
4,096 parameter sets are certified bistable.

This script re-exports the curves with the protocol that actually certified the
behaviour: at each input level, both initial conditions are integrated to a
certified steady state independently.  The continuation sweeps are exported
alongside, so the difference between the two protocols is visible rather than
hidden.

The example is chosen from the 16,384-set higher-order sweep, which used the same
sampler and seed, so the parameter indices are directly comparable.

Out: results/v3/dynamics_composite_hysteresis_curves.csv  (overwrites)
     results/v3/dynamics_composite_bistability_example.csv
"""
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
DT = 0.01
SEED = 20260908
N_SET = 16384
TOL, CHUNK, TMAX = 1e-5, 25.0, 500.0
S_SCAN = np.logspace(np.log10(0.02), np.log10(20.0), 31)

A = pd.read_csv(f"{RES}/dynamics_higher_order_persetset.csv.gz",
                usecols=["family", "module", "param_set", "is_bistable", "hysteresis"])
g = A[(A.family == "COMP_C2_toggle") & (A.module == "n3") & (A.is_bistable > 0)]
print(f"certified-bistable COMP_C2 3-node cores in the 16,384-set sweep: {len(g)} "
      f"({100*len(g)/N_SET:.2f}%)")
g = g.sort_values("hysteresis", ascending=False)
top = g.head(max(1, len(g) // 10))
pick = int(top.iloc[len(top) // 2].param_set)
print(f"example parameter set: {pick} (relative branch gap in the sweep = "
      f"{float(g[g.param_set == pick].hysteresis.iloc[0]):.3f})")

P0 = {k: np.asarray(v, float) for k, v in SMP.sample(N_SET, seed=SEED).items()}
P = {k: v[[pick]] for k, v in P0.items()}
reg = D.build_registry()


def settle(topo, Pp, Z0, S, act):
    Z = Z0.copy(); t = 0.0; r = np.inf
    rhs = D.make_rhs(topo, Pp)
    while t < TMAX:
        Z = D.integrate(rhs, Z, lambda tt: S, CHUNK, DT, active=act)
        t += CHUNK
        d = rhs(0.0, Z, S)[act]
        r = float((np.abs(d) / np.maximum(np.abs(Z[act]), 1e-3)).max())
        if r < TOL:
            break
    return Z, r, t


rows = []
t0 = time.time()
for nm in ("COMP_C2_AND", "C2_MIR_AND", "C1_TXN_AND"):
    topo = reg[nm]; act = topo.active_states()
    Zu = np.zeros((D.NSTATE, 1))                       # continuation, up
    up = []
    for S in S_SCAN:
        Zu, _, _ = settle(topo, P, Zu, S, act)
        up.append(float(Zu[D.IP1, 0]))
    Zd = np.zeros((D.NSTATE, 1)); Zd[act] = 3.0        # continuation, down
    dn = []
    for S in S_SCAN[::-1]:
        Zd, _, _ = settle(topo, P, Zd, S, act)
        dn.append(float(Zd[D.IP1, 0]))
    dn = dn[::-1]
    lo, hi, rl, rh = [], [], [], []
    for S in S_SCAN:                                   # independent two-IC test
        Z0 = np.zeros((D.NSTATE, 1))
        Zl, res_l, _ = settle(topo, P, Z0, S, act)
        Zh0 = np.zeros((D.NSTATE, 1)); Zh0[act] = 3.0
        Zh, res_h, _ = settle(topo, P, Zh0, S, act)
        lo.append(float(Zl[D.IP1, 0])); hi.append(float(Zh[D.IP1, 0]))
        rl.append(res_l); rh.append(res_h)
    for S, u, d_, l_, h_, a_, b_ in zip(S_SCAN, up, dn, lo, hi, rl, rh):
        rows.append(dict(topology=nm, param_set=pick, S=S,
                         protein_low_IC=l_, protein_high_IC=h_,
                         protein_up_sweep=u, protein_down_sweep=d_,
                         residual_low_IC=a_, residual_high_IC=b_,
                         both_branches_certified=int((a_ < TOL) and (b_ < TOL))))
    print(f"  {nm:12s} done ({time.time()-t0:.0f}s)", flush=True)

H = pd.DataFrame(rows)
H.to_csv(f"{RES}/dynamics_composite_hysteresis_curves.csv", index=False)

summ = []
for nm, d in H.groupby("topology", sort=False):
    den = np.maximum(np.maximum(np.abs(d.protein_low_IC), np.abs(d.protein_high_IC)), 1e-9)
    gap_ic = float((np.abs(d.protein_low_IC - d.protein_high_IC) / den).max())
    den2 = np.maximum(np.maximum(np.abs(d.protein_up_sweep), np.abs(d.protein_down_sweep)), 1e-9)
    gap_cont = float((np.abs(d.protein_up_sweep - d.protein_down_sweep) / den2).max())
    summ.append(dict(topology=nm, param_set=pick,
                     has_reciprocal_miRNA_to_TF_arm=reg[nm].mir_to_TF,
                     max_relative_gap_two_initial_conditions=gap_ic,
                     max_relative_gap_continuation_sweep=gap_cont,
                     n_input_levels_with_both_branches_certified=int(d.both_branches_certified.sum()),
                     n_input_levels=len(d)))
S = pd.DataFrame(summ)
S.to_csv(f"{RES}/dynamics_composite_bistability_example.csv", index=False)
pd.set_option("display.width", 250)
print()
print(S.round(4).to_string(index=False))
