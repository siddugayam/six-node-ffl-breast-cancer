#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06b_composite_showcase.py
-------------------------
Concrete, certified examples of the two behaviours that the COMPOSITE FFL has and no
acyclic (transcriptional or non-composite miRNA) FFL can have:

  (A) HYSTERESIS / BISTABILITY, from the double-negative TF -| miRNA -| TF arm
      (451 of the 1,434 composite cores; 292 of them sign-resolved).
  (B) RINGING / DAMPED OSCILLATION, from the TF -> miRNA -| TF delayed negative
      feedback arm (791 composite cores; 473 sign-resolved).

For each we search the Sobol ensemble for parameter sets that show the behaviour,
verify that both branches are at certified steady states, and export the curves.
The SAME parameter vector is then run through the matched acyclic circuit
(non-composite core and Alon transcriptional FFL) to show the behaviour disappears.

Out: dynamics_composite_hysteresis_curves.csv
     dynamics_composite_ringing_timecourses.csv
     dynamics_composite_showcase_summary.csv
"""
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
DT = 0.01
N = 4096
SEED = 20260908
reg = D.build_registry()
P = {k: np.asarray(v, float) for k, v in SMP.sample(16384, seed=SEED).items()}
P = {k: v[:N] for k, v in P.items()}


def settle(topo, Pp, Z0, S, act, tol=1e-7, chunk=25.0, tmax=400.0):
    n = Z0.shape[1]
    Zf = Z0.copy(); res = np.full(n, np.inf)
    idx = np.arange(n); Z = Z0.copy(); t = 0.0
    while len(idx) and t < tmax:
        Ps = {k: (v[idx] if np.ndim(v) else v) for k, v in Pp.items()}
        rhs = D.make_rhs(topo, Ps)
        Z = D.integrate(rhs, Z, lambda tt: S, chunk, DT, active=act)
        t += chunk
        d = rhs(0.0, Z, S)[act]
        sc = np.maximum(np.abs(Z[act]), 1e-3)
        r = (np.abs(d) / sc).max(axis=0)
        Zf[:, idx] = Z; res[idx] = r
        keep = r >= tol
        idx = idx[keep]; Z = Z[:, keep]
    return Zf, res


# ============================================================ (A) hysteresis ==
S_SCAN = np.logspace(np.log10(0.02), np.log10(20.0), 31)
topo = reg["COMP_C2_AND"]
act = topo.active_states()
t0 = time.time()
lowZ = np.zeros((D.NSTATE, N))
highZ = np.zeros((D.NSTATE, N)); highZ[act] = 3.0
gap = np.zeros(N); conv = np.ones(N, bool); best_S = np.zeros(N)
for S in S_SCAN[::3]:
    Zl, rl = settle(topo, P, lowZ.copy(), S, act)
    Zh, rh = settle(topo, P, highZ.copy(), S, act)
    den = np.maximum(np.maximum(np.abs(Zl[D.IP1]), np.abs(Zh[D.IP1])), 1e-6)
    g = np.abs(Zl[D.IP1] - Zh[D.IP1]) / den
    ok = (rl < 1e-7) & (rh < 1e-7)
    upd = (g > gap) & ok
    best_S = np.where(upd, S, best_S); gap = np.where(upd, g, gap)
    conv &= ok
print(f"  hysteresis scan {time.time()-t0:.0f}s; certified-bistable sets: "
      f"{int(((gap > 0.05) & conv).sum())} / {N}", flush=True)

cand = np.where((gap > 0.20) & conv & (P["ng12"] > 0))[0]
if len(cand) == 0:
    cand = np.where((gap > 0.05) & conv)[0]
# pick the median-hysteresis certified example (not the most extreme one)
pick = cand[np.argsort(gap[cand])[len(cand) // 2]] if len(cand) else None
rows_h = []
summ = []
if pick is not None:
    for nm in ("COMP_C2_AND", "C2_MIR_AND", "C1_TXN_AND"):
        tp = reg[nm]; a = tp.active_states()
        P1 = {k: v[[pick]] for k, v in P.items()}
        Zu = np.zeros((D.NSTATE, 1)); Zd = np.zeros((D.NSTATE, 1)); Zd[a] = 3.0
        up, dn = [], []
        for S in S_SCAN:                       # up sweep (continuation)
            Zu, _ = settle(tp, P1, Zu, S, a)
            up.append(float(Zu[D.IP1, 0]))
        for S in S_SCAN[::-1]:                 # down sweep
            Zd, _ = settle(tp, P1, Zd, S, a)
            dn.append(float(Zd[D.IP1, 0]))
        dn = dn[::-1]
        for S, u, d_ in zip(S_SCAN, up, dn):
            rows_h.append(dict(topology=nm, param_set=int(pick), S=S,
                               protein_up_sweep=u, protein_down_sweep=d_))
        den = np.maximum(np.maximum(np.abs(np.array(up)), np.abs(np.array(dn))), 1e-9)
        summ.append(dict(panel="A_hysteresis", topology=nm, param_set=int(pick),
                         has_reciprocal_arm=tp.mir_to_TF,
                         max_relative_gap_between_branches=float(
                             (np.abs(np.array(up) - np.array(dn)) / den).max()),
                         note="same parameter vector in every row"))
pd.DataFrame(rows_h).to_csv(f"{RES}/dynamics_composite_hysteresis_curves.csv", index=False)

# ============================================================ (B) ringing =====
topo = reg["COMP_I1_AND"]
act = topo.active_states()
Z, res = settle(topo, P, np.zeros((D.NSTATE, N)), 2.0, act)
rhs = D.make_rhs(topo, P)
k = len(act)
f0 = rhs(0.0, Z, 2.0)[act]
J = np.empty((N, k, k))
for j, idx in enumerate(act):
    Zp = Z.copy(); dh = 1e-6 * np.maximum(np.abs(Z[idx]), 1.0)
    Zp[idx] = Zp[idx] + dh
    J[:, :, j] = ((rhs(0.0, Zp, 2.0)[act] - f0) / dh[None, :]).T
ev = np.linalg.eigvals(J)
im = np.abs(ev.imag).max(axis=1)
re_at_max_im = np.take_along_axis(ev.real, np.argmax(np.abs(ev.imag), axis=1)[:, None], 1)[:, 0]
cand = np.where((im > 1e-3) & (res < 1e-7))[0]
print(f"  complex eigenvalues in COMP_I1_AND: {len(cand)} / {N} "
      f"({100*len(cand)/N:.1f}%)", flush=True)

rows_r = []
if len(cand):
    # the most under-damped certified example: largest |Im| / |Re|
    q = im[cand] / np.maximum(np.abs(re_at_max_im[cand]), 1e-12)
    pick2 = cand[int(np.argmax(q))]
    for nm in ("COMP_I1_AND", "I1_MIR_AND", "I1_TXN_AND", "CASCADE"):
        tp = reg[nm]; a = tp.active_states()
        P1 = {kk: v[[pick2]] for kk, v in P.items()}
        Zlo, _ = settle(tp, P1, np.zeros((D.NSTATE, 1)), 0.10, a)
        rr = D.make_rhs(tp, P1)
        _, ts, tr = D.integrate(rr, Zlo, lambda t: 2.0, 60.0, DT, record_every=5,
                                active=a, record_states=(D.IP1,))
        ss = float(tr[-1, 0, 0])
        for tt, v in zip(ts, tr[:, 0, 0]):
            rows_r.append(dict(topology=nm, param_set=int(pick2), t=tt,
                               protein=float(v), protein_over_steady_state=float(v) / max(ss, 1e-12)))
        summ.append(dict(panel="B_ringing", topology=nm, param_set=int(pick2),
                         has_reciprocal_arm=tp.mir_to_TF,
                         max_relative_gap_between_branches=np.nan,
                         note=f"|Im|/|Re| of the leading complex pair in COMP_I1_AND = "
                              f"{q.max():.2f}; period = {2*np.pi/im[pick2]:.1f} tau"))
pd.DataFrame(rows_r).to_csv(f"{RES}/dynamics_composite_ringing_timecourses.csv", index=False)
pd.DataFrame(summ).to_csv(f"{RES}/dynamics_composite_showcase_summary.csv", index=False)
print("\n", pd.DataFrame(summ).to_string(index=False))
print("\nwritten to", RES)
