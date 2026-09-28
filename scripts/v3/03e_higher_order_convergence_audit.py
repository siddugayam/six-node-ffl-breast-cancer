#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03e_higher_order_convergence_audit.py
=====================================
HOW MUCH OF THE HIGHER-ORDER "GAIN" IS JUST AN UNCONVERGED INTEGRATION?

03_higher_order_sweep.py computes the steady-state dose-response by continuation,
integrating T_DOSE = 20 tau at each of 11 input levels.  The well-posedness check
(dynamics_model_wellposedness_check.csv) showed that with the saturating mutual
stabilisation 7-10% of parameter sets at n >= 4 are still relaxing at 250 tau.
For those sets a 20-tau continuation step is not a steady state, and both the
effective Hill coefficient and the dynamic range would be biased.

This script re-derives the dose-response for a random subsample of parameter sets
per module using CERTIFIED settling (integrate in 25-tau chunks until
max|dZ/dt| / max(|Z|, 1e-3) < 1e-5, giving up at 500 tau) and reports:
  * the fraction of sets that fail to certify at each input level,
  * the change in n_eff and in the dynamic range,
  * how many sets change their ultrasensitivity call (n_eff > 2).

Out: results/v3/dynamics_higher_order_convergence_audit.csv
"""
import sys, os, time, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
N_SET = 16384
SEED = 20260908
SUB = 2048
DT = 0.01
S_DOSE = np.logspace(np.log10(0.02), np.log10(20.0), 11)
TOL, CHUNK, TMAX = 1e-5, 25.0, 500.0
NPROC = 5


def eff_hill(S_grid, R):
    n = R.shape[1]; out = np.full(n, np.nan); logS = np.log(S_grid)
    lo = R.min(axis=0); hi = R.max(axis=0)
    ok = (hi - lo) > 1e-6
    for j in np.where(ok)[0]:
        r = R[:, j]
        i0, i1 = int(np.argmin(r)), int(np.argmax(r))
        a, b = (i0, i1) if i0 < i1 else (i1, i0)
        if b - a < 2:
            continue
        seg = r[a:b + 1]; sl = logS[a:b + 1]
        if seg[-1] < seg[0]:
            seg = seg[::-1]; sl = sl[::-1]
        t10 = seg[0] + 0.1 * (seg[-1] - seg[0]); t90 = seg[0] + 0.9 * (seg[-1] - seg[0])
        s10 = np.interp(t10, seg, sl); s90 = np.interp(t90, seg, sl)
        d = abs(s90 - s10)
        if d > 1e-9:
            out[j] = np.log(81.0) / d
    return out


def settle_certified(topo, P, Z0, S, act):
    N = Z0.shape[1]
    Zf = Z0.copy(); res = np.full(N, np.inf)
    idx = np.arange(N); Z = Z0.copy(); t = 0.0
    while len(idx) and t < TMAX:
        Ps = {k: (v[idx] if np.ndim(v) else v) for k, v in P.items()}
        rhs = D.make_rhs(topo, Ps)
        Z = D.integrate(rhs, Z, lambda tt: S, CHUNK, DT, active=act)
        t += CHUNK
        d = rhs(0.0, Z, S)[act]
        r = (np.abs(d) / np.maximum(np.abs(Z[act]), 1e-3)).max(axis=0)
        Zf[:, idx] = Z; res[idx] = r
        keep = r >= TOL
        idx = idx[keep]; Z = Z[:, keep]
    return Zf, res


def run(args):
    fam, size, topo = args
    t0 = time.time()
    P0 = SMP.sample(N_SET, seed=SEED)
    rng = np.random.default_rng(99)
    sub = np.sort(rng.choice(N_SET, SUB, replace=False))
    P = {k: np.asarray(v, float)[sub] for k, v in P0.items()}
    act = topo.active_states()
    N = SUB
    Z = np.zeros((D.NSTATE, N))
    Ps, unc = [], []
    for S in S_DOSE:
        Z, res = settle_certified(topo, P, Z, S, act)
        Ps.append(Z[D.IP1].copy())
        unc.append(res >= TOL)
    Ps = np.array(Ps); unc = np.array(unc)
    neff = eff_hill(S_DOSE, Ps)
    hi = Ps.max(axis=0); lo = Ps.min(axis=0)
    dr_fl = (hi + 1e-3) / (lo + 1e-3)
    out = pd.DataFrame(dict(family=fam, module=size, param_set=sub,
                            n_eff_certified=neff, dyn_range_floored_certified=dr_fl,
                            n_input_levels_uncertified=unc.sum(axis=0)))
    print(f"  {fam}/{size} done ({time.time()-t0:.0f}s)", flush=True)
    return out


if __name__ == "__main__":
    t0 = time.time()
    fams = {
        "I1_miRNA_FFL": D.higher_order_registry("AND", "I1_MIR", composite=False),
        "COMP_C2_toggle": D.higher_order_registry("AND", "C2_MIR", composite=True),
    }
    tasks = [(f, s, reg[s]) for f, reg in fams.items() for s in ("n3", "n4", "n5", "n6", "n4d")]
    with Pool(processes=NPROC) as pool:
        parts = pool.map(run, tasks)
    C = pd.concat(parts, ignore_index=True)

    A = pd.read_csv(f"{RES}/dynamics_higher_order_persetset.csv.gz",
                    usecols=["family", "module", "param_set", "n_eff", "dyn_range_floored"])
    M = C.merge(A, on=["family", "module", "param_set"], how="left")
    rows = []
    for (f, s), g in M.groupby(["family", "module"]):
        both = np.isfinite(g.n_eff) & np.isfinite(g.n_eff_certified)
        rows.append(dict(
            family=f, module=s, n_subsample=len(g),
            pct_sets_with_any_uncertified_input_level=100 * float((g.n_input_levels_uncertified > 0).mean()),
            median_n_input_levels_uncertified=float(g.n_input_levels_uncertified.median()),
            n_eff_median_20tau=float(np.nanmedian(g.n_eff)),
            n_eff_median_certified=float(np.nanmedian(g.n_eff_certified)),
            pct_ultrasensitive_20tau=100 * float(np.nanmean(g.n_eff > 2)),
            pct_ultrasensitive_certified=100 * float(np.nanmean(g.n_eff_certified > 2)),
            pct_ultrasensitivity_call_changed=100 * float(np.nanmean(
                (g.n_eff[both] > 2) != (g.n_eff_certified[both] > 2))) if both.any() else np.nan,
            dyn_range_floored_median_20tau=float(np.nanmedian(g.dyn_range_floored)),
            dyn_range_floored_median_certified=float(np.nanmedian(g.dyn_range_floored_certified)),
            median_abs_log2_ratio_dyn_range=float(np.nanmedian(np.abs(np.log2(
                (g.dyn_range_floored_certified + 1e-12) / (g.dyn_range_floored + 1e-12))))),
        ))
    S = pd.DataFrame(rows)
    S.to_csv(f"{RES}/dynamics_higher_order_convergence_audit.csv", index=False)
    pd.set_option("display.width", 260)
    print()
    print(S.round(3).to_string(index=False))
    print(f"\nDONE in {(time.time()-t0)/60:.1f} min")
