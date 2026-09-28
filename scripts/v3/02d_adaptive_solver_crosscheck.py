#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02d_adaptive_solver_crosscheck.py
=================================
INDEPENDENT NUMERICAL VERIFICATION of the dynamical battery.

01_integrator_validation.py compared the fixed-step vectorised RK4 against
scipy solve_ivp at STEADY STATE only.  The metrics that the paper reports
(response half-times, pulses, filter index) are properties of TRANSIENTS, so
this script re-derives them for a random subsample of the actual Sobol
parameter sets using scipy.integrate.solve_ivp (LSODA, rtol=1e-8, atol=1e-11,
adaptive step, dense output) and compares them, set by set, with the values in
dynamics_3node_metrics_ensemble.csv.

It also re-derives the steady-state dose-response and its effective Hill
coefficient with the adaptive solver.

Out: results/v3/dynamics_solver_crosscheck.csv
     results/v3/dynamics_solver_crosscheck_summary.csv
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_metrics as M
from scipy.integrate import solve_ivp

RES = "/path/to/revision/results/v3"
RTOL, ATOL = 1e-8, 1e-11
N_SETS = 40
TOPOS = ["CASCADE", "C1_TXN_AND", "C1_TXN_OR", "I1_TXN_AND", "I1_MIR_AND",
         "COMP_I1_AND", "COMP_C2_AND", "C3_MIR_AND"]

reg = D.build_registry()
Z = np.load(f"{RES}/dynamics_3node_sobol_params.npz")
Pall = {k: Z[k] for k in Z.files}
NT = len(Pall["gx"])
rng = np.random.default_rng(11)
# IMPORTANT INDEXING CONVENTION.  Column 0 of dynamics_3node_sobol_params.npz is the
# REFERENCE parameter set, which the ensemble file labels param_set = -1.  Ensemble
# parameter set j therefore lives in column j+1.  Getting this wrong silently compares
# one parameter vector's metrics with another's.
sets = np.sort(rng.choice(np.arange(0, NT - 1), size=N_SETS, replace=False))

E = pd.read_csv(f"{RES}/dynamics_3node_metrics_ensemble.csv")


def scalar_rhs(topo, Pj):
    rhs = D.make_rhs(topo, Pj)
    act = topo.active_states()

    def f(t, z, S):
        dZ = rhs(t, z.reshape(D.NSTATE, 1), S)
        out = np.zeros(D.NSTATE)
        out[act] = dZ[act, 0]
        return out
    return f, act


def settle(f, z0, S, t_end=200.0):
    s = solve_ivp(f, (0, t_end), z0, args=(S,), method="LSODA", rtol=RTOL, atol=ATOL)
    return s.y[:, -1]


def traj(f, z0, S_of_t, t_end, n=1201):
    ts = np.linspace(0, t_end, n)
    s = solve_ivp(lambda t, z: f(t, z, S_of_t(t)), (0, t_end), z0, method="LSODA",
                  rtol=RTOL, atol=ATOL, t_eval=ts, max_step=0.25)
    return s.t, s.y


def half_time_1d(ts, y, y_pre, y_ss):
    lvl = y_pre + 0.5 * (y_ss - y_pre)
    if y_ss >= y_pre:
        hit = np.where(y >= lvl)[0]
    else:
        hit = np.where(y <= lvl)[0]
    if len(hit) == 0:
        return np.nan
    i = hit[0]
    if i == 0:
        return 0.0
    t0, t1 = ts[i - 1], ts[i]
    y0, y1 = y[i - 1], y[i]
    if abs(y1 - y0) < 1e-15:
        return t1
    return t0 + (lvl - y0) / (y1 - y0) * (t1 - t0)


rows = []
S_grid = np.logspace(np.log10(0.02), np.log10(20.0), 41)
for nm in TOPOS:
    topo = reg[nm]
    Esub = E[E.topology == nm].set_index("param_set")
    for j in sets:
        Pj = {k: np.array([v[j + 1]]) for k, v in Pall.items()}
        f, act = scalar_rhs(topo, Pj)
        z0 = np.zeros(D.NSTATE)
        zlo = settle(f, z0, M.S_LO)
        p_pre = zlo[D.IP1]
        ts, Y = traj(f, zlo, lambda t: M.S_HI, M.T_OBS)
        p_on = Y[D.IP1]
        zon = settle(f, zlo, M.S_HI)
        p_ss = zon[D.IP1]
        T_on = half_time_1d(ts, p_on, p_pre, p_ss)

        zhi = settle(f, z0, M.S_HI)
        p_pre2 = zhi[D.IP1]
        ts2, Y2 = traj(f, zhi, lambda t: M.S_LO, M.T_OBS)
        zoff = settle(f, zhi, M.S_LO)
        T_off = half_time_1d(ts2, Y2[D.IP1], p_pre2, zoff[D.IP1])

        # brief pulse vs persistent
        ts3, Y3 = traj(f, zlo, lambda t: (M.S_HI if t < 0.5 else M.S_LO), M.T_OBS)
        amp_short = np.abs(Y3[D.IP1] - p_pre).max()
        amp_long = np.abs(p_on - p_pre).max()
        filt = 1.0 - (amp_short / amp_long if amp_long > 1e-9 else np.nan)

        # overshoot
        dyn = p_on - p_pre
        fin = p_ss - p_pre
        sg = np.sign(dyn[np.abs(dyn).argmax()]) or 1.0
        ext = (dyn * sg).max() * sg
        ovs = ext / fin if abs(fin) > 1e-9 else np.nan

        # dose-response (continuation, adaptive)
        zz = np.zeros(D.NSTATE)
        Ps = []
        for S in S_grid:
            zz = settle(f, zz, S, 120.0)
            Ps.append(zz[D.IP1])
        Ps = np.array(Ps)[:, None]
        neff = M.eff_hill(S_grid, Ps)[0]

        r = Esub.loc[j]
        rows.append(dict(topology=nm, param_set=int(j),
                         T_half_ON_rk4=r.T_half_ON, T_half_ON_lsoda=T_on,
                         T_half_OFF_rk4=r.T_half_OFF, T_half_OFF_lsoda=T_off,
                         filter_index_rk4=r.filter_index, filter_index_lsoda=filt,
                         overshoot_rk4=r.overshoot_ratio, overshoot_lsoda=ovs,
                         n_eff_rk4=r.n_eff_protein, n_eff_lsoda=neff,
                         p_lo_rk4=r.p_lo, p_lo_lsoda=p_pre,
                         p_hi_rk4=r.p_hi, p_hi_lsoda=p_ss))
    print(f"  {nm} done", flush=True)

R = pd.DataFrame(rows)
R.to_csv(f"{RES}/dynamics_solver_crosscheck.csv", index=False)

pairs = [("T_half_ON", "tau"), ("T_half_OFF", "tau"), ("filter_index", "-"),
         ("overshoot", "-"), ("n_eff", "-"), ("p_lo", "-"), ("p_hi", "-")]
sm = []
for base, unit in pairs:
    a = R[f"{base}_rk4"].to_numpy(float)
    b = R[f"{base}_lsoda"].to_numpy(float)
    ok = np.isfinite(a) & np.isfinite(b)
    den = np.maximum(np.abs(a[ok]), 1e-3)
    sm.append(dict(metric=base, unit=unit, n=int(ok.sum()),
                   n_finiteness_mismatch=int((np.isfinite(a) != np.isfinite(b)).sum()),
                   max_abs_diff=float(np.abs(a[ok] - b[ok]).max()),
                   median_abs_diff=float(np.median(np.abs(a[ok] - b[ok]))),
                   max_rel_diff=float((np.abs(a[ok] - b[ok]) / den).max()),
                   median_rel_diff=float(np.median(np.abs(a[ok] - b[ok]) / den)),
                   pearson_r=float(np.corrcoef(a[ok], b[ok])[0, 1])))
S = pd.DataFrame(sm)
S.to_csv(f"{RES}/dynamics_solver_crosscheck_summary.csv", index=False)
pd.set_option("display.width", 220)
print(S.to_string(index=False))
print("\nrows:", len(R), " topologies:", R.topology.nunique(), " param sets each:", N_SETS)
