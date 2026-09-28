#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02f_pulse_longwindow.py
=======================
FIX FOR A REPORTING ARTEFACT IN THE PULSE STATISTICS.

The main battery (dyn_metrics.py) observes the ON step for 60 tau.  A pulse
whose excursion has not come back down inside that window has its full width at
half maximum CENSORED at the window edge, and 31-100% of the pulses detected per
topology were censored in that way.  A median "pulse width" of ~57 tau in a
60-tau window is a statement about the window, not about the circuit.

This script re-integrates the ON step for T_LONG = 400 tau (dt = 0.02, which the
integrator validation showed is accurate to ~1e-8 against LSODA) for every
parameter set that the 60-tau screen flagged as pulsing, and recomputes

    overshoot ratio, time to peak, true FWHM, adaptation precision

with the censoring rate at 400 tau reported explicitly.  It also re-tests a
random subsample of NON-flagged sets over the long window, to estimate how many
pulses the short window missed (a false-negative rate).

Out: results/v3/dynamics_pulse_longwindow.csv          (per parameter set)
     results/v3/dynamics_pulse_longwindow_summary.csv  (per topology)
"""
import sys, os, time, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_metrics as M

RES = "/path/to/revision/results/v3"
T_LONG = 400.0
DT = 0.02
REC = 25                      # record every 0.5 tau
N_NEG_SUBSAMPLE = 512
NPROC = 6

E = pd.read_csv(f"{RES}/dynamics_3node_metrics_ensemble.csv")
E = E[E.param_set >= 0]
Z = np.load(f"{RES}/dynamics_3node_sobol_params.npz")
Pall = {k: Z[k] for k in Z.files}     # column 0 is the reference set
reg = D.build_registry()
rng = np.random.default_rng(3)


def pulse_stats(ts, p, p_pre, p_ss):
    """Vectorised overshoot / t_peak / FWHM / adaptation over a long window."""
    N = p.shape[1]
    dyn = p - p_pre[None, :]
    fin = p_ss - p_pre
    sg = np.sign(dyn[np.abs(dyn).argmax(axis=0), np.arange(N)])
    sg[sg == 0] = 1.0
    ext = (dyn * sg[None, :]).max(axis=0)
    i_pk = (dyn * sg[None, :]).argmax(axis=0)
    peak = ext * sg
    with np.errstate(divide="ignore", invalid="ignore"):
        ovs = np.where(np.abs(fin) > 1e-9, peak / fin, np.nan)
        adapt = 1.0 - np.abs(fin) / np.maximum(np.abs(peak), 1e-12)
    half = p_pre + 0.5 * peak
    up = ((p - half[None, :]) * sg[None, :]) >= 0
    dt_rec = ts[1] - ts[0]
    n_t = p.shape[0]
    fwhm = np.full(N, np.nan); cens = np.zeros(N)
    for j in range(N):
        k = int(i_pk[j])
        if not up[k, j]:
            continue
        l = k
        while l > 0 and up[l - 1, j]:
            l -= 1
        r = k
        while r < n_t - 1 and up[r + 1, j]:
            r += 1
        fwhm[j] = (r - l) * dt_rec
        cens[j] = 1.0 if r == n_t - 1 else 0.0
    is_pulse = ((np.abs(ext) > 1e-6) & (np.abs(ext) > 1.10 * np.abs(fin))).astype(float)
    return dict(overshoot_ratio_long=ovs, t_peak_long=ts[i_pk], pulse_fwhm_long=fwhm,
                fwhm_censored_at_400tau=cens, adaptation_long=adapt,
                is_pulse_long=is_pulse, peak_excursion_long=peak, final_change_long=fin)


def run_one(nm):
    t0 = time.time()
    topo = reg[nm]
    g = E[E.topology == nm]
    pos = g.loc[g.is_pulse > 0, "param_set"].to_numpy(int)
    neg = g.loc[g.is_pulse == 0, "param_set"].to_numpy(int)
    if len(neg) > N_NEG_SUBSAMPLE:
        neg = np.sort(rng.choice(neg, N_NEG_SUBSAMPLE, replace=False))
    sets = np.concatenate([pos, neg]).astype(int)
    flag = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    if len(sets) == 0:
        return pd.DataFrame()
    cols = sets + 1                     # column 0 of Pall is the reference set
    P = {k: v[cols] for k, v in Pall.items()}
    N = len(sets)
    act = topo.active_states()
    rhs = D.make_rhs(topo, P)
    Zlo = D.integrate(rhs, np.zeros((D.NSTATE, N)), lambda t: M.S_LO, 80.0, DT, active=act)
    p_pre = Zlo[D.IP1].copy()
    Zon, ts, tr = D.integrate(rhs, Zlo, lambda t: M.S_HI, T_LONG, DT, record_every=REC,
                              active=act, record_states=(D.IP1,))
    p = tr[:, 0, :]
    st = pulse_stats(ts, p, p_pre, Zon[D.IP1].copy())
    out = pd.DataFrame(dict(topology=nm, param_set=sets,
                            flagged_pulse_in_60tau_window=flag, **st))
    print(f"  {nm:16s} n={N:5d}  ({time.time()-t0:6.1f}s)", flush=True)
    return out


if __name__ == "__main__":
    t0 = time.time()
    names = [n for n in reg if (E[E.topology == n].is_pulse > 0).any()
             or n in ("CASCADE", "C1_TXN_AND", "C1_TXN_OR")]
    with Pool(processes=NPROC) as pool:
        parts = pool.map(run_one, names)
    A = pd.concat([p for p in parts if len(p)], ignore_index=True)
    A.to_csv(f"{RES}/dynamics_pulse_longwindow.csv", index=False)

    rows = []
    for nm, g in A.groupby("topology", sort=False):
        pos = g[g.flagged_pulse_in_60tau_window > 0]
        neg = g[g.flagged_pulse_in_60tau_window == 0]
        unc = pos[pos.fwhm_censored_at_400tau == 0]
        rows.append(dict(
            topology=nm,
            n_flagged_pulse_60tau=len(pos),
            pct_of_flagged_still_pulsing_at_400tau=100 * float(np.nanmean(pos.is_pulse_long)) if len(pos) else np.nan,
            pct_of_flagged_censored_at_400tau=100 * float(np.nanmean(pos.fwhm_censored_at_400tau)) if len(pos) else np.nan,
            fwhm_median_tau_uncensored=float(np.nanmedian(unc.pulse_fwhm_long)) if len(unc) else np.nan,
            fwhm_q25_tau_uncensored=float(np.nanpercentile(unc.pulse_fwhm_long, 25)) if len(unc) else np.nan,
            fwhm_q75_tau_uncensored=float(np.nanpercentile(unc.pulse_fwhm_long, 75)) if len(unc) else np.nan,
            n_uncensored=len(unc),
            overshoot_median_long=float(np.nanmedian(pos.overshoot_ratio_long)) if len(pos) else np.nan,
            t_peak_median_tau_long=float(np.nanmedian(pos.t_peak_long)) if len(pos) else np.nan,
            adaptation_median_long=float(np.nanmedian(pos.adaptation_long)) if len(pos) else np.nan,
            n_nonflagged_subsample=len(neg),
            pct_nonflagged_that_pulse_at_400tau=100 * float(np.nanmean(neg.is_pulse_long)) if len(neg) else np.nan,
        ))
    S = pd.DataFrame(rows)
    S.to_csv(f"{RES}/dynamics_pulse_longwindow_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print()
    print(S.round(3).to_string(index=False))
    print(f"\nDONE in {time.time()-t0:.1f}s")
