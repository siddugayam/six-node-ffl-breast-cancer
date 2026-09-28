#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dyn_metrics.py
Dynamical characterisation of the FFL topologies defined in dyn_models.py.

All metrics follow Alon's definitions where one exists:
  * response time  = time to reach half of the FINAL steady-state change
                     (Alon 2006 doi:10.1201/9781420011432, ch. 2;
                      Mangan & Alon 2003 doi:10.1073/pnas.2133841100)
  * sign-sensitive delay = delay present on one direction of the step only
                     (Mangan, Zaslaver & Alon 2003 doi:10.1016/j.jmb.2003.09.049)
  * pulse / adaptation    = I1-FFL signature (Mangan & Alon 2003;
                            Bleris et al. 2011 doi:10.1038/msb.2011.49)
  * fold-change detection = invariance of the whole trajectory to the absolute
                            input level at fixed fold change
                            (Goentoro et al. 2009 doi:10.1016/j.molcel.2009.11.018)
  * persistence detection / noise filtering = the classical FFL function
                            (Mangan & Alon 2003; Alon 2007 doi:10.1038/nrg2102)
Time unit throughout = the target-mRNA lifetime tau = 1/alpha_Y.
"""
from __future__ import annotations
import numpy as np
import dyn_models as D

DT       = 0.01
DR_FLOOR = 1e-3              # detection floor for the steady-state dynamic range
T_SETTLE = 80.0
T_OBS    = 60.0
REC      = 5                 # record every 5 steps -> resolution 0.05 tau
S_LO, S_HI = 0.10, 2.00


# ---------------------------------------------------------------- helpers ---
def _first_crossing(ts, Y, level, rising):
    """First time each column of Y (n_t, N) crosses `level` (N,). NaN if never."""
    if rising:
        hit = Y >= level[None, :]
    else:
        hit = Y <= level[None, :]
    idx = np.argmax(hit, axis=0)
    ok = hit.any(axis=0)
    out = np.full(Y.shape[1], np.nan)
    i = idx[ok]
    cols = np.where(ok)[0]
    t1 = ts[i]
    y1 = Y[i, cols]
    prev = np.maximum(i - 1, 0)
    t0 = ts[prev]
    y0 = Y[prev, cols]
    denom = (y1 - y0)
    frac = np.where(np.abs(denom) > 1e-15, (level[cols] - y0) / np.where(np.abs(denom) > 1e-15, denom, 1.0), 0.0)
    frac = np.clip(frac, 0.0, 1.0)
    out[cols] = t0 + frac * (t1 - t0)
    return out


def response_half_time(ts, Y, y_pre, y_ss):
    """Time to reach y_pre + 0.5*(y_ss - y_pre)."""
    lvl = y_pre + 0.5 * (y_ss - y_pre)
    rising = (y_ss >= y_pre)
    t = np.full(Y.shape[1], np.nan)
    if rising.any():
        cols = np.where(rising)[0]
        t[cols] = _first_crossing(ts, Y[:, cols], lvl[cols], True)
    if (~rising).any():
        cols = np.where(~rising)[0]
        t[cols] = _first_crossing(ts, Y[:, cols], lvl[cols], False)
    return t


# ---------------------------------------------------------------- protocols -
def _rhsact(topo, P):
    return D.make_rhs(topo, P), topo.active_states()


def equilibrate(topo, P, S, N, Z0=None, t=T_SETTLE):
    rhs, act = _rhsact(topo, P)
    if Z0 is None:
        Z0 = np.zeros((D.NSTATE, N))
    return D.integrate(rhs, Z0, lambda tt: S, t, DT, active=act)


def step_response(topo, P, N, Zpre, S_new, t_obs=T_OBS, rec=REC, states=(D.IY1, D.IP1)):
    rhs, act = _rhsact(topo, P)
    Zf, ts, tr = D.integrate(rhs, Zpre, lambda tt: S_new, t_obs, DT, record_every=rec,
                             active=act, record_states=states)
    return Zf, ts, tr


def dose_response(topo, P, N, S_grid, t_each=40.0):
    """Continuation up the S grid; returns (n_S, N) arrays of y1 and p1 steady states."""
    rhs, act = _rhsact(topo, P)
    Z = np.zeros((D.NSTATE, N))
    ys, ps, xs, ms = [], [], [], []
    for S in S_grid:
        Z = D.integrate(rhs, Z, lambda tt: S, t_each, DT, active=act)
        ys.append(Z[D.IY1].copy()); ps.append(Z[D.IP1].copy())
        xs.append(Z[D.IX1].copy()); ms.append(Z[D.IM1].copy())
    return np.array(ys), np.array(ps), np.array(xs), np.array(ms)


# ------------------------------------------------------------ metric blocks -
def eff_hill(S_grid, R):
    """Effective Hill coefficient n_eff = ln(81)/ln(S90/S10) of a monotone dose-response."""
    N = R.shape[1]
    out = np.full(N, np.nan)
    lo = R.min(axis=0); hi = R.max(axis=0)
    rng = hi - lo
    ok = rng > 1e-6
    logS = np.log(S_grid)
    for j in np.where(ok)[0]:
        r = R[:, j]
        # use the monotone segment between the global min and global max
        i0, i1 = int(np.argmin(r)), int(np.argmax(r))
        a, b = (i0, i1) if i0 < i1 else (i1, i0)
        if b - a < 2:
            continue
        seg = r[a:b+1]; sl = logS[a:b+1]
        if seg[-1] < seg[0]:
            seg = seg[::-1]; sl = sl[::-1]
        t10 = seg[0] + 0.1 * (seg[-1] - seg[0])
        t90 = seg[0] + 0.9 * (seg[-1] - seg[0])
        try:
            s10 = np.interp(t10, seg, sl); s90 = np.interp(t90, seg, sl)
        except Exception:
            continue
        d = abs(s90 - s10)
        if d > 1e-9:
            out[j] = np.log(81.0) / d
    return out


def characterise(topo, P, N, S_grid=None, t_obs=T_OBS):
    """Full metric battery for one topology across N parameter sets (vectorised)."""
    M = {}
    states = (D.IY1, D.IP1)
    iY, iP = 0, 1

    # ---- ON step  S_LO -> S_HI ------------------------------------------
    Zlo = equilibrate(topo, P, S_LO, N)
    p_pre = Zlo[D.IP1].copy()
    Zon, ts, tr_on = step_response(topo, P, N, Zlo, S_HI, t_obs, states=states)
    p_on = tr_on[:, iP, :]
    p_ss_on = Zon[D.IP1].copy()
    M["p_lo"] = p_pre
    M["p_hi"] = p_ss_on
    M["T_half_ON"] = response_half_time(ts, p_on, p_pre, p_ss_on)

    # pulse / overshoot on the ON step
    dyn = p_on - p_pre[None, :]
    fin = p_ss_on - p_pre
    # peak in the direction of the eventual change (or the dominant excursion)
    sgn = np.where(np.abs(dyn).max(axis=0) > 0, np.sign(dyn[np.abs(dyn).argmax(axis=0),
                                                          np.arange(N)]), 1.0)
    sgn[sgn == 0] = 1.0
    ext = (dyn * sgn[None, :]).max(axis=0)                # signed peak excursion
    i_pk = (dyn * sgn[None, :]).argmax(axis=0)
    M["peak_excursion"] = ext * sgn
    M["t_peak"] = ts[i_pk]
    with np.errstate(divide="ignore", invalid="ignore"):
        M["overshoot_ratio"] = np.where(np.abs(fin) > 1e-9, (ext * sgn) / fin, np.nan)
        # adaptation precision: 1 = perfect adaptation (returns to baseline)
        M["adaptation"] = 1.0 - np.abs(fin) / np.maximum(np.abs(ext), 1e-12)
    M["is_pulse"] = ((np.abs(ext) > 1e-6) &
                     (np.abs(ext) > 1.10 * np.abs(fin)) &
                     (ts[i_pk] < 0.95 * t_obs)).astype(float)
    # pulse width.  Two numbers, because they mean different things:
    #   pulse_width_above_half_peak  total time spent above the half-peak level
    #   pulse_fwhm                   width of the CONTIGUOUS excursion containing the
    #                                peak, i.e. the true full width at half maximum.
    # A trajectory that rises and never comes back down has a half-peak "width" equal
    # to the observation window; pulse_fwhm_censored flags exactly those, so a censored
    # width is never reported as if it were a measured pulse duration.
    half = p_pre + 0.5 * (ext * sgn)
    up = ((p_on - half[None, :]) * sgn[None, :]) >= 0
    dt_rec = ts[1] - ts[0]
    M["pulse_width_above_half_peak"] = np.where(M["is_pulse"] > 0,
                                                up.sum(axis=0) * dt_rec, np.nan)
    n_t = p_on.shape[0]
    fwhm = np.full(N, np.nan); cens = np.full(N, np.nan)
    for j in np.where(M["is_pulse"] > 0)[0]:
        u = up[:, j]; k = int(i_pk[j])
        if not u[k]:
            continue
        l = k
        while l > 0 and u[l - 1]:
            l -= 1
        r = k
        while r < n_t - 1 and u[r + 1]:
            r += 1
        fwhm[j] = (r - l) * dt_rec
        cens[j] = 1.0 if r == n_t - 1 else 0.0
    M["pulse_fwhm"] = fwhm
    M["pulse_fwhm_censored"] = cens

    # ---- OFF step  S_HI -> S_LO -----------------------------------------
    Zhi = equilibrate(topo, P, S_HI, N)
    p_pre2 = Zhi[D.IP1].copy()
    Zoff, ts2, tr_off = step_response(topo, P, N, Zhi, S_LO, t_obs, states=states)
    p_off = tr_off[:, iP, :]
    p_ss_off = Zoff[D.IP1].copy()
    M["T_half_OFF"] = response_half_time(ts2, p_off, p_pre2, p_ss_off)

    # ---- ringing / damped oscillation on the ON step ---------------------
    d = np.diff(p_on, axis=0)
    sgn_changes = (np.diff(np.sign(d + 0.0), axis=0) != 0).sum(axis=0)
    M["n_turning_points"] = sgn_changes.astype(float)

    # ---- fold-change detection ------------------------------------------
    F = 3.0
    bases = (0.25, 0.5, 1.0)
    curves = []
    for S0 in bases:
        Z0 = equilibrate(topo, P, S0, N)
        b0 = Z0[D.IP1].copy()
        _, tsf, trf = step_response(topo, P, N, Z0, F * S0, 30.0, states=states)
        c = trf[:, iP, :] / np.maximum(b0[None, :], 1e-9)
        curves.append(c)
    C = np.array(curves)                                   # (3, n_t, N)
    span = C.max(axis=0) - C.min(axis=0)
    scale = np.maximum(np.abs(C.mean(axis=0)), 1e-9)
    M["FCD_error"] = (span / scale).max(axis=0)
    # A circuit that simply does not respond to the input is trivially "invariant".
    # Fold-change detection is only scored where the fold input actually moves the
    # output by at least 10% of its pre-step level.
    M["FCD_response_amplitude"] = np.abs(C[:, -1, :].mean(axis=0) - 1.0)
    M["is_FCD"] = ((M["FCD_error"] < 0.10) &
                   (M["FCD_response_amplitude"] > 0.10)).astype(float)

    # ---- persistence detection / noise filtering -------------------------
    # brief spurious input (0.5 tau) vs persistent input, same amplitude
    Zlo2 = equilibrate(topo, P, S_LO, N)
    rhs, act = _rhsact(topo, P)
    d_short = 0.5
    S_short = lambda tt: (S_HI if tt < d_short else S_LO)
    _, ts3, tr_s = D.integrate(rhs, Zlo2, S_short, t_obs, DT, record_every=REC,
                               active=act, record_states=states)
    p_short = tr_s[:, iP, :]
    base = Zlo2[D.IP1]
    amp_short = np.abs(p_short - base[None, :]).max(axis=0)
    amp_long = np.abs(p_on - p_pre[None, :]).max(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        M["noise_transmission"] = np.where(amp_long > 1e-9, amp_short / amp_long, np.nan)
    M["filter_index"] = 1.0 - M["noise_transmission"]

    # ---- steady-state input-output ---------------------------------------
    if S_grid is None:
        S_grid = np.logspace(np.log10(0.02), np.log10(20.0), 21)
    Ys, Ps, Xs, Ms = dose_response(topo, P, N, S_grid)
    M["n_eff_protein"] = eff_hill(S_grid, Ps)
    M["n_eff_mRNA"] = eff_hill(S_grid, Ys)
    lo = Ps.min(axis=0); hi = Ps.max(axis=0)
    # The raw ratio diverges as the OFF state approaches zero and is then a statement
    # about floating point, not about biology; dyn_range_protein_floored adds an
    # explicit detection floor (0.1% of the maximal production level).
    M["dyn_range_protein"] = (hi + 1e-9) / (lo + 1e-9)
    M["dyn_range_protein_floored"] = (hi + DR_FLOOR) / (lo + DR_FLOOR)
    M["ss_max_protein"] = hi
    M["ss_min_protein"] = lo
    # non-monotonicity must exceed 2% of the response range, not machine epsilon
    dP = np.diff(Ps, axis=0)
    thr = 0.02 * np.maximum(hi - lo, 1e-12)
    M["nonmonotone_dose"] = ((dP > thr[None, :]).any(axis=0) &
                             (dP < -thr[None, :]).any(axis=0)).astype(float)
    return M, dict(S_grid=S_grid, Ps=Ps, Ys=Ys, ts=ts, p_on=p_on, ts2=ts2, p_off=p_off,
                   ts3=ts3, p_short=p_short)
