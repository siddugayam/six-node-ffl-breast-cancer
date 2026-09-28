#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_higher_order_sweep.py  --  Part C: DOES ORDER ADD FUNCTION?

For two circuit families (non-composite I1 miRNA-FFL, and composite C2 with the
reciprocal miRNA -| TF arm) we build the 3 -> 4 -> 5 -> 6 node escalation that the
manuscript defines:

   n3  core FFL                     TF -> miR -| G1 ; TF -> G1
   n4  + gene-gene edge             G1 <-> G2 (STRING edges are undirected, so the
                                    generous reading is mutual protein stabilisation;
                                    a directed G1 -> G2 cascade variant is also run)
   n5  + miRNA-miRNA coupling       miR2 co-transcribed with miR1 from the shared
                                    cluster promoter, both competing for one RISC pool
   n6  + TF-TF edge                 TF1 -> TF2, TF2 also regulating G1 and the miRNA

Every larger module CONTAINS the n = 3 core and is evaluated at the SAME parameter
vector, so a behaviour present at n = k and absent at n = 3 for that same vector is
an unambiguous gain of function from the added nodes.

16,384 scrambled-Sobol parameter sets (>= the 10,000 required).  Readout is always
the shared target protein p1.

Behaviours scored per parameter set:
   bistability      two CERTIFIED steady states at one input level (low-IC vs high-IC),
                    each integrated until max|dZ/dt|/max(|Z|,1e-3) < 1e-5; sets that
                    are still moving after 500 tau are reported as unresolved, not
                    as bistable
   sustained osc.   unstable focus at the operating point + non-decaying oscillation
   damped osc.      complex dominant eigenvalue with negative real part (ringing)
   ultrasensitivity effective Hill coefficient of the steady-state dose-response > 2
   memory           OFF-step half-time > 5 tau
   noise rejection  filter index > 0.5 (brief input rejected relative to a persistent one)
   pulse            non-monotone ON response with >10% overshoot
   dynamic range    max/min of the steady-state dose-response

Outputs: dynamics_higher_order_persetset.csv(.gz), dynamics_higher_order_fractions_compI1.csv,
         dynamics_higher_order_gain_loss_compI1.csv
"""
import sys, os, time, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP
from multiprocessing import Pool

RES   = "/path/to/revision/results/v3"
N_SET = 16384
DT    = 0.01
SEED  = 20260908
S_OP  = 2.0                      # operating input level
S_BIS = np.array([0.05, 0.3, 1.0, 3.0, 10.0])
S_DOSE = np.logspace(np.log10(0.02), np.log10(20.0), 11)
T_SETTLE = 50.0
T_DOSE = 20.0
T_BIS = 30.0
T_OSC = 150.0


# ------------------------------------------------------------------ helpers --
def jac_eigs(rhs, Z, S, active, h=1e-6):
    """Numerical Jacobian eigenvalues on the active states, vectorised over N."""
    k = len(active); N = Z.shape[1]
    f0 = rhs(0.0, Z, S)[active]                       # (k,N)
    J = np.empty((N, k, k))
    for j, idx in enumerate(active):
        Zp = Z.copy(); dh = h * np.maximum(np.abs(Z[idx]), 1.0)
        Zp[idx] = Zp[idx] + dh
        fp = rhs(0.0, Zp, S)[active]
        J[:, :, j] = ((fp - f0) / dh[None, :]).T
    ev = np.linalg.eigvals(J)                         # (N,k)
    order = np.argsort(-ev.real, axis=1)
    lead = np.take_along_axis(ev, order[:, :1], axis=1)[:, 0]
    return lead.real, np.abs(lead.imag), ev.real.max(axis=1)


def half_time(ts, Y, y_pre, y_ss):
    lvl = y_pre + 0.5 * (y_ss - y_pre)
    rising = y_ss >= y_pre
    n = Y.shape[1]
    out = np.full(n, np.nan)
    hit = np.where(rising[None, :], Y >= lvl[None, :], Y <= lvl[None, :])
    any_hit = hit.any(axis=0)
    idx = np.argmax(hit, axis=0)
    cols = np.where(any_hit)[0]
    out[cols] = ts[idx[cols]]
    return out


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
        seg = r[a:b+1]; sl = logS[a:b+1]
        if seg[-1] < seg[0]:
            seg = seg[::-1]; sl = sl[::-1]
        t10 = seg[0] + 0.1 * (seg[-1] - seg[0]); t90 = seg[0] + 0.9 * (seg[-1] - seg[0])
        s10 = np.interp(t10, seg, sl); s90 = np.interp(t90, seg, sl)
        d = abs(s90 - s10)
        if d > 1e-9:
            out[j] = np.log(81.0) / d
    return out


BIS_TOL   = 1e-5      # relative residual max_i |dZ_i/dt| / max(|Z_i|, 1e-3)
BIS_CHUNK = 25.0      # tau per settling chunk
BIS_TMAX  = 500.0     # tau; sets still moving after this are called UNRESOLVED
DR_FLOOR  = 1e-3      # detection floor for the dynamic range (fraction of max production)


def settle_adaptive(topo, P, Z0, S, act, tol=BIS_TOL, chunk=BIS_CHUNK, tmax=BIS_TMAX):
    """Integrate to a certified steady state, dropping columns as they converge.

    Returns (Z_final, relative_residual, time_used).  A column whose residual is
    still above `tol` at `tmax` has NOT reached a steady state; it is reported as
    unresolved rather than silently treated as a fixed point.
    """
    N = Z0.shape[1]
    Zf = Z0.copy()
    res = np.full(N, np.inf)
    tu = np.zeros(N)
    idx = np.arange(N)
    Z = Z0.copy()
    t = 0.0
    while len(idx) and t < tmax:
        Ps = {k: (v[idx] if np.ndim(v) else v) for k, v in P.items()}
        rhs = D.make_rhs(topo, Ps)
        Z = D.integrate(rhs, Z, lambda tt: S, chunk, DT, active=act)
        t += chunk
        d = rhs(0.0, Z, S)[act]
        sc = np.maximum(np.abs(Z[act]), 1e-3)
        r = (np.abs(d) / sc).max(axis=0)
        Zf[:, idx] = Z
        res[idx] = r
        tu[idx] = t
        keep = r >= tol
        idx = idx[keep]
        Z = Z[:, keep]
    return Zf, res, tu


# ------------------------------------------------------------------- worker --
def run_module(args):
    family, size, topo = args
    P = SMP.sample(N_SET, seed=SEED)
    P = {k: np.asarray(v, float) for k, v in P.items()}
    N = N_SET
    rhs = D.make_rhs(topo, P); act = topo.active_states()
    t0 = time.time()

    # --- settle low / high --------------------------------------------------
    Zlo = D.integrate(rhs, np.zeros((D.NSTATE, N)), lambda t: 0.10, T_SETTLE, DT, active=act)
    Zhi = D.integrate(rhs, np.zeros((D.NSTATE, N)), lambda t: S_OP, T_SETTLE, DT, active=act)
    p_lo = Zlo[D.IP1].copy(); p_hi = Zhi[D.IP1].copy()

    # --- ON step -------------------------------------------------------------
    Zon, ts, tr = D.integrate(rhs, Zlo, lambda t: S_OP, 50.0, DT, record_every=10,
                              active=act, record_states=(D.IP1,))
    p_on = tr[:, 0, :]; p_ss_on = Zon[D.IP1].copy()
    T_ON = half_time(ts, p_on, p_lo, p_ss_on)
    dyn = p_on - p_lo[None, :]; fin = p_ss_on - p_lo
    sg = np.sign(dyn[np.abs(dyn).argmax(axis=0), np.arange(N)]); sg[sg == 0] = 1
    ext = (dyn * sg[None, :]).max(axis=0) * sg
    with np.errstate(divide="ignore", invalid="ignore"):
        overshoot = np.where(np.abs(fin) > 1e-9, ext / fin, np.nan)
    is_pulse = ((np.abs(ext) > 1e-6) & (np.abs(ext) > 1.10 * np.abs(fin))).astype(int)
    d1 = np.diff(p_on, axis=0)
    n_turn = (np.diff(np.sign(d1), axis=0) != 0).sum(axis=0)

    # --- OFF step ------------------------------------------------------------
    Zoff, ts2, tr2 = D.integrate(rhs, Zhi, lambda t: 0.10, 50.0, DT, record_every=10,
                                 active=act, record_states=(D.IP1,))
    T_OFF = half_time(ts2, tr2[:, 0, :], p_hi, Zoff[D.IP1])

    # --- brief spurious input vs persistent ----------------------------------
    Sf = lambda t: (S_OP if t < 0.5 else 0.10)
    _, ts3, tr3 = D.integrate(rhs, Zlo, Sf, 50.0, DT, record_every=10,
                              active=act, record_states=(D.IP1,))
    amp_short = np.abs(tr3[:, 0, :] - p_lo[None, :]).max(axis=0)
    amp_long = np.abs(p_on - p_lo[None, :]).max(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        noise_tx = np.where(amp_long > 1e-9, amp_short / amp_long, np.nan)
    filt = 1.0 - noise_tx

    # --- steady-state dose-response (continuation) ---------------------------
    Z = np.zeros((D.NSTATE, N)); Ps = []
    for S in S_DOSE:
        Z = D.integrate(rhs, Z, lambda t: S, T_DOSE, DT, active=act)
        Ps.append(Z[D.IP1].copy())
    Ps = np.array(Ps)
    neff = eff_hill(S_DOSE, Ps)
    hiP = Ps.max(axis=0); loP = Ps.min(axis=0)
    dr = (hiP + 1e-9) / (loP + 1e-9)                      # raw (unbounded as loP -> 0)
    dr_floored = (hiP + DR_FLOOR) / (loP + DR_FLOOR)      # with a detection floor
    rngP = np.maximum(hiP - loP, 1e-12)
    dP = np.diff(Ps, axis=0)
    thr = 0.02 * rngP                                     # ignore steps < 2% of the range
    nonmono = ((dP > thr[None, :]).any(axis=0) &
               (dP < -thr[None, :]).any(axis=0)).astype(int)

    # --- bistability: two initial conditions per input level, each integrated to a
    #     CERTIFIED steady state.  The earlier continuation-based test could not tell
    #     a second attractor from a branch that had simply not finished relaxing, so
    #     every call is now conditioned on both branches converging.
    Zhi_ic = np.zeros((D.NSTATE, N)); Zhi_ic[act] = 3.0
    hyst = np.zeros(N); absgap = np.zeros(N)
    conv_all = np.ones(N, bool)
    worst_res = np.zeros(N)
    max_settle_t = np.zeros(N)
    for S in S_BIS:
        Zl, rl, tl = settle_adaptive(topo, P, np.zeros((D.NSTATE, N)), S, act)
        Zh, rh, th = settle_adaptive(topo, P, Zhi_ic.copy(), S, act)
        pl = Zl[D.IP1]; ph = Zh[D.IP1]
        den = np.maximum(np.maximum(np.abs(pl), np.abs(ph)), 1e-6)
        hyst = np.maximum(hyst, np.abs(pl - ph) / den)
        absgap = np.maximum(absgap, np.abs(pl - ph))
        conv_all &= (rl < BIS_TOL) & (rh < BIS_TOL)
        worst_res = np.maximum(worst_res, np.maximum(rl, rh))
        max_settle_t = np.maximum(max_settle_t, np.maximum(tl, th))
    is_bistable = ((hyst > 0.05) & (absgap > 1e-3) & conv_all).astype(int)
    is_unresolved = (~conv_all).astype(int)
    resid = worst_res
    # --- oscillation: eigenvalues at the operating point + long run ----------
    lead_re, lead_im, max_re = jac_eigs(rhs, Zhi, S_OP, act)
    Zl, tsl, trl = D.integrate(rhs, Zlo, lambda t: S_OP, T_OSC, DT, record_every=20,
                               active=act, record_states=(D.IP1,))
    tail = trl[int(0.75 * trl.shape[0]):, 0, :]
    pp = tail.max(axis=0) - tail.min(axis=0)
    rel_pp = pp / np.maximum(np.abs(tail.mean(axis=0)), 1e-6)
    is_sustained_osc = ((rel_pp > 0.02) & (max_re > 1e-4) & (lead_im > 1e-3)).astype(int)
    is_damped_osc = ((lead_im > 1e-3) & (lead_re < 0) & (n_turn >= 2)).astype(int)

    out = pd.DataFrame(dict(
        family=family, module=size, param_set=np.arange(N),
        p_lo=p_lo, p_hi=p_hi, T_half_ON=T_ON, T_half_OFF=T_OFF,
        overshoot_ratio=overshoot, is_pulse=is_pulse, n_turning_points=n_turn,
        noise_transmission=noise_tx, filter_index=filt,
        n_eff=neff, dyn_range=dr, dyn_range_floored=dr_floored,
        ss_min_protein=loP, ss_max_protein=hiP, nonmonotone_dose=nonmono,
        hysteresis=hyst, is_bistable=is_bistable,
        bistability_settle_residual=resid,
        bistability_unresolved=is_unresolved,
        bistability_max_settle_time_tau=max_settle_t,
        lead_eig_re=lead_re, lead_eig_im=lead_im, max_eig_re=max_re,
        rel_peak_to_peak=rel_pp,
        is_sustained_osc=is_sustained_osc, is_damped_osc=is_damped_osc,
        is_ultrasensitive=(neff > 2).astype(float),
        is_very_ultrasensitive=(neff > 4).astype(float),
        # a set whose OFF half-time exceeds the 50-tau window returns NaN; that is
        # the SLOWEST behaviour, so it counts as memory rather than as a failure
        is_memory=((T_OFF > 5.0) | (~np.isfinite(T_OFF))).astype(float),
        is_noise_rejecting=(filt > 0.5).astype(float),
    ))
    print(f"    {family}/{size} finished in {time.time()-t0:.0f}s", flush=True)
    return out


# --------------------------------------------------------------------- main --
if __name__ == "__main__":
    t0 = time.time()
    # Idempotence guard: this is a ~2 h sweep and it is launched from two places
    # (directly, and from 99_dynamics_driver.sh).  If the finished output is already
    # on disk, do nothing rather than recompute it.
    _done = f"{RES}/dynamics_higher_order_gain_loss_compI1.csv"
    if os.path.exists(_done):
        print(f"already complete: {_done} exists; nothing to do")
        sys.exit(0)
    tasks = []
    fams = {
        "COMP_I1_negfeedback": D.higher_order_registry("AND", "I1_MIR", composite=True),
        
    }
    for fam, reg in fams.items():
        for size in ("n3", "n4", "n5", "n6", "n4d"):
            tasks.append((fam, size, reg[size]))
    meta = [dict(family=f, module=s, topology=t.name, label=t.label, n_nodes=t.n_nodes,
                 gate=t.gate, gene_gene=t.gene_gene, mir_mir=t.mir_mir, tf_tf=t.tf_tf,
                 mir_to_TF=t.mir_to_TF, s_TM=t.s_TM, s_MY=t.s_MY, s_TY=t.s_TY)
            for f, s, t in tasks]
    pd.DataFrame(meta).to_csv(f"{RES}/dynamics_higher_order_modules_compI1.csv", index=False)

    with Pool(processes=min(4, len(tasks))) as pool:
        parts = pool.map(run_module, tasks)
    A = pd.concat(parts, ignore_index=True)
    A.to_csv(f"{RES}/dynamics_higher_order_persetset_compI1.csv.gz", index=False,
             compression="gzip")

    BIN = ["is_bistable", "bistability_unresolved", "is_sustained_osc", "is_damped_osc", "is_ultrasensitive",
           "is_very_ultrasensitive", "is_memory", "is_noise_rejecting", "is_pulse",
           "nonmonotone_dose"]
    CON = ["T_half_ON", "T_half_OFF", "n_eff", "dyn_range", "dyn_range_floored",
           "filter_index", "hysteresis", "overshoot_ratio", "rel_peak_to_peak"]

    frac = []
    for (fam, mod), g in A.groupby(["family", "module"]):
        r = dict(family=fam, module=mod, n_param_sets=len(g))
        for b in BIN:
            r[f"pct_{b}"] = 100 * np.nanmean(g[b].astype(float))
        for c in CON:
            r[f"{c}_median"] = np.nanmedian(g[c].astype(float))
            r[f"{c}_q90"] = np.nanpercentile(g[c].astype(float), 90)
        frac.append(r)
    pd.DataFrame(frac).to_csv(f"{RES}/dynamics_higher_order_fractions_compI1.csv", index=False)

    # ---- matched gain / loss versus the embedded n=3 core --------------------
    from scipy.stats import wilcoxon, binomtest
    gl = []
    for fam in A.family.unique():
        base = A[(A.family == fam) & (A.module == "n3")].set_index("param_set")
        for mod in ("n4", "n5", "n6", "n4d"):
            sub = A[(A.family == fam) & (A.module == mod)].set_index("param_set")
            for b in BIN:
                a = base[b].astype(float).values; c = sub[b].astype(float).values
                gain = int(((c > 0.5) & (a < 0.5)).sum())
                loss = int(((c < 0.5) & (a > 0.5)).sum())
                both = int(((c > 0.5) & (a > 0.5)).sum())
                nn = len(a)
                # exact McNemar
                p = binomtest(gain, gain + loss, 0.5).pvalue if (gain + loss) > 0 else 1.0
                gl.append(dict(family=fam, module=mod, behaviour=b, n=nn,
                               pct_core=100 * np.nanmean(a), pct_module=100 * np.nanmean(c),
                               n_gain=gain, n_loss=loss, n_both=both,
                               pct_gain=100 * gain / nn, pct_loss=100 * loss / nn,
                               mcnemar_p=p, metric_type="binary"))
            for c_ in CON:
                a = base[c_].astype(float).values; c = sub[c_].astype(float).values
                ok = np.isfinite(a) & np.isfinite(c)
                if ok.sum() > 20 and np.any(a[ok] != c[ok]):
                    st = wilcoxon(c[ok], a[ok], zero_method="zsplit")
                    pv = st.pvalue
                else:
                    pv = np.nan
                gl.append(dict(family=fam, module=mod, behaviour=c_, n=int(ok.sum()),
                               pct_core=np.nanmedian(a), pct_module=np.nanmedian(c),
                               n_gain=int((c[ok] > a[ok]).sum()), n_loss=int((c[ok] < a[ok]).sum()),
                               n_both=np.nan,
                               pct_gain=100 * np.nanmean(c[ok] > a[ok]),
                               pct_loss=100 * np.nanmean(c[ok] < a[ok]),
                               mcnemar_p=pv, metric_type="continuous_wilcoxon"))
    pd.DataFrame(gl).to_csv(f"{RES}/dynamics_higher_order_gain_loss_compI1.csv", index=False)
    print(f"\nDONE in {(time.time()-t0)/60:.1f} min")
    print(pd.DataFrame(frac)[["family", "module", "pct_is_bistable", "pct_bistability_unresolved",
                              "pct_is_sustained_osc",
                              "pct_is_ultrasensitive", "pct_is_memory",
                              "pct_is_noise_rejecting"]].to_string(index=False))
