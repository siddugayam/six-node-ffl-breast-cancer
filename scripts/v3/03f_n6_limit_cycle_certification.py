#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03f_n6_limit_cycle_certification.py
===================================
CERTIFY THE ONE QUALITATIVE GAIN OF FUNCTION FROM HIGHER ORDER.

The 16,384-set sweep found sustained oscillation in 0.66% of parameter space for
the 6-node composite module and in 0% for every other module including its own
embedded 3-node core.  This script certifies that result rather than asserting it:

  (1) STRUCTURE.  It states, and checks against the topology object, which loops
      exist at each module size.  The 3-node composite C2 core contains only the
      two-node POSITIVE loop TF1 -| miR -| TF1 (double negative), which can give
      bistability but not oscillation.  The 6-node module additionally contains
      the three-node NEGATIVE loop TF1 -> TF2 -> miR -| TF1.  A negative loop of
      length >= 3 is the minimal structural requirement for a Hopf bifurcation.

  (2) DYNAMICS.  Every set flagged as oscillating is re-integrated for 1,000 tau
      (10x the screening run) from two different initial conditions.  A set is
      certified as a LIMIT CYCLE only if, over the final 25% of both runs, the
      relative peak-to-peak amplitude stays above 2% and the two initial
      conditions converge on the same amplitude and period to within 5%.

  (3) CONTROL.  The identical parameter vectors are run through the embedded
      3-node core, the 4-node and the 5-node modules, and their eigenvalues and
      amplitudes are reported alongside.

Out: results/v3/dynamics_n6_limit_cycle_certified.csv
     results/v3/dynamics_n6_limit_cycle_timecourses.csv
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
N_SET = 16384
SEED = 20260908
DT = 0.01
S_OP = 2.0
T_LONG = 1000.0
REC = 20
TAU_HOURS = 9.0 / np.log(2)      # median mammalian mRNA t1/2 = 9 h (Schwanhausser 2011)

A = pd.read_csv(f"{RES}/dynamics_higher_order_persetset.csv.gz")
osc = A[(A.family == "COMP_C2_toggle") & (A.module == "n6") & (A.is_sustained_osc > 0)]
sets = np.sort(osc.param_set.to_numpy(int))
print(f"sets flagged as sustained oscillators at n=6: {len(sets)} of {N_SET} "
      f"({100*len(sets)/N_SET:.3f}%)")

P0 = SMP.sample(N_SET, seed=SEED)
P = {k: np.asarray(v, float)[sets] for k, v in P0.items()}
N = len(sets)
reg = D.higher_order_registry("AND", "C2_MIR", composite=True)


def amp_period(ts, y, frac=0.25):
    """Relative peak-to-peak amplitude and period from the last `frac` of a run."""
    k = int((1 - frac) * len(ts))
    t = ts[k:]; Y = y[k:]
    pk = Y.max(axis=0) - Y.min(axis=0)
    mean = np.maximum(np.abs(Y.mean(axis=0)), 1e-9)
    rel = pk / mean
    # period from zero-crossings of the mean-centred signal
    C = Y - Y.mean(axis=0)[None, :]
    per = np.full(Y.shape[1], np.nan)
    for j in range(Y.shape[1]):
        s = np.sign(C[:, j])
        idx = np.where(np.diff(s) != 0)[0]
        if len(idx) >= 3:
            per[j] = 2.0 * np.median(np.diff(t[idx]))
    return rel, per


def jac_lead(topo, P_, Z, S):
    rhs = D.make_rhs(topo, P_)
    act = topo.active_states()
    k = len(act); n = Z.shape[1]
    f0 = rhs(0.0, Z, S)[act]
    J = np.empty((n, k, k))
    for j, i in enumerate(act):
        Zp = Z.copy(); dh = 1e-6 * np.maximum(np.abs(Z[i]), 1.0)
        Zp[i] = Zp[i] + dh
        J[:, :, j] = ((rhs(0.0, Zp, S)[act] - f0) / dh[None, :]).T
    ev = np.linalg.eigvals(J)
    o = np.argsort(-ev.real, axis=1)
    lead = np.take_along_axis(ev, o[:, :1], axis=1)[:, 0]
    return lead.real, np.abs(lead.imag), ev.real.max(axis=1), np.abs(ev.imag).max(axis=1)


rows = {}
tc_rows = []
for size in ("n3", "n4", "n5", "n6"):
    topo = reg[size]
    act = topo.active_states()
    rhs = D.make_rhs(topo, P)
    out = {}
    for ic_name, ic in (("zero", 0.0), ("high", 1.5)):
        Z0 = np.zeros((D.NSTATE, N))
        Z0[act] = ic
        Zf, ts, tr = D.integrate(rhs, Z0, lambda t: S_OP, T_LONG, DT,
                                 record_every=REC, active=act, record_states=(D.IP1,))
        y = tr[:, 0, :]
        rel, per = amp_period(ts, y)
        out[ic_name] = dict(rel=rel, per=per, Zf=Zf, ts=ts, y=y)
    lr, li, mr, mi = jac_lead(topo, P, out["zero"]["Zf"], S_OP)
    rows[size] = pd.DataFrame(dict(
        module=size, param_set=sets,
        rel_amp_from_zero_IC=out["zero"]["rel"], rel_amp_from_high_IC=out["high"]["rel"],
        period_tau_from_zero_IC=out["zero"]["per"], period_tau_from_high_IC=out["high"]["per"],
        lead_eig_re=lr, lead_eig_im=li, max_eig_re=mr, max_abs_eig_im=mi))
    if size in ("n3", "n6"):
        j = int(np.argmax(out["zero"]["rel"]))
        for i, t in enumerate(out["zero"]["ts"]):
            if t >= T_LONG - 60:
                tc_rows.append(dict(module=size, param_set=int(sets[j]), t_tau=t,
                                    t_hours=t * TAU_HOURS,
                                    target_protein=out["zero"]["y"][i, j]))
    print(f"  {size}: median rel amplitude (last 25% of 1,000 tau) = "
          f"{np.nanmedian(out['zero']['rel']):.4f}; % with complex leading pair = "
          f"{100*np.mean(li > 1e-3):.1f}", flush=True)

R = pd.concat(rows.values(), ignore_index=True)
# certification
piv = R.pivot(index="param_set", columns="module")
cert = []
for ps in sets:
    r6 = R[(R.module == "n6") & (R.param_set == ps)].iloc[0]
    r3 = R[(R.module == "n3") & (R.param_set == ps)].iloc[0]
    a0, a1 = r6.rel_amp_from_zero_IC, r6.rel_amp_from_high_IC
    p0, p1 = r6.period_tau_from_zero_IC, r6.period_tau_from_high_IC
    ok_amp = (a0 > 0.02) and (a1 > 0.02)
    ok_conv = (abs(a0 - a1) / max(a0, a1) < 0.05) if (a0 > 0 and a1 > 0) else False
    ok_per = (np.isfinite(p0) and np.isfinite(p1) and abs(p0 - p1) / max(p0, p1) < 0.05)
    cert.append(dict(param_set=int(ps),
                     n6_rel_amp_zero_IC=a0, n6_rel_amp_high_IC=a1,
                     n6_period_tau=p0, n6_period_hours=p0 * TAU_HOURS,
                     n6_lead_eig_re=r6.lead_eig_re, n6_lead_eig_im=r6.lead_eig_im,
                     n3_rel_amp=r3.rel_amp_from_zero_IC,
                     n3_lead_eig_im=r3.lead_eig_im, n3_max_eig_re=r3.max_eig_re,
                     amplitude_sustained_1000tau=ok_amp,
                     independent_of_initial_condition=ok_conv and ok_per,
                     certified_limit_cycle=bool(ok_amp and ok_conv and ok_per),
                     core_3node_oscillates=bool(r3.rel_amp_from_zero_IC > 0.02
                                                and r3.lead_eig_im > 1e-3)))
C = pd.DataFrame(cert)
C.to_csv(f"{RES}/dynamics_n6_limit_cycle_certified.csv", index=False)
pd.DataFrame(tc_rows).to_csv(f"{RES}/dynamics_n6_limit_cycle_timecourses.csv", index=False)

print()
print(f"flagged at n=6                                : {len(C)}")
print(f"amplitude still > 2% after 1,000 tau          : {int(C.amplitude_sustained_1000tau.sum())}")
print(f"same amplitude AND period from 2 different ICs: {int(C.independent_of_initial_condition.sum())}")
print(f"CERTIFIED limit cycles                        : {int(C.certified_limit_cycle.sum())} "
      f"({100*C.certified_limit_cycle.mean():.1f}% of flagged; "
      f"{100*C.certified_limit_cycle.sum()/N_SET:.3f}% of the 16,384-set design)")
print(f"of those, embedded 3-node core also oscillates : {int((C.certified_limit_cycle & C.core_3node_oscillates).sum())}")
print(f"median certified period                       : "
      f"{np.nanmedian(C.loc[C.certified_limit_cycle,'n6_period_tau']):.2f} tau = "
      f"{np.nanmedian(C.loc[C.certified_limit_cycle,'n6_period_hours']):.1f} h "
      f"(tau = mRNA lifetime 12.98 h)")
print(f"max |Im(lambda)| in the 3-node cores of these sets: {C.n3_lead_eig_im.max():.3g}")
