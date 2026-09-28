#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09b_composite_oscillation.py
----------------------------
WHY A COMPOSITE FFL CAN DO SOMETHING NO TRANSCRIPTIONAL FFL CAN.

Alon's feed-forward loop is a directed ACYCLIC graph: X -> Y, X -> Z, Y -> Z.  Order
the states (X, Y, Z); the Jacobian is lower triangular, so its eigenvalues are exactly
the diagonal entries -- the (negative, real) decay rates.  A transcriptional FFL
therefore CANNOT ring, cannot oscillate and cannot be multistable, whatever its
parameters and whatever the shape of its input functions.  This is not a numerical
result, it is a structural one.

Two features of the miRNA circuits break the triangular structure:

  (1) the reciprocal miRNA -| TF arm of the COMPOSITE FFL (1,434 of our 1,649 cores)
      closes a two-node loop: TF -> miRNA -| TF (delayed negative feedback, 791 cores)
      or TF -| miRNA -| TF (double-negative toggle, 451 cores);

  (2) stoichiometric TITRATION.  When miRNA and target are removed in a 1:1 complex
      (-theta*m*y in BOTH equations, Levine et al. 2007 doi:10.1371/journal.pbio.0050229;
      Mukherji et al. 2011 doi:10.1038/ng.905) the target feeds back on the miRNA even
      when the miRNA arm is the only connection.  A catalytic miRNA (theta = 0) leaves
      the loop acyclic; a stoichiometric one does not.

This script measures, per topology and per parameter set:
  * whether the leading eigenvalue of the Jacobian at the operating point is complex
    (capacity for ringing / damped oscillation),
  * whether it is complex with positive real part (Hopf / sustained oscillation),
  * the oscillation period implied by the imaginary part,
  * the relative peak-to-peak amplitude of the target protein over the last quarter of
    a 300-tau run (numerical confirmation),
and splits the miRNA circuits by whether titration is switched on (theta > 0).

Out: results/v3/dynamics_3node_oscillation.csv
     results/v3/dynamics_3node_oscillation_by_titration.csv
"""
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D
from multiprocessing import Pool

RES = "/path/to/revision/results/v3"
DT = 0.01
S_OP = 2.0
T_RUN = 300.0
TAU_HOURS = 9.0 / np.log(2)      # median mammalian mRNA t1/2 9 h (Schwanhausser 2011)

Pall = dict(np.load(f"{RES}/dynamics_3node_sobol_params.npz"))
NT = len(Pall["gx"])
reg = D.build_registry()


def settle_adaptive(topo, P, Z0, S, act, tol=1e-7, chunk=25.0, tmax=400.0):
    N = Z0.shape[1]
    Zf = Z0.copy(); res = np.full(N, np.inf)
    idx = np.arange(N); Z = Z0.copy(); t = 0.0
    while len(idx) and t < tmax:
        Ps = {k: (v[idx] if np.ndim(v) else v) for k, v in P.items()}
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


def jac_eigs(rhs, Z, S, act, h=1e-6):
    k = len(act); N = Z.shape[1]
    f0 = rhs(0.0, Z, S)[act]
    J = np.empty((N, k, k))
    for j, idx in enumerate(act):
        Zp = Z.copy(); dh = h * np.maximum(np.abs(Z[idx]), 1.0)
        Zp[idx] = Zp[idx] + dh
        fp = rhs(0.0, Zp, S)[act]
        J[:, :, j] = ((fp - f0) / dh[None, :]).T
    ev = np.linalg.eigvals(J)
    order = np.argsort(-ev.real, axis=1)
    lead = np.take_along_axis(ev, order[:, :1], axis=1)[:, 0]
    return lead.real, np.abs(lead.imag), ev.real.max(axis=1), np.abs(ev.imag).max(axis=1)


def run_one(args):
    """mode = 'sampled'  theta drawn from the Sobol design (stoichiometric titration on)
       mode = 'zero'     theta forced to 0 (purely catalytic RISC: the ONLY difference)"""
    nm, mode = args
    t0 = time.time()
    topo = reg[nm]
    act = topo.active_states()
    P = dict(Pall)
    if mode == "zero":
        P["theta"] = np.zeros(NT)
    Z, res = settle_adaptive(topo, P, np.zeros((D.NSTATE, NT)), S_OP, act)
    rhs = D.make_rhs(topo, P)
    lead_re, lead_im, max_re, max_im = jac_eigs(rhs, Z, S_OP, act)
    _, ts, tr = D.integrate(rhs, np.zeros((D.NSTATE, NT)), lambda t: S_OP, T_RUN, DT,
                            record_every=20, active=act, record_states=(D.IP1,))
    tail = tr[int(0.75 * tr.shape[0]):, 0, :]
    pp = tail.max(axis=0) - tail.min(axis=0)
    rel_pp = pp / np.maximum(np.abs(tail.mean(axis=0)), 1e-6)
    with np.errstate(divide="ignore", invalid="ignore"):
        period = np.where(lead_im > 1e-6, 2 * np.pi / np.maximum(lead_im, 1e-12), np.nan)
    df = pd.DataFrame(dict(
        topology=nm, titration_mode=mode, census_class=topo.core_class,
        n_sign_resolved_cores_in_BRCA_network=topo.n_cores_in_network,
        acyclic_without_titration=not topo.mir_to_TF,
        param_set=np.arange(NT) - 1,
        theta=P["theta"], settle_residual=res,
        lead_eig_re=lead_re, lead_eig_im=lead_im,
        max_eig_re=max_re, max_abs_eig_im=max_im,
        complex_leading_pair=(lead_im > 1e-6).astype(int),
        any_complex_eigenvalue=(max_im > 1e-6).astype(int),
        hopf_unstable_focus=((max_re > 1e-6) & (max_im > 1e-6)).astype(int),
        oscillation_period_tau=period,
        oscillation_period_hours=period * TAU_HOURS,
        rel_peak_to_peak_last_quarter=rel_pp,
        sustained_oscillation=((rel_pp > 0.02) & (max_re > 1e-6) & (max_im > 1e-6)).astype(int),
    ))
    print(f"  {nm:16s} theta={mode:7s} done ({time.time()-t0:5.1f}s)  complex={100*df.any_complex_eigenvalue.mean():5.1f}%"
          f"  hopf={100*df.hopf_unstable_focus.mean():4.1f}%", flush=True)
    return df


if __name__ == "__main__":
    t0 = time.time()
    tasks = [(nm, mode) for mode in ("sampled", "zero") for nm in reg]
    with Pool(processes=8) as pool:
        parts = pool.map(run_one, tasks)
    A = pd.concat(parts, ignore_index=True)
    A.to_csv(f"{RES}/dynamics_3node_oscillation.csv", index=False)

    E = A[A.param_set >= 0].copy()
    E["titration"] = np.where(E.titration_mode == "zero",
                              "theta = 0 (catalytic RISC)",
                              "theta sampled 0-10 (stoichiometric titration)")
    rows = []
    for (nm, tit), g in E.groupby(["topology", "titration"]):
        t = reg[nm]
        rows.append(dict(
            topology=nm, census_class=t.core_class,
            n_sign_resolved_cores_in_BRCA_network=t.n_cores_in_network,
            has_reciprocal_miRNA_to_TF_arm=t.mir_to_TF,
            miRNA_arm_post_transcriptional=t.mir_arm,
            titration=tit, n_param_sets=len(g),
            pct_any_complex_eigenvalue=100 * g.any_complex_eigenvalue.mean(),
            pct_complex_leading_pair=100 * g.complex_leading_pair.mean(),
            pct_hopf_unstable_focus=100 * g.hopf_unstable_focus.mean(),
            pct_sustained_oscillation=100 * g.sustained_oscillation.mean(),
            median_period_tau=np.nanmedian(g.oscillation_period_tau),
            median_period_hours=np.nanmedian(g.oscillation_period_hours),
            median_rel_peak_to_peak=np.nanmedian(g.rel_peak_to_peak_last_quarter),
            pct_settle_unresolved=100 * np.mean(g.settle_residual > 1e-7)))
    B = pd.DataFrame(rows).sort_values(["topology", "titration"])
    B.to_csv(f"{RES}/dynamics_3node_oscillation_by_titration.csv", index=False)
    print(f"\nDONE in {time.time()-t0:.0f}s")
    print(B[["topology", "titration", "pct_any_complex_eigenvalue", "pct_hopf_unstable_focus",
             "pct_sustained_oscillation", "median_period_hours"]].to_string(index=False))
