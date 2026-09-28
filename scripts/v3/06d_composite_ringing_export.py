#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06d_composite_ringing_export.py
===============================
CORRECTION TO THE RINGING FIGURE.

06b_composite_showcase.py picked its ringing example by the largest
|Im| / |Re| among ALL eigenvalues.  A complex pair that is not the LEADING
eigenvalue decays before it can be seen, and the exported trajectory for the
chosen set turned out to be monotone (0 turning points, peak 2.26 falling
straight to the steady state) - an overshoot, not ringing.

Here the example is chosen on the correct criterion: the LEADING (slowest)
eigenvalue must be complex, and among those sets we take the one with the
largest |Im| / |Re|, i.e. the most under-damped.  Selection uses the
already-computed eigenvalues in dynamics_3node_oscillation.csv, restricted to
the catalytic-RISC case (theta = 0) so that the demonstration does not depend on
stoichiometric titration.  The number of turning points in the exported
trajectory is written out, so the figure caption can be checked against it.

Out: results/v3/dynamics_composite_ringing_timecourses.csv  (overwrites)
     results/v3/dynamics_composite_ringing_example.csv
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D

RES = "/path/to/revision/results/v3"
DT = 0.01
TAU_HOURS = 9.0 / np.log(2)

O = pd.read_csv(f"{RES}/dynamics_3node_oscillation.csv")
g = O[(O.topology == "COMP_I1_AND") & (O.titration_mode == "zero") &
      (O.param_set >= 0) & (O.lead_eig_im > 1e-3) & (O.settle_residual < 1e-7)].copy()
print(f"COMP_I1_AND sets with a COMPLEX LEADING eigenvalue at theta = 0: {len(g)} "
      f"of {int((O.topology=='COMP_I1_AND').sum()/2)}")
g["underdamping"] = g.lead_eig_im / np.abs(g.lead_eig_re)
g = g.sort_values("underdamping", ascending=False)
pick = int(g.iloc[0].param_set)
print(f"most under-damped set: {pick}  |Im|/|Re| = {float(g.iloc[0].underdamping):.2f}  "
      f"period = {2*np.pi/float(g.iloc[0].lead_eig_im):.2f} tau "
      f"= {2*np.pi/float(g.iloc[0].lead_eig_im)*TAU_HOURS:.1f} h")

Z = np.load(f"{RES}/dynamics_3node_sobol_params.npz")
P = {k: Z[k][[pick + 1]] for k in Z.files}     # column 0 is the reference set
P["theta"] = np.zeros(1)                        # catalytic RISC, as selected
reg = D.build_registry()

rows, summ = [], []
for nm in ("COMP_I1_AND", "I1_MIR_AND", "I1_TXN_AND", "CASCADE"):
    tp = reg[nm]; act = tp.active_states()
    rhs = D.make_rhs(tp, P)
    Zlo = D.integrate(rhs, np.zeros((D.NSTATE, 1)), lambda t: 0.10, 120.0, DT, active=act)
    Zf, ts, tr = D.integrate(rhs, Zlo, lambda t: 2.0, 60.0, DT, record_every=5,
                             active=act, record_states=(D.IP1,))
    y = tr[:, 0, 0]
    ss = float(Zf[D.IP1, 0])
    rel = y / max(ss, 1e-12)
    # Count turning points only where the step is larger than 1e-6 of the response
    # range; without that, a flat trajectory registers dozens of sign flips of
    # machine-epsilon differences.
    d = np.diff(y)
    thr = 1e-6 * max(y.max() - y.min(), 1e-12)
    dsig = np.where(np.abs(d) > thr, np.sign(d), 0.0)
    nz = dsig[dsig != 0]
    turns = int((np.diff(nz) != 0).sum()) if len(nz) > 1 else 0
    for t, v, r in zip(ts, y, rel):
        rows.append(dict(topology=nm, param_set=pick, t=t, t_hours=t * TAU_HOURS,
                         protein=float(v), protein_over_steady_state=float(r)))
    summ.append(dict(topology=nm, param_set=pick,
                     has_reciprocal_miRNA_to_TF_arm=tp.mir_to_TF,
                     n_turning_points_in_ON_response=turns,
                     peak_over_steady_state=float(rel.max()),
                     trough_after_peak_over_steady_state=float(
                         rel[int(np.argmax(rel)):].min()),
                     undershoots_below_steady_state=bool(rel[int(np.argmax(rel)):].min() < 0.999),
                     steady_state=ss))
pd.DataFrame(rows).to_csv(f"{RES}/dynamics_composite_ringing_timecourses.csv", index=False)
S = pd.DataFrame(summ)
S.to_csv(f"{RES}/dynamics_composite_ringing_example.csv", index=False)
pd.set_option("display.width", 250)
print()
print(S.round(4).to_string(index=False))
