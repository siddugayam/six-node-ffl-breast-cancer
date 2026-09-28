#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03c_wellposedness_check.py
--------------------------
The 4-, 5- and 6-node modules add an undirected STRING gene-gene edge, which we read
generously as mutual protein stabilisation (complex formation protects both partners
from degradation).  Written in the obvious linear way,

        d p1/dt = gp * ( y1 - p1 /(1 + kappa*p2) )
        d p2/dt = gp * ( y2 - p2 /(1 + kappa*p1) )

the pair has NO positive steady state whenever  kappa^2 * y1 * y2 > 1  (solve the two
steady-state equations: p1*(1 - kappa^2*y1*y2) = y1 + kappa*y1*y2).  In that regime the
module runs away to infinity rather than settling, and any two-initial-condition
"bistability" test will report a large gap between branches purely because neither
branch has converged.

This script quantifies that failure over the sampled parameter space and shows that the
saturating form actually used in dyn_models.py,

        d p1/dt = gp * ( y1 - p1 /(1 + kappa*f+(p2; Kg12, ng12)) )

is bounded (protection at most 1+kappa) everywhere.

Out: results/v3/dynamics_model_wellposedness_check.csv
"""
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
N = 4096
DT = 0.01
S_OP = 2.0

P = SMP.sample(16384, seed=20260908)
P = {k: np.asarray(v, float)[:N] for k, v in P.items()}

rows = []
for variant in ("linear_legacy", "saturating"):
    D.MUTUAL_STABILISATION = variant
    for fam, base, comp in (("I1_miRNA_FFL", "I1_MIR", False),
                            ("COMP_C2_toggle", "C2_MIR", True)):
        reg = D.higher_order_registry("AND", base, composite=comp)
        for size in ("n3", "n4", "n4d", "n5", "n6"):
            topo = reg[size]
            act = topo.active_states()
            rhs = D.make_rhs(topo, P)
            t0 = time.time()
            Z1 = D.integrate(rhs, np.zeros((D.NSTATE, N)), lambda t: S_OP, 50.0, DT, active=act)
            Z2 = D.integrate(rhs, Z1, lambda t: S_OP, 200.0, DT, active=act)
            p1a, p1b = Z1[D.IP1].copy(), Z2[D.IP1].copy()
            resid = np.abs(rhs(0.0, Z2, S_OP)[act]).max(axis=0)
            scale = np.maximum(np.abs(Z2[act]), 1e-3)
            rel_resid = (np.abs(rhs(0.0, Z2, S_OP)[act]) / scale).max(axis=0)
            grow = p1b / np.maximum(p1a, 1e-12)
            rows.append(dict(
                mutual_stabilisation=variant, family=fam, module=size,
                n_param_sets=N,
                p1_median=float(np.median(p1b)), p1_max=float(np.nanmax(p1b)),
                pct_p1_gt_10=100 * float(np.mean(p1b > 10)),
                pct_p1_gt_100=100 * float(np.mean(p1b > 100)),
                pct_still_growing_1p5x_from_50_to_250tau=100 * float(np.mean(grow > 1.5)),
                median_abs_residual_at_250tau=float(np.median(resid)),
                pct_rel_residual_gt_1e4=100 * float(np.mean(rel_resid > 1e-4)),
                analytic_runaway_condition_kappa2_y1_y2_gt_1=(
                    100 * float(np.mean(P["kappa"] ** 2 * Z2[D.IY1] * Z2[D.IY2] > 1))
                    if topo.gene_gene == "mutual" else np.nan),
                seconds=round(time.time() - t0, 1)))
            print(f"  {variant:14s} {fam:15s} {size:4s} "
                  f"p1max={rows[-1]['p1_max']:.3g} "
                  f"unconverged={rows[-1]['pct_rel_residual_gt_1e4']:.1f}%", flush=True)

D.MUTUAL_STABILISATION = "saturating"
W = pd.DataFrame(rows)
W.to_csv(f"{RES}/dynamics_model_wellposedness_check.csv", index=False)
print()
print(W[["mutual_stabilisation", "family", "module", "p1_max",
         "pct_still_growing_1p5x_from_50_to_250tau", "pct_rel_residual_gt_1e4",
         "analytic_runaway_condition_kappa2_y1_y2_gt_1"]].to_string(index=False))
print("\nwritten:", f"{RES}/dynamics_model_wellposedness_check.csv")
