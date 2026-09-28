#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_analytic_bounds.py
Three analytic consequences of putting the repressing arm POST-transcriptionally,
each derived in closed form and then checked numerically against the simulator.
Out: dynamics_analytic_bounds.csv
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
DT, N = 0.01, 512
P = {k: np.asarray(v, float) for k, v in SMP.sample(N, seed=99).items()}
P["theta"] = np.zeros(N)          # bounds below are derived for the catalytic regime
rows = []

# ---------------------------------------------------------------------------
# BOUND 1.  Maximum repression of the target by the miRNA arm.
#   OR (SUM) gate:   y_ss = [(u_TF + b)/(1+b)] / (1 + lam_y * q(m))
#   q in [0,1]  =>   y_ss(m=0) / y_ss(m->inf) = 1 + lam_y   EXACTLY.
#   A promoter-bound repressor with a Hill function drives production to 0, so its
#   repression is unbounded.  The miRNA arm is an ATTENUATOR, not a logic NOT.
# ---------------------------------------------------------------------------
topo = D.build_registry()["I1_MIR_OR"]
# For B1 and B2 the shared-RISC capacity is removed (Kc -> inf) so the bound is the
# clean 1 + lambda_y.  With a FINITE capacity Kc the effective repressor pool
# saturates at Kc and the achievable repression falls to 1 + lambda_y*q(Kc); that is
# reported separately below.
P1 = dict(P); P1["Kc"] = np.full(N, 1e9)
rhs = D.make_rhs(topo, P1); act = topo.active_states()
Z = np.zeros((D.NSTATE, N)); Z[D.IM1] = 0.0
# clamp the miRNA at 0 and at a saturating level by disabling its own dynamics
class Clamp:
    def __init__(self, base, mval): self.base, self.m = base, mval
    def __call__(self, t, Zz, S):
        Zz = Zz.copy(); Zz[D.IM1] = self.m
        dZ = self.base(t, Zz, S); dZ[D.IM1] = 0.0
        return dZ
lo = D.integrate(Clamp(rhs, np.zeros(N)), np.zeros((D.NSTATE, N)), lambda t: 2.0, 80, DT, active=act)
hi = D.integrate(Clamp(rhs, np.full(N, 1e4)), np.zeros((D.NSTATE, N)), lambda t: 2.0, 80, DT, active=act)
obs = lo[D.IY1] / np.maximum(hi[D.IY1], 1e-12)
pred = 1.0 + P1["lam_y"]
err = np.abs(obs - pred) / pred
rows.append(dict(bound="B1_max_repression_by_miRNA_arm",
                 statement="OR/SUM gate: y_ss(m=0)/y_ss(m->inf) = 1 + lambda_y exactly; the "
                           "miRNA arm attenuates by at most (1+lambda_y)-fold and can never "
                           "drive the target to zero, unlike a promoter-bound repressor",
                 n_checked=N, median_predicted=float(np.median(pred)),
                 median_observed=float(np.median(obs)),
                 max_relative_error=float(np.nanmax(err)),
                 median_relative_error=float(np.nanmedian(err))))

# ---------------------------------------------------------------------------
# BOUND 2.  Relaxation rate of the target transcript.
#   dy/dt = prod - (1 + lam_y*q(m)) y   =>  the target's relaxation RATE is
#   (1 + lam_y*q(m)) in units of its own basal decay rate, so its response
#   half-time is  ln2 / (1 + lam_y*q(m))  -- STRICTLY FASTER than the unregulated
#   transcript whenever the miRNA is present.  Post-transcriptional repression
#   accelerates the target REGARDLESS of the loop's coherence.  Alon's C1 delays.
# ---------------------------------------------------------------------------
mfix = np.full(N, 1.0)
class Clamp2:
    """Clamp both the miRNA and the TF so the target sees an INSTANTANEOUS production
    step; the measured half-time is then the target's own relaxation time alone."""
    def __init__(self, base, mval, xval): self.base, self.m, self.x = base, mval, xval
    def __call__(self, t, Zz, S):
        Zz = Zz.copy(); Zz[D.IM1] = self.m; Zz[D.IX1] = self.x
        dZ = self.base(t, Zz, S); dZ[D.IM1] = 0.0; dZ[D.IX1] = 0.0
        return dZ
xfix = np.full(N, 1.0)
rhs2 = Clamp2(rhs, mfix, xfix)
Z0 = np.zeros((D.NSTATE, N))
Zf, ts, tr = D.integrate(rhs2, Z0, lambda t: 2.0, 40, DT, record_every=2, active=act,
                         record_states=(D.IY1,))
Y = tr[:, 0, :]; y0 = np.zeros(N); yss = Zf[D.IY1]
lvl = y0 + 0.5 * (yss - y0)
hitmask = Y >= lvl[None, :]
idx = np.argmax(hitmask, axis=0)
t_half = np.where(hitmask.any(axis=0), ts[idx], np.nan)
q = 1.0 - (P1["Kmy"] ** P1["nmy"]) / (P1["Kmy"] ** P1["nmy"] + mfix ** P1["nmy"])
pred2 = np.log(2) / (1.0 + P1["lam_y"] * q)
err2 = np.abs(t_half - pred2) / pred2
rows.append(dict(bound="B2_target_relaxation_rate",
                 statement="target-mRNA response half-time = ln2 / (1 + lambda_y*q(m)); the "
                           "post-transcriptional arm accelerates the target for ANY loop sign, "
                           "whereas Alon's coherent C1 delays it",
                 n_checked=N, median_predicted=float(np.nanmedian(pred2)),
                 median_observed=float(np.nanmedian(t_half)),
                 max_relative_error=float(np.nanmax(err2)),
                 median_relative_error=float(np.nanmedian(err2))))

# ---------------------------------------------------------------------------
# BOUND 3.  Maximum overshoot of the miRNA-mediated I1 pulse.
#   Before the miRNA rises the target relaxes towards prod/1; afterwards towards
#   prod/(1+lam_y*q_inf).  The peak therefore cannot exceed (1+lam_y*q_inf) times
#   the final steady state, i.e. the pulse amplitude is CAPPED BY lambda_y.  A
#   transcriptional I1 has no such cap: its repressor can shut the promoter, so its
#   overshoot ratio diverges as the repression deepens.
# ---------------------------------------------------------------------------
t3 = D.build_registry()["I1_MIR_OR"]
r3 = D.make_rhs(t3, P); a3 = t3.active_states()
Zl = D.integrate(r3, np.zeros((D.NSTATE, N)), lambda t: 0.05, 80, DT, active=a3)
Zf3, ts3, tr3 = D.integrate(r3, Zl, lambda t: 3.0, 60, DT, record_every=2, active=a3,
                            record_states=(D.IY1,))
Yv = tr3[:, 0, :]
over = Yv.max(axis=0) / np.maximum(Zf3[D.IY1], 1e-12)
qinf = 1.0 - (P["Kmy"] ** P["nmy"]) / (P["Kmy"] ** P["nmy"] + np.maximum(Zf3[D.IM1], 1e-12) ** P["nmy"])
cap = 1.0 + P["lam_y"] * qinf
rows.append(dict(bound="B3_pulse_amplitude_cap",
                 statement="miRNA-mediated I1 overshoot ratio is capped at 1 + lambda_y*q(m_ss); "
                           "a transcriptional I1 has no such cap",
                 n_checked=N, median_predicted=float(np.median(cap)),
                 median_observed=float(np.median(over)),
                 max_relative_error=float(np.nanmax(np.maximum(over - cap, 0) / cap)),
                 median_relative_error=float(np.nanmedian(np.maximum(over - cap, 0) / cap))))
rows[-1]["fraction_violating_the_cap"] = float(np.mean(over > cap * 1.001))

# corollary to B1: with a finite shared RISC capacity Kc the pool saturates at Kc
qKc = 1.0 - (P["Kmy"] ** P["nmy"]) / (P["Kmy"] ** P["nmy"] + P["Kc"] ** P["nmy"])
rows.append(dict(bound="B1b_finite_RISC_capacity_corollary",
                 statement="with a finite shared RISC/AGO capacity Kc the effective repressor "
                           "pool saturates at Kc, so the achievable repression falls from "
                           "1+lambda_y to 1+lambda_y*q(Kc): competition for RISC caps how much "
                           "any single miRNA can repress",
                 n_checked=N, median_predicted=float(np.median(1 + P["lam_y"] * qKc)),
                 median_observed=np.nan, max_relative_error=np.nan,
                 median_relative_error=np.nan))

B = pd.DataFrame(rows)
B.to_csv(f"{RES}/dynamics_analytic_bounds.csv", index=False)
print(B.to_string(index=False))
