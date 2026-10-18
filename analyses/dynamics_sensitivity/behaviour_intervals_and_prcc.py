#!/usr/bin/env python3
"""Behaviour frequencies of the higher-order module sweep (Table 4, Fig. 5): exact Clopper-Pearson intervals, PRCC of
every acting parameter against memory, pulse and oscillation, the parameters chosen for the one-at-a-time excursions,
and the parameter table.

Reads the public repository (root: environment variable FFL_ROOT):
  results/systems_and_dynamics/dynamics_higher_order_persetset.csv.gz, ..._persetset_compI1.csv.gz (per-point flags),
  dynamics_higher_order_gain_loss*.csv, dynamics_higher_order_fractions*.csv (stored summaries, for checks),
  dynamics_n6_limit_cycle_certified.csv, analyses/six_node_pattern/S2/s2_limit_cycle_certification.csv,
  dynamics_parameter_definitions.csv, dynamics_3node_sobol_params.npz,
  scripts/11_dynamics/dyn_sample.py and dyn_models.py (imported unchanged).
Writes (this folder): behaviour_frequencies_clopper_pearson.csv, parameters_acting_by_module.csv, prcc_binary.csv,
prcc_continuous.csv, prcc_top5.csv, excursion_plan.csv, parameter_table.csv,
integration_settings.csv (constants read from 03_higher_order_sweep.py).
usage: FFL_ROOT=<repository checkout> python3 behaviour_intervals_and_prcc.py"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta, rankdata, t as tdist

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results", "systems_and_dynamics")
sys.path.insert(0, os.path.join(ROOT, "scripts", "11_dynamics"))
import dyn_models as D, dyn_sample as SMP

N_SET, SEED = 16384, 20260908
FAMILIES = {"I1_miRNA_FFL": ("I1_MIR", False), "COMP_C2_toggle": ("C2_MIR", True),
            "COMP_I1_negfeedback": ("I1_MIR", True)}
MODULES = ["n3", "n4", "n5", "n6", "n4d"]
BIN = ["is_memory", "is_pulse", "is_sustained_osc", "is_damped_osc", "is_bistable", "bistability_unresolved",
       "is_ultrasensitive", "is_very_ultrasensitive", "is_noise_rejecting", "nonmonotone_dose"]
PRIMARY = {"memory": "is_memory", "pulse": "is_pulse", "oscillation": "is_sustained_osc"}
PROXY = {"memory": "T_half_OFF", "pulse": "overshoot_ratio", "oscillation": "rel_peak_to_peak"}
MIN_EACH = 10


def cp(k, n):
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return lo, hi


# ---- the parameter sample, checked against the stored three-node sample ------------------------------------------
P = {k: np.asarray(v, float) for k, v in SMP.sample(N_SET, seed=SEED).items()}
npz = np.load(os.path.join(RES, "dynamics_3node_sobol_params.npz"))
for k in SMP.NAMES:
    assert np.array_equal(npz[k][1:1025], P[k][:1024]), f"sample check failed for {k}"
print("parameter sample: points 0-1023 equal the stored npz (all 35 parameters)")

# ---- per-point flags ------------------------------------------------------------------------------------------------
A = pd.concat([pd.read_csv(os.path.join(RES, "dynamics_higher_order_persetset.csv.gz")),
               pd.read_csv(os.path.join(RES, "dynamics_higher_order_persetset_compI1.csv.gz"))], ignore_index=True)
GL = pd.concat([pd.read_csv(os.path.join(RES, "dynamics_higher_order_gain_loss.csv")),
                pd.read_csv(os.path.join(RES, "dynamics_higher_order_gain_loss_compI1.csv"))], ignore_index=True)
FR = pd.concat([pd.read_csv(os.path.join(RES, "dynamics_higher_order_fractions.csv")),
                pd.read_csv(os.path.join(RES, "dynamics_higher_order_fractions_compI1.csv"))], ignore_index=True)

# ---- (a) Clopper-Pearson --------------------------------------------------------------------------------------------
rows = []
for fam in FAMILIES:
    for mod in MODULES:
        g = A[(A.family == fam) & (A.module == mod)].sort_values("param_set")
        assert len(g) == N_SET and (g.param_set.values == np.arange(N_SET)).all(), (fam, mod)
        fr = FR[(FR.family == fam) & (FR.module == mod)].iloc[0]
        for b in BIN:
            v = g[b].astype(float).values
            assert np.isfinite(v).all()
            k = int((v > 0.5).sum()); n = N_SET
            assert abs(100 * k / n - fr[f"pct_{b}"]) < 1e-9, (fam, mod, b)
            if mod == "n3":
                chk = GL[(GL.family == fam) & (GL.module == "n4") & (GL.behaviour == b)].iloc[0]
                assert k == int(chk.n_loss + chk.n_both)
            else:
                chk = GL[(GL.family == fam) & (GL.module == mod) & (GL.behaviour == b)].iloc[0]
                assert k == int(chk.n_gain + chk.n_both)
            lo, hi = cp(k, n)
            rows.append(dict(family=fam, module=mod, behaviour=b, k=k, n=n, pct=100 * k / n,
                             cp95_lo_pct=100 * lo, cp95_hi_pct=100 * hi, flag="uncertified (as stored)"))
cert_c2 = pd.read_csv(os.path.join(RES, "dynamics_n6_limit_cycle_certified.csv"))
k_c2 = int(cert_c2.certified_limit_cycle.astype(str).str.upper().eq("TRUE").sum())
s2 = pd.read_csv(os.path.join(ROOT, "analyses", "six_node_pattern", "S2", "s2_limit_cycle_certification.csv"))
cert = lambda fam: int(s2[(s2.family == fam) & (s2.module == "GG+MM+TT")].certified_limit_cycle.astype(str)
                       .str.upper().eq("TRUE").sum())
assert cert("COMP_C2_toggle") == k_c2, "certified counts disagree"
for fam, mod, k in [("COMP_C2_toggle", "n6", k_c2), ("COMP_I1_negfeedback", "n6", cert("COMP_I1_negfeedback"))]:
    lo, hi = cp(k, N_SET)
    rows.append(dict(family=fam, module=mod, behaviour="is_sustained_osc", k=k, n=N_SET, pct=100 * k / N_SET,
                     cp95_lo_pct=100 * lo, cp95_hi_pct=100 * hi, flag="certified limit cycle"))
CPDF = pd.DataFrame(rows)

# ---- acting parameters ----------------------------------------------------------------------------------------------
rng = np.random.default_rng(1)
pts = rng.choice(N_SET, 16, replace=False)
topos = {(fam, mod): D.higher_order_registry("AND", base, composite=comp)[mod]
         for fam, (base, comp) in FAMILIES.items() for mod in MODULES}
Z = rng.uniform(0.0, 2.0, size=(D.NSTATE, 16 * 64))
Prep = {k: np.repeat(v[pts], 64) for k, v in P.items()}
act_rows, ACT = [], {}
for (fam, mod), topo in topos.items():
    act = topo.active_states()
    Za = np.zeros_like(Z); Za[act] = Z[act]          # inactive states stay 0, as in every integration
    acting = []
    for S in (0.10, 2.0):
        f0 = D.make_rhs(topo, Prep)(0.0, Za, S)[act]
        for k in SMP.NAMES:
            P2 = dict(Prep); P2[k] = Prep[k] * 2.0
            f1 = D.make_rhs(topo, P2)(0.0, Za, S)[act]
            if np.any(np.abs(f1 - f0) > 1e-12 * (1.0 + np.abs(f0))) and k not in acting: acting.append(k)
    acting = [k for k in SMP.NAMES if k in acting]
    ACT[(fam, mod)] = acting
    act_rows.append(dict(family=fam, module=mod, n_acting=len(acting), acting=";".join(acting),
                         not_acting=";".join(k for k in SMP.NAMES if k not in acting)))
pd.DataFrame(act_rows).to_csv(os.path.join(HERE, "parameters_acting_by_module.csv"), index=False)


# ---- (b) PRCC -------------------------------------------------------------------------------------------------------
def prcc(X, y):
    """X: (N, p) parameter values, y: (N,) output; returns PRCC and p per column."""
    R = np.column_stack([rankdata(X[:, j]) for j in range(X.shape[1])])
    ry = rankdata(y)
    N, p = R.shape
    out_r, out_p = np.empty(p), np.empty(p)
    for j in range(p):
        O = np.column_stack([np.ones(N), np.delete(R, j, axis=1)])
        bx, *_ = np.linalg.lstsq(O, R[:, j], rcond=None)
        by, *_ = np.linalg.lstsq(O, ry, rcond=None)
        ex, ey = R[:, j] - O @ bx, ry - O @ by
        r = float(np.corrcoef(ex, ey)[0, 1])
        dfree = N - 2 - (p - 1)
        tt = r * np.sqrt(dfree / max(1e-300, 1 - r * r))
        out_r[j], out_p[j] = r, 2 * tdist.sf(abs(tt), dfree)
    return out_r, out_p


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


pr_bin, pr_con, top = [], [], []
for fam in FAMILIES:
    for mod in MODULES:
        g = A[(A.family == fam) & (A.module == mod)].sort_values("param_set")
        names = ACT[(fam, mod)]
        X = np.column_stack([P[k] for k in names])
        for beh, col in PRIMARY.items():
            y = g[col].astype(float).values
            k1 = int((y > 0.5).sum())
            if k1 < MIN_EACH or N_SET - k1 < MIN_EACH:
                pr_bin.append(dict(family=fam, module=mod, behaviour=beh, flag=col, k=k1, parameter="",
                                   prcc=np.nan, p=np.nan, q_bh=np.nan, rank_abs=np.nan,
                                   note=f"not estimable ({k1} of {N_SET} points with the behaviour)"))
                continue
            r, p = prcc(X, y); q = bh(p)
            order = np.argsort(-np.abs(r)); rk = np.empty(len(r), int); rk[order] = np.arange(1, len(r) + 1)
            for j, k in enumerate(names):
                pr_bin.append(dict(family=fam, module=mod, behaviour=beh, flag=col, k=k1, parameter=k, prcc=r[j],
                                   p=p[j], q_bh=q[j], rank_abs=rk[j], note=""))
            for j in order[:5]:
                top.append(dict(family=fam, module=mod, behaviour=beh, k=k1, rank=int(rk[j]), parameter=names[j],
                                prcc=r[j], p=p[j], q_bh=q[j]))
            # secondary: continuous proxy
            pc = PROXY[beh]
            z = g[pc].astype(float).values
            if pc == "T_half_OFF":
                z = np.where(np.isfinite(z), z, np.inf); ok = np.ones(N_SET, bool)
                z = np.where(np.isinf(z), np.nanmax(z[np.isfinite(z)]) + 1.0, z)
            else:
                ok = np.isfinite(z)
            r2, p2 = prcc(X[ok], z[ok]); q2 = bh(p2)
            o2 = np.argsort(-np.abs(r2)); rk2 = np.empty(len(r2), int); rk2[o2] = np.arange(1, len(r2) + 1)
            for j, k in enumerate(names):
                pr_con.append(dict(family=fam, module=mod, behaviour=beh, proxy=pc, n_points=int(ok.sum()),
                                   parameter=k, prcc=r2[j], p=p2[j], q_bh=q2[j], rank_abs=rk2[j]))
PB, PC, TOP = pd.DataFrame(pr_bin), pd.DataFrame(pr_con), pd.DataFrame(top)
PB.to_csv(os.path.join(HERE, "prcc_binary.csv"), index=False)
PC.to_csv(os.path.join(HERE, "prcc_continuous.csv"), index=False)
TOP.to_csv(os.path.join(HERE, "prcc_top5.csv"), index=False)
CPDF.to_csv(os.path.join(HERE, "behaviour_frequencies_clopper_pearson.csv"), index=False)

# ---- excursion plan (c): union of the top-5 over the estimable behaviours, modules n3-n6 -----------------------------
ORDER = [("n6", f) for f in ("COMP_C2_toggle", "COMP_I1_negfeedback", "I1_miRNA_FFL")] + \
        [("n5", f) for f in ("COMP_C2_toggle", "COMP_I1_negfeedback", "I1_miRNA_FFL")] + \
        [("n3", f) for f in ("COMP_C2_toggle", "COMP_I1_negfeedback", "I1_miRNA_FFL")] + \
        [("n4", f) for f in ("COMP_C2_toggle", "COMP_I1_negfeedback", "I1_miRNA_FFL")]
plan, prio = [], 0
for mod, fam in ORDER:
    sel = TOP[(TOP.family == fam) & (TOP.module == mod)]
    params = [k for k in ACT[(fam, mod)] if k in set(sel.parameter)]
    for k in params:
        for f in (0.5, 2.0):
            prio += 1
            plan.append(dict(priority=prio, family=fam, module=mod, parameter=k, factor=f,
                             chosen_for=";".join(sorted(set(sel[sel.parameter == k].behaviour)))))
pd.DataFrame(plan).to_csv(os.path.join(HERE, "excursion_plan.csv"), index=False)

# ---- (d) parameter table --------------------------------------------------------------------------------------------
defs = pd.read_csv(os.path.join(RES, "dynamics_parameter_definitions.csv")).set_index("symbol")
tab = []
for k, (lo, hi, sc, why) in SMP.RANGES.items():
    used = [f"{fam}:{mod}" for (fam, mod), a in ACT.items() if k in a]
    tab.append(dict(parameter=k, quantity=defs.loc[k, "quantity"] if k in defs.index else "",
                    meaning=defs.loc[k, "meaning"] if k in defs.index else "", low=lo, high=hi,
                    distribution="log-uniform" if sc == "log" else "uniform", justification=why,
                    acts_in_modules=";".join(used) if used else "none of the modules (OR-gate parameter)",
                    sampler="scrambled Sobol (scipy.stats.qmc.Sobol, d=35, scramble=True, seed=20260908), "
                            "first 16,384 points (random_base2(14)); u in [0,1) mapped as lo*(hi/lo)^u (log) or "
                            "lo+u*(hi-lo) (uniform)"))
pd.DataFrame(tab).to_csv(os.path.join(HERE, "parameter_table.csv"), index=False)
import importlib.util
spec = importlib.util.spec_from_file_location("higher_order_sweep", os.path.join(ROOT, "scripts", "11_dynamics", "03_higher_order_sweep.py"))
hos = importlib.util.module_from_spec(spec); spec.loader.exec_module(hos)
integ = [("integrator", "fixed-step fourth-order Runge-Kutta (dyn_models.integrate), states clipped at 0"), ("dt", hos.DT),
         ("parameter sets", hos.N_SET), ("Sobol seed", hos.SEED), ("operating input S", hos.S_OP),
         ("low input", 0.10), ("settling time (tau)", hos.T_SETTLE), ("ON and OFF steps (tau)", 50.0),
         ("short input pulse (tau)", 0.5), ("dose-response levels", ";".join(f"{x:.4g}" for x in hos.S_DOSE)),
         ("time per dose level (tau)", hos.T_DOSE), ("oscillation run (tau)", hos.T_OSC),
         ("bistability input levels", ";".join(f"{x:g}" for x in hos.S_BIS)), ("bistability residual tolerance", hos.BIS_TOL),
         ("bistability settling chunk (tau)", hos.BIS_CHUNK), ("bistability maximum time (tau)", hos.BIS_TMAX)]
pd.DataFrame(integ, columns=["setting", "value"]).to_csv(os.path.join(HERE, "integration_settings.csv"), index=False)
print(CPDF[CPDF.behaviour.isin(["is_memory", "is_pulse", "is_sustained_osc"])].to_string(index=False))
print(TOP.to_string(index=False))
print(f"excursion plan: {len(plan)} runs")
