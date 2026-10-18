#!/usr/bin/env python3
"""Behaviour frequencies of the one-at-a-time parameter excursions against the published sweep.

For every row of excursion_plan.csv: the frequency of each binary behaviour over the 16,384 points with the parameter
multiplied by the factor (excursion_runs/), its exact 95 % Clopper-Pearson interval, the frequency at the sampled values
(the stored per-point flags of results/systems_and_dynamics/dynamics_higher_order_persetset*.csv.gz), and the matched
change point by point (gained, lost, exact McNemar p). Runs without a file are reported as NOT RUN.
Writes excursion_frequencies.csv, and excursion_range_by_module.csv: per module and behaviour, the baseline frequency
and the lowest and highest excursion frequency with the excursion that gave each (descriptive only, no test).
usage: FFL_ROOT=<repository checkout> python3 summarise_excursions.py"""
import os
import numpy as np, pandas as pd
from scipy.stats import beta, binomtest

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results", "systems_and_dynamics")
N = 16384
BIN = ["is_memory", "is_pulse", "is_sustained_osc", "is_damped_osc", "is_bistable", "bistability_unresolved",
       "is_ultrasensitive", "is_very_ultrasensitive", "is_noise_rejecting", "nonmonotone_dose"]


def cp(k, n):
    return (0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1), 1.0 if k == n else beta.ppf(0.975, k + 1, n - k))


A = pd.concat([pd.read_csv(os.path.join(RES, "dynamics_higher_order_persetset.csv.gz")),
               pd.read_csv(os.path.join(RES, "dynamics_higher_order_persetset_compI1.csv.gz"))], ignore_index=True)
plan = pd.read_csv(os.path.join(HERE, "excursion_plan.csv")).sort_values("priority")
out = []
for r in plan.itertuples():
    base = A[(A.family == r.family) & (A.module == r.module)].sort_values("param_set")
    fn = os.path.join(HERE, "excursion_runs", f"{r.family}__{r.module}__{r.parameter}__x{r.factor:g}.csv.gz")
    ex = pd.read_csv(fn).sort_values("param_set") if os.path.exists(fn) else None
    if ex is not None: assert (ex.param_set.values == base.param_set.values).all() and len(ex) == N
    for b in BIN:
        a = base[b].astype(float).values > 0.5; kb = int(a.sum())
        row = dict(priority=r.priority, family=r.family, module=r.module, parameter=r.parameter, factor=r.factor,
                   chosen_for=r.chosen_for, behaviour=b, n=N, baseline_k=kb, baseline_pct=100 * kb / N)
        if ex is None:
            row.update(status="NOT RUN"); out.append(row); continue
        c = ex[b].astype(float).values > 0.5; k = int(c.sum()); lo, hi = cp(k, N)
        g = int((c & ~a).sum()); l = int((~c & a).sum())
        row.update(status="run", k=k, pct=100 * k / N, cp95_lo_pct=100 * lo, cp95_hi_pct=100 * hi, gained=g, lost=l,
                   mcnemar_p=binomtest(g, g + l, 0.5).pvalue if g + l > 0 else 1.0)
        out.append(row)
O = pd.DataFrame(out)
O.to_csv(os.path.join(HERE, "excursion_frequencies.csv"), index=False)
R = O[O.status == "run"].copy(); R["change_pp"] = R.pct - R.baseline_pct; rng = []
for (fam, mod, b), d in R.groupby(["family", "module", "behaviour"], sort=False):
    lo, hi = d.loc[d.pct.idxmin()], d.loc[d.pct.idxmax()]
    rng.append(dict(family=fam, module=mod, behaviour=b, excursions=len(d), baseline_pct=d.baseline_pct.iloc[0],
                    min_pct=lo.pct, min_at=f"{lo.parameter} x{lo.factor:g}", max_pct=hi.pct, max_at=f"{hi.parameter} x{hi.factor:g}",
                    median_abs_change_pp=d.change_pp.abs().median(), max_abs_change_pp=d.change_pp.abs().max(),
                    appears=bool(d.baseline_k.iloc[0] == 0 and (d.k > 0).any()), vanishes=bool(d.baseline_k.iloc[0] > 0 and (d.k == 0).any())))
pd.DataFrame(rng).to_csv(os.path.join(HERE, "excursion_range_by_module.csv"), index=False)
done = O[O.status == "run"][["family", "module", "parameter", "factor"]].drop_duplicates()
print(f"{len(done)} of {len(plan)} excursion runs present")
