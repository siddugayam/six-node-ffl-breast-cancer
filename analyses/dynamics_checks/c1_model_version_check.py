#!/usr/bin/env python3
"""
C1 check: which version of dyn_models.py produced results/v3/dynamics_higher_order_persetset.csv.gz?

03_higher_order_sweep.py started at 2026-09-09 08:19:12; dyn_models.py was last modified at 08:27:30
(the 08:04 backup dyn_models.py.bak has the unbounded mutual-stabilisation term 1/(1+kappa*p2), the
current file the saturating 1/(1+kappa*f+(p2))).  The n3 core does not use that term; n4, n5 and n6 do.
This script re-runs the first 64 parameter sets of the stored composite n4 and n6 modules with the
CURRENT code, and also with the current code switched to MUTUAL_STABILISATION = "linear_legacy",
and reports which reproduces the stored per-set values.  Reads the project read-only.
"""
import sys, os, time, importlib.util
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
from multiprocessing import Pool

REV = "/path/to/revision"
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{REV}/scripts/11_dynamics")
import dyn_models as D, dyn_sample as SMP

spec = importlib.util.spec_from_file_location("hos", f"{REV}/scripts/11_dynamics/03_higher_order_sweep.py")
hos = importlib.util.module_from_spec(spec); spec.loader.exec_module(hos)
REG = D.higher_order_registry("AND", "C2_MIR", composite=True)
NSUB = 64


class _Sub:
    def sample(self, n, seed):
        full = SMP.sample(16384, seed=seed)          # the full design, then its first NSUB sets
        return {k: np.asarray(v)[:NSUB] for k, v in full.items()}


def run(args):
    module, variant = args
    D.MUTUAL_STABILISATION = variant
    hos.SMP = _Sub(); hos.N_SET = NSUB
    t0 = time.time()
    df = hos.run_module(("COMP_C2_toggle", module, REG[module]))
    df["variant"] = variant
    return df, time.time() - t0


if __name__ == "__main__":
    stored = pd.read_csv(f"{REV}/results/v3/dynamics_higher_order_persetset.csv.gz")
    tasks = [(m, v) for m in ("n4", "n6") for v in ("saturating", "linear_legacy")]
    with Pool(4) as pool:
        res = pool.map(run, tasks)
    rows = []
    for (module, variant), (df, secs) in zip(tasks, res):
        old = stored[(stored.family == "COMP_C2_toggle") & (stored.module == module)].set_index("param_set").loc[:NSUB - 1]
        new = df.set_index("param_set")
        num = [c for c in new.columns if c not in ("family", "module", "variant")]
        a = old[num].to_numpy(float); b = new[num].to_numpy(float)
        same = np.isclose(a, b, rtol=1e-9, atol=1e-12, equal_nan=True)
        rel = np.nanmax(np.abs(a - b) / np.maximum(np.abs(a), 1e-12))
        flags = ["is_bistable", "is_sustained_osc", "is_damped_osc", "is_ultrasensitive", "is_memory", "is_pulse"]
        nflag = int((old[flags].to_numpy(float) != new[flags].to_numpy(float)).sum())
        rows.append(dict(module=module, model_variant=variant, n_sets=NSUB, seconds=round(secs),
                         cells_matching=int(same.sum()), cells_total=int(same.size),
                         max_rel_diff=float(rel), behaviour_flags_differing=nflag,
                         p_hi_first5_stored=list(np.round(old["p_hi"].to_numpy()[:5], 6)),
                         p_hi_first5_rerun=list(np.round(new["p_hi"].to_numpy()[:5], 6))))
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/c1_model_version_check.csv", index=False)
    print(R.to_string(index=False))
