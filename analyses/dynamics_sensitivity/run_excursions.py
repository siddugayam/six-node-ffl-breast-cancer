#!/usr/bin/env python3
"""One-at-a-time parameter excursions of the higher-order module sweep.

For every row of excursion_plan.csv (module, parameter, factor 0.5 or 2), all 16,384 sampled values of that parameter
are multiplied by the factor and run_module() of scripts/11_dynamics/03_higher_order_sweep.py is run unchanged, on
chunks of 2,048 points (the sampler is replaced by one that returns the scaled chunk, as
analyses/six_node_pattern/S2/s2_run_modules.py does). Before any excursion, an unscaled run of points 0-2,047 of
COMP_C2_toggle n6 must reproduce the stored flags of those points exactly (gate).
Writes excursion_runs/<family>__<module>__<parameter>__x<factor>.csv.gz (param_set and the 10 binary flags) and
excursion_runs/gate.txt; runs whose file exists are skipped, so the script can be restarted.
usage: FFL_ROOT=<repository checkout> OMP_NUM_THREADS=1 python3 run_excursions.py [workers, default 10]"""
import importlib.util, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "excursion_runs")
DYN = os.path.join(ROOT, "scripts", "11_dynamics")
sys.path.insert(0, DYN)
import dyn_models as D, dyn_sample as SMP
spec = importlib.util.spec_from_file_location("higher_order_sweep", os.path.join(DYN, "03_higher_order_sweep.py"))
hos = importlib.util.module_from_spec(spec); spec.loader.exec_module(hos)

N_ALL, SEED, CHUNK = 16384, 20260908, 2048
FAMILIES = {"I1_miRNA_FFL": ("I1_MIR", False), "COMP_C2_toggle": ("C2_MIR", True),
            "COMP_I1_negfeedback": ("I1_MIR", True)}
BIN = ["is_memory", "is_pulse", "is_sustained_osc", "is_damped_osc", "is_bistable", "bistability_unresolved",
       "is_ultrasensitive", "is_very_ultrasensitive", "is_noise_rejecting", "nonmonotone_dose"]


class ScaledChunkSampler:
    """Stands in for dyn_sample inside run_module: one chunk of the full sample, one parameter multiplied by a factor."""
    def __init__(self, lo, hi, param, factor):
        self.lo, self.hi, self.param, self.factor = lo, hi, param, factor

    def sample(self, n, seed):
        full = SMP.sample(N_ALL, seed=seed)
        out = {k: np.asarray(v, float)[self.lo:self.hi].copy() for k, v in full.items()}
        if self.param is not None:
            out[self.param] = out[self.param] * self.factor
        return out


def run_chunk(args):
    fam, mod, param, factor, lo, hi = args
    base, comp = FAMILIES[fam]
    topo = D.higher_order_registry("AND", base, composite=comp)[mod]
    hos.SMP = ScaledChunkSampler(lo, hi, param, factor)
    hos.N_SET = hi - lo
    t0 = time.time()
    df = hos.run_module((fam, mod, topo))
    df["param_set"] = df["param_set"] + lo
    return args, df[["param_set"] + BIN], time.time() - t0


def fname(fam, mod, param, factor):
    return os.path.join(OUT, f"{fam}__{mod}__{param}__x{factor:g}.csv.gz")


if __name__ == "__main__":
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    gate_file = os.path.join(OUT, "gate.txt")
    if not os.path.exists(gate_file):
        _, df, dt = run_chunk(("COMP_C2_toggle", "n6", None, 1.0, 0, CHUNK))
        st = pd.read_csv(os.path.join(ROOT, "results", "systems_and_dynamics", "dynamics_higher_order_persetset.csv.gz"))
        st = st[(st.family == "COMP_C2_toggle") & (st.module == "n6") & (st.param_set < CHUNK)].sort_values("param_set")
        diff = {b: int((df[b].astype(float).values != st[b].astype(float).values).sum()) for b in BIN}
        ok = all(v == 0 for v in diff.values())
        with open(gate_file + (".tmp" if ok else ".failed"), "w") as fh:
            fh.write(f"gate: unscaled COMP_C2_toggle n6, points 0-{CHUNK - 1}, {dt:.0f} s\n")
            fh.write("differing points per flag: " + ", ".join(f"{b} {v}" for b, v in diff.items()) + "\n")
            fh.write("PASS\n" if ok else "FAIL\n")
        if not ok:
            print("gate FAILED", diff, flush=True); sys.exit(2)
        os.replace(gate_file + ".tmp", gate_file)
        print(f"gate passed in {dt:.0f} s", flush=True)
    plan = pd.read_csv(os.path.join(HERE, "excursion_plan.csv")).sort_values("priority")
    tasks, need, nch = [], {}, {}
    for r in plan.itertuples():
        key = (r.family, r.module, r.parameter, float(r.factor))
        if os.path.exists(fname(*key)):
            continue
        need[key] = []
        for lo in range(0, N_ALL, CHUNK):
            tasks.append(key + (lo, lo + CHUNK))
        nch[key] = N_ALL // CHUNK
    print(f"{len(need)} excursion runs to do, {len(tasks)} chunk tasks, {workers} workers", flush=True)
    with Pool(processes=workers) as pool:
        for args, df, dt in pool.imap_unordered(run_chunk, tasks, chunksize=1):
            key = args[:4]
            need[key].append(df)
            print(f"[{(time.time() - t0) / 60:7.1f} min] {key[0]} {key[1]} {key[2]} x{key[3]:g} points "
                  f"{args[4]}-{args[5] - 1} in {dt:.0f} s ({len(need[key])}/{nch[key]})", flush=True)
            if len(need[key]) == nch[key]:
                A = pd.concat(need[key], ignore_index=True).sort_values("param_set")
                assert (A.param_set.values == np.arange(N_ALL)).all()
                fn = fname(*key)
                A.to_csv(fn + ".tmp", index=False, compression="gzip"); os.replace(fn + ".tmp", fn)
                need[key] = []
                print(f"  written {os.path.basename(fn)}", flush=True)
    print(f"ALL_DONE in {(time.time() - t0) / 3600:.2f} h", flush=True)
