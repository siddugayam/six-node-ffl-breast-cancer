#!/usr/bin/env python3
"""
C3 control for the first re-runs.

Question: do the six-node composite behaviours need the whole six-node module, or only the TF-TF
layer?  Builds a four-node composite module = composite C2 core + the second TF only
(TF1 -> TF2, TF2 -> target G1, TF2 -> miRNA1; no gene-gene and no miRNA-miRNA layer), taken from
the six-node topology of dyn_models.higher_order_registry("AND", "C2_MIR", composite=True)["n6"]
with gene_gene="none" and mir_mir=False and nothing else changed.

Everything else is reused unchanged and read-only from scripts/v3:
  dyn_models.py, dyn_sample.py (same scrambled-Sobol sample: 16,384 sets, seed 20260908) and
  run_module() of 03_higher_order_sweep.py (same integration, same scoring, same thresholds).
Each parameter set is integrated independently (the model is vectorised column-wise), so the
16,384 sets are split into chunks run in parallel; the three-node core is re-run the same way and
compared set-by-set with the stored results/v3/dynamics_higher_order_persetset.csv.gz as a
reproducibility check.

usage: c3_tf2_only_control.py test   (first 64 sets of the core, compared with the stored run)
       c3_tf2_only_control.py full   (all 16,384 sets, core re-run + four-node TF2 module)
"""
import sys, os, time, importlib.util, dataclasses
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
from multiprocessing import Pool

REV = "/path/to/revision"
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{REV}/scripts/v3")
import dyn_models as D, dyn_sample as SMP

spec = importlib.util.spec_from_file_location("hos", f"{REV}/scripts/v3/03_higher_order_sweep.py")
hos = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hos)          # defines functions only; its __main__ block does not run

N_ALL, SEED = hos.N_SET, hos.SEED
FAMILY = "COMP_C2_toggle"
REG = D.higher_order_registry("AND", "C2_MIR", composite=True)
TOPO_N4TF = dataclasses.replace(REG["n6"], name="COMP_C2_MIR_n4tf",
                                label="4-node: composite core + TF2 only (TF1->TF2, TF2->G1, TF2->miR1)",
                                n_nodes=4, gene_gene="none", mir_mir=False)
BIN = ["is_bistable", "bistability_unresolved", "is_sustained_osc", "is_damped_osc", "is_ultrasensitive",
       "is_very_ultrasensitive", "is_memory", "is_noise_rejecting", "is_pulse", "nonmonotone_dose"]


class _ChunkSampler:
    """Stands in for dyn_sample inside run_module: returns the chunk of the full 16,384-set sample."""
    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def sample(self, n, seed):
        full = SMP.sample(N_ALL, seed=seed)
        return {k: np.asarray(v)[self.lo:self.hi] for k, v in full.items()}


def run_chunk(args):
    module, topo, lo, hi = args
    hos.SMP = _ChunkSampler(lo, hi)
    hos.N_SET = hi - lo
    t0 = time.time()
    df = hos.run_module((FAMILY, module, topo))
    df["param_set"] = df["param_set"] + lo
    print(f"  {module} sets {lo}-{hi-1} done in {time.time()-t0:.0f}s", flush=True)
    return df


def compare_with_stored(new_core):
    A = pd.read_csv(f"{REV}/results/v3/dynamics_higher_order_persetset.csv.gz")
    old = A[(A.family == FAMILY) & (A.module == "n3")].set_index("param_set").loc[new_core.param_set]
    new = new_core.set_index("param_set")
    rows = []
    for c in new.columns:
        if c in ("family", "module"):
            continue
        a = old[c].to_numpy(float); b = new[c].to_numpy(float)
        both_nan = np.isnan(a) & np.isnan(b)
        diff = np.where(both_nan, 0.0, np.abs(a - b))
        rows.append(dict(column=c, n_sets=len(a), n_differing=int((diff > 1e-9).sum()),
                         max_abs_diff=float(np.nanmax(diff)) if len(diff) else 0.0))
    return pd.DataFrame(rows)


def gain_loss(base, sub, label_core, label_mod):
    rows = []
    b = base.set_index("param_set").sort_index(); s = sub.set_index("param_set").sort_index()
    assert (b.index == s.index).all()
    from scipy.stats import binomtest
    for beh in BIN:
        a = b[beh].astype(float).values; c = s[beh].astype(float).values
        gain = int(((c > 0.5) & (a < 0.5)).sum()); loss = int(((c < 0.5) & (a > 0.5)).sum())
        both = int(((c > 0.5) & (a > 0.5)).sum()); nn = len(a)
        p = binomtest(gain, gain + loss, 0.5).pvalue if (gain + loss) > 0 else 1.0
        rows.append(dict(core=label_core, module=label_mod, behaviour=beh, n=nn,
                         pct_core=100 * np.nanmean(a), pct_module=100 * np.nanmean(c),
                         n_gain=gain, n_loss=loss, n_both=both,
                         pct_gain=100 * gain / nn, pct_loss=100 * loss / nn,
                         net_pct=100 * (gain - loss) / nn, mcnemar_p=p))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    mode = sys.argv[1]
    t0 = time.time()
    if mode == "test":
        core = run_chunk(("n3", REG["n3"], 0, 64))
        cmp_ = compare_with_stored(core)
        cmp_.to_csv(f"{OUT}/c3_test_core64_vs_stored.csv", index=False)
        print(cmp_.to_string(index=False))
        print(f"test done in {time.time()-t0:.0f}s")
        sys.exit(0)

    chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    tasks = []
    for module, topo in (("n4tf", TOPO_N4TF), ("n3", REG["n3"])):
        for lo in range(0, N_ALL, chunk):
            tasks.append((module, topo, lo, min(lo + chunk, N_ALL)))
    print(f"topology n4tf: {TOPO_N4TF}", flush=True)
    print(f"active states n4tf: {[D.STATE_NAMES[i] for i in TOPO_N4TF.active_states()]}", flush=True)
    with Pool(processes=workers) as pool:
        parts = pool.map(run_chunk, tasks, chunksize=1)
    A = pd.concat(parts, ignore_index=True).sort_values(["module", "param_set"])
    A.to_csv(f"{OUT}/c3_persetset.csv.gz", index=False, compression="gzip")
    core_new = A[A.module == "n3"].reset_index(drop=True)
    mod = A[A.module == "n4tf"].reset_index(drop=True)

    cmp_ = compare_with_stored(core_new)
    cmp_.to_csv(f"{OUT}/c3_core_rerun_vs_stored.csv", index=False)

    S = pd.read_csv(f"{REV}/results/v3/dynamics_higher_order_persetset.csv.gz")
    core_stored = S[(S.family == FAMILY) & (S.module == "n3")].reset_index(drop=True)
    gl = pd.concat([gain_loss(core_stored, mod, "n3 stored (09-09 sweep)", "n4tf (this run)"),
                    gain_loss(core_new, mod, "n3 re-run (this run)", "n4tf (this run)")], ignore_index=True)
    gl.to_csv(f"{OUT}/c3_gain_loss_vs_core.csv", index=False)

    frac = []
    for (m), g in A.groupby("module"):
        r = dict(family=FAMILY, module=m, n_param_sets=len(g))
        for b in BIN:
            r[f"pct_{b}"] = 100 * np.nanmean(g[b].astype(float))
        frac.append(r)
    pd.DataFrame(frac).to_csv(f"{OUT}/c3_fractions.csv", index=False)
    print(cmp_.to_string(index=False))
    print(gl.to_string(index=False))
    print(f"\nDONE in {(time.time()-t0)/60:.1f} min")
