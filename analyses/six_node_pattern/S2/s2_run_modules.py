#!/usr/bin/env python3
"""
S2 (layer-factorial dynamics) and S3 (TF-TF sign controls) for analyses/six_node_pattern.

Runs every module whose per-set flags are NOT already stored.  Everything is reused unchanged and
read-only from scripts/11_dynamics: dyn_models.py, dyn_sample.py (the same scrambled-Sobol design: 16,384 sets,
seed 20260908) and run_module() of 03_higher_order_sweep.py (same integration, scoring and thresholds).
Each parameter set is integrated independently (the model is vectorised column-wise), so the 16,384
sets are split into chunks run in parallel, as in the C3 control of the first re-runs.

Factorial modules = the three-node core plus any combination of the three layers of
dyn_models.higher_order_registry():
  GG  gene_gene="mutual"   (G1 <-> G2 mutual protein stabilisation, G2 co-regulated)
  MM  mir_mir=True         (miRNA2 co-transcribed with miRNA1, shared RISC pool)
  TT  tf_tf=True           (TF1 -> TF2, TF2 -> G1, TF2 -> miRNA1)
Each is built from the stored n3/n4/n5/n6 topology with layer switches only (dataclasses.replace);
nothing else changes.  Stored per-set flags reused (not re-run):
  core=n3, GG=n4, GG+MM=n5, GG+MM+TT=n6  results/v3/dynamics_higher_order_persetset.csv.gz
                                          (I1_miRNA_FFL, COMP_C2_toggle) and ..._compI1.csv.gz
                                          (COMP_I1_negfeedback)
  COMP_C2_toggle TT                       analyses/dynamics_controls/c3_persetset.csv.gz (n4tf)
The stored COMP_I1_negfeedback modules are verified first by re-running sets 0-2047 of n3 and n6.

S3 variants (COMP_C2_toggle six-node module, one sign changed each):
  S3a TF2 -| miRNA (s_T2M=-1); S3b TF2 -| G1 (s_T2Y=-1); S3c TF1 -| TF2 (s_T12=-1); S3c on the TT-only module.

usage: s2_run_modules.py [chunk=2048] [workers=10]
Output: perset/<family>__<module>.csv.gz, written as soon as every chunk of that module is done.
"""
import sys, os, time, importlib.util, dataclasses
sys.dont_write_bytecode = True
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np, pandas as pd
from multiprocessing import Pool

REV = "/path/to/revision"
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{REV}/scripts/11_dynamics")
import dyn_models as D, dyn_sample as SMP

spec = importlib.util.spec_from_file_location("hos", f"{REV}/scripts/11_dynamics/03_higher_order_sweep.py")
hos = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hos)          # defines functions only; its __main__ block does not run
N_ALL, SEED = hos.N_SET, hos.SEED

FAM = {"COMP_C2_toggle": D.higher_order_registry("AND", "C2_MIR", composite=True),
       "COMP_I1_negfeedback": D.higher_order_registry("AND", "I1_MIR", composite=True),
       "I1_miRNA_FFL": D.higher_order_registry("AND", "I1_MIR", composite=False)}


def factorial(reg, tag, mod):
    """module label -> Topology, built from the stored topologies by layer switches only."""
    r = dataclasses.replace
    return {
        "core": reg["n3"], "GG": reg["n4"], "GG+MM": reg["n5"], "GG+MM+TT": reg["n6"],
        "MM": r(reg["n5"], name=f"{tag}_MM", label="core + miRNA-miRNA only", n_nodes=4, gene_gene="none"),
        "TT": r(reg["n6"], name=f"{tag}_TT", label="core + TF2 only", n_nodes=4, gene_gene="none", mir_mir=False),
        "GG+TT": r(reg["n6"], name=f"{tag}_GG+TT", label="core + gene-gene + TF2", n_nodes=5, mir_mir=False),
        "MM+TT": r(reg["n6"], name=f"{tag}_MM+TT", label="core + miRNA-miRNA + TF2", n_nodes=5, gene_gene="none"),
    }[mod]


def s3(variant):
    reg = FAM["COMP_C2_toggle"]; r = dataclasses.replace
    tt = factorial(reg, "COMP_C2_MIR", "TT")
    return {"GG+MM+TT_S3a_TF2-|miR": r(reg["n6"], name="COMP_C2_MIR_n6_S3a", s_T2M=-1),
            "GG+MM+TT_S3b_TF2-|G1": r(reg["n6"], name="COMP_C2_MIR_n6_S3b", s_T2Y=-1),
            "GG+MM+TT_S3c_TF1-|TF2": r(reg["n6"], name="COMP_C2_MIR_n6_S3c", s_T12=-1),
            "TT_S3c_TF1-|TF2": r(tt, name="COMP_C2_MIR_TT_S3c", s_T12=-1)}[variant]


TAG = {"COMP_C2_toggle": "COMP_C2_MIR", "COMP_I1_negfeedback": "COMP_I1_MIR", "I1_miRNA_FFL": "I1_MIR"}
# priority order of the analysis plan: (1) C2 GG+TT, MM+TT, MM; (2) COMP_I1 core, TT, GG, GG+TT, GG+MM+TT
# (core, GG, GG+MM+TT are stored: verified here on sets 0-2047); (3) I1 TT, GG+TT; S3 (limit 1.5 h);
# (4) the rest.
QUEUE = [("COMP_C2_toggle", "GG+TT"), ("COMP_C2_toggle", "MM+TT"), ("COMP_C2_toggle", "MM"),
         ("COMP_I1_negfeedback", "core@verify"), ("COMP_I1_negfeedback", "TT"),
         ("COMP_I1_negfeedback", "GG+TT"), ("COMP_I1_negfeedback", "GG+MM+TT@verify"),
         ("I1_miRNA_FFL", "TT"), ("I1_miRNA_FFL", "GG+TT"),
         ("COMP_C2_toggle", "GG+MM+TT_S3c_TF1-|TF2"), ("COMP_C2_toggle", "GG+MM+TT_S3a_TF2-|miR"),
         ("COMP_C2_toggle", "GG+MM+TT_S3b_TF2-|G1"),
         ("COMP_I1_negfeedback", "MM"), ("COMP_I1_negfeedback", "MM+TT"),
         ("I1_miRNA_FFL", "MM"), ("I1_miRNA_FFL", "MM+TT"),
         ("COMP_C2_toggle", "TT_S3c_TF1-|TF2")]


def topo_of(fam, mod):
    base = mod.split("@")[0]
    return s3(base) if "_S3" in base else factorial(FAM[fam], TAG[fam], base)


class _ChunkSampler:
    """Stands in for dyn_sample inside run_module: returns one chunk of the full 16,384-set sample."""
    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def sample(self, n, seed):
        full = SMP.sample(N_ALL, seed=seed)
        return {k: np.asarray(v)[self.lo:self.hi] for k, v in full.items()}


def run_chunk(args):
    fam, mod, lo, hi = args
    hos.SMP = _ChunkSampler(lo, hi)
    hos.N_SET = hi - lo
    t0 = time.time()
    df = hos.run_module((fam, mod, topo_of(fam, mod)))
    df["param_set"] = df["param_set"] + lo
    return fam, mod, lo, hi, df, time.time() - t0


if __name__ == "__main__":
    chunk = int(sys.argv[1]) if len(sys.argv) > 1 else 2048
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    t0 = time.time()
    tasks, need = [], {}
    for fam, mod in QUEUE:
        fn = f"{OUT}/perset/{fam}__{mod.replace('|', '')}.csv.gz"
        if os.path.exists(fn):
            print(f"skip (exists): {fn}", flush=True); continue
        hi_all = 2048 if mod.endswith("@verify") else N_ALL
        t = topo_of(fam, mod)
        print(f"{fam:22s} {mod:26s} {t.name:22s} gene_gene={t.gene_gene:6s} mir_mir={t.mir_mir!s:5s} "
              f"tf_tf={t.tf_tf!s:5s} s_T12={t.s_T12:+d} s_T2Y={t.s_T2Y:+d} s_T2M={t.s_T2M:+d} "
              f"active={[D.STATE_NAMES[i] for i in t.active_states()]}", flush=True)
        need[(fam, mod)] = []
        for lo in range(0, hi_all, chunk):
            tasks.append((fam, mod, lo, min(lo + chunk, hi_all)))
    print(f"{len(tasks)} chunk tasks, chunk={chunk}, workers={workers}", flush=True)
    nchunks = {k: sum(1 for t in tasks if (t[0], t[1]) == k) for k in need}
    with Pool(processes=workers) as pool:
        for fam, mod, lo, hi, df, dt in pool.imap_unordered(run_chunk, tasks, chunksize=1):
            need[(fam, mod)].append(df)
            print(f"  [{(time.time()-t0)/60:6.1f} min] {fam} {mod} sets {lo}-{hi-1} in {dt:.0f}s "
                  f"({len(need[(fam, mod)])}/{nchunks[(fam, mod)]})", flush=True)
            if len(need[(fam, mod)]) == nchunks[(fam, mod)]:
                A = pd.concat(need[(fam, mod)], ignore_index=True).sort_values("param_set")
                fn = f"{OUT}/perset/{fam}__{mod.replace('|', '')}.csv.gz"
                A.to_csv(fn + ".tmp", index=False, compression="gzip"); os.replace(fn + ".tmp", fn)
                print(f"  WROTE {fn} ({len(A)} sets)", flush=True)
                need[(fam, mod)] = []
    print(f"ALL DONE in {(time.time()-t0)/60:.1f} min", flush=True)
