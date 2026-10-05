#!/usr/bin/env python3
"""A1 of SETTINGS.md: runs modules of the extended model (dyn_models_bhat.py) with the paper's run_module()
(scripts/11_dynamics/03_higher_order_sweep.py: integration, scoring and thresholds unchanged) and the chunked runner of
analyses/six_node_pattern S2. The 16,384-set design (dyn_sample.py, seed 20260908) is unchanged; Kx21 and nx21 (TF2 -> TF1)
come from a separate 2-dimensional scrambled Sobol sequence (seed 20260909) with the ranges of Kx12 and nx12, paired
by index. Each parameter set is integrated independently (the model is vectorised column-wise), so chunks of sets
run in parallel.

Module list (CSV) columns:
  family, module          labels of the output file
  registry                'C2_MIR:composite', 'I1_MIR:composite' or 'I1_MIR:plain' for a module built from
                          dyn_models.higher_order_registry() (the paper's six-node modules; nested check), else empty
  size_key                n3/n4/n5/n6 (registry modules) or the layers: core, GG, GG+MM, GG+MM+TT
  mir_to_TF, input_node, s_TM, s_TY, s_T12, s_T21, s_T2Y, s_T2M   topology fields (empty = the registry's or default)
usage: python3 a1_run.py <analysis root> <module list csv> <output folder> [workers=10] [chunk=2048]
Output: <output folder>/perset/<family>__<module>.csv.gz
"""
import sys, os, time, importlib.util, dataclasses
sys.dont_write_bytecode = True
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np, pandas as pd
from scipy.stats import qmc
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = sys.argv[1] if len(sys.argv) > 1 else None
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(ROOT, "scripts", "11_dynamics"))
import dyn_models_bhat as D
import dyn_sample as SMP
spec = importlib.util.spec_from_file_location("hos", os.path.join(ROOT, "scripts", "11_dynamics", "03_higher_order_sweep.py"))
hos = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hos)            # defines functions only; its __main__ block does not run
hos.D = D                               # the extended model; with s_T21 = 0 it is the original operation for operation
N_ALL, SEED, SEED21 = hos.N_SET, hos.SEED, 20260909


def extra_params(n):
    """Kx21, nx21: separate scrambled Sobol sequence, ranges of Kx12 (log) and nx12 (linear)."""
    U = qmc.Sobol(d=2, scramble=True, seed=SEED21).random_base2(int(np.ceil(np.log2(n))))[:n]
    lo, hi = SMP.RANGES["Kx12"][:2]; lo2, hi2 = SMP.RANGES["nx12"][:2]
    return {"Kx21": np.exp(np.log(lo) + U[:, 0] * (np.log(hi) - np.log(lo))), "nx21": lo2 + U[:, 1] * (hi2 - lo2)}


class _ChunkSampler:
    """Stands in for dyn_sample inside run_module: one chunk of the full 16,384-set sample, plus Kx21 and nx21."""
    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def sample(self, n, seed):
        full = SMP.sample(N_ALL, seed=seed); full.update(extra_params(N_ALL))
        return {k: np.asarray(v)[self.lo:self.hi] for k, v in full.items()}


LAYERS = {"core": dict(gene_gene="none", mir_mir=False, tf_tf=False), "GG": dict(gene_gene="mutual", mir_mir=False, tf_tf=False),
          "GG+MM": dict(gene_gene="mutual", mir_mir=True, tf_tf=False), "GG+MM+TT": dict(gene_gene="mutual", mir_mir=True, tf_tf=True)}
NN = {"core": 3, "GG": 4, "GG+MM": 5, "GG+MM+TT": 6}


def topo_of(row):
    ov = {k: row[k] for k in ("s_TM", "s_TY", "s_T12", "s_T21", "s_T2Y", "s_T2M") if str(row.get(k, "")) not in ("", "nan")}
    ov = {k: int(float(v)) for k, v in ov.items()}
    if str(row.get("registry", "")) not in ("", "nan"):
        base, kind = row["registry"].split(":")
        t = D.higher_order_registry("AND", base, composite=(kind == "composite"))[row["size_key"]]
        return dataclasses.replace(t, **ov)
    lay = LAYERS[row["size_key"]]
    return D.Topology(name=f"{row['family']}_{row['module']}", label=f"Bhat {row['family']} {row['module']}", n_nodes=NN[row["size_key"]],
                      gate="AND", mir_arm=True, s_MY=-1, mir_to_TF=str(row["mir_to_TF"]).lower() in ("true", "1"),
                      input_node=row["input_node"], **lay, **{"s_TM": 0, "s_TY": 1, "s_T12": 1, "s_T2Y": 1, "s_T2M": 1, "s_T21": 0, **ov})


def run_chunk(args):
    row, lo, hi = args
    hos.SMP = _ChunkSampler(lo, hi)
    hos.N_SET = hi - lo
    t0 = time.time()
    df = hos.run_module((row["family"], row["module"], topo_of(row)))
    df["param_set"] = df["param_set"] + lo
    return row["family"], row["module"], lo, hi, df, time.time() - t0


def main(root, mlist, out, workers=10, chunk=2048):
    workers, chunk = int(workers), int(chunk)
    os.makedirs(os.path.join(out, "perset"), exist_ok=True)
    M = pd.read_csv(mlist, dtype=str).fillna("")
    t0 = time.time(); tasks, need = [], {}
    for row in M.to_dict("records"):
        fn = os.path.join(out, "perset", f"{row['family']}__{row['module'].replace('|', '')}.csv.gz")
        if os.path.exists(fn): print(f"skip (exists): {os.path.basename(fn)}", flush=True); continue
        t = topo_of(row)
        print(f"{row['family']:24s} {row['module']:34s} n={t.n_nodes} input={t.input_node} mir_to_TF={t.mir_to_TF!s:5s} s_TM={t.s_TM:+d} "
              f"s_TY={t.s_TY:+d} GG={t.gene_gene} MM={t.mir_mir!s:5s} TT={t.tf_tf!s:5s} s_T12={t.s_T12:+d} s_T21={t.s_T21:+d} "
              f"s_T2Y={t.s_T2Y:+d} s_T2M={t.s_T2M:+d}", flush=True)
        need[(row["family"], row["module"])] = []
        for lo in range(0, N_ALL, chunk): tasks.append((row, lo, min(lo + chunk, N_ALL)))
    nch = {k: sum(1 for t in tasks if (t[0]["family"], t[0]["module"]) == k) for k in need}
    print(f"{len(tasks)} chunk tasks, chunk={chunk}, workers={workers}", flush=True)
    with Pool(processes=workers) as pool:
        for fam, mod, lo, hi, df, dt in pool.imap_unordered(run_chunk, tasks, chunksize=1):
            need[(fam, mod)].append(df)
            print(f"  [{(time.time() - t0) / 60:6.1f} min] {fam} {mod} sets {lo}-{hi - 1} in {dt:.0f}s ({len(need[(fam, mod)])}/{nch[(fam, mod)]})", flush=True)
            if len(need[(fam, mod)]) == nch[(fam, mod)]:
                A = pd.concat(need[(fam, mod)], ignore_index=True).sort_values("param_set")
                fn = os.path.join(out, "perset", f"{fam}__{mod.replace('|', '')}.csv.gz")
                A.to_csv(fn + ".tmp", index=False, compression="gzip"); os.replace(fn + ".tmp", fn)
                print(f"  WROTE {os.path.basename(fn)} ({len(A)} sets)", flush=True); need[(fam, mod)] = []
    print(f"ALL DONE in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:6])
