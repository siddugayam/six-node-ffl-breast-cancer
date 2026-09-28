#!/usr/bin/env python3
"""
S6 (analyses/six_node_pattern): model what the network contains.  Condition (brief): fully signed six-node
configurations not covered by S2-S3 make up >= 10 % of the fully signed MODEL6-architecture instances.
S5(a) (S5/s5a_top10_fully_signed_configurations.csv): 16 MODEL6-architecture instances are fully signed;
12 (75 %) carry one configuration not covered by S2-S3, 4 carry the S3a configuration (TF2 -| miRNA):
    TF1 -| miR1, TF1 -| miR2, TF1 -| G1, TF1 -| G2, TF1 -| TF2, TF2 -| G1, TF2 -| G2, TF2 -| miR1, TF2 -| miR2
i.e. the composite COMP_I3 core of dyn_models.py (s_TM = -1, s_TY = -1, miRNA -| TF1) with every TF2 arc
negative (s_T12 = s_T2Y = s_T2M = -1).  It is the only uncovered fully signed configuration, so it is the
only one run ("the three most frequent (or fewer)").
Modules: the factorial of S2 on that core; the analysis plan's minimum (core, TT, GG+TT, GG+MM+TT) first, then the
other four.  Everything else is exactly the S2 pipeline (S2/s2_run_modules.py: the same 16,384 Sobol sets,
seed 20260908, run_module() of 03_higher_order_sweep.py, chunked).
BHAT6 configurations are not run: their compulsory TF arcs are T1 -> T2 AND T2 -> T1 and T2 -> G2 (Bhat Fig.
1d), which have no counterpart in dyn_models.py (its TF-TF layer is TF1 -> TF2 only) without changing the
model's structure.
usage: s6_run_modules.py [chunk=2048] [workers=3]
(chunks are checkpointed in S6/chunks/ with claim files, so a second process can join the work)
"""
import sys, os, time, dataclasses
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'S2'))
import s2_run_modules as RM
import pandas as pd
from multiprocessing import Pool

FAMILY = "COMP_I3_allneg_TT"
base = RM.FAM["COMP_C2_toggle"]                     # composite C2 registry; re-signed below
r = dataclasses.replace
core = r(base["n3"], name="COMP_I3_MIR_n3", s_TM=-1, s_TY=-1)
signs = dict(s_TM=-1, s_TY=-1, s_T12=-1, s_T2Y=-1, s_T2M=-1)
TOP = {
    "core": core,
    "TT": r(base["n6"], name="COMP_I3neg_TT", n_nodes=4, gene_gene="none", mir_mir=False, **signs),
    "GG+TT": r(base["n6"], name="COMP_I3neg_GG+TT", n_nodes=5, mir_mir=False, **signs),
    "GG+MM+TT": r(base["n6"], name="COMP_I3neg_n6", **signs),
    "GG": r(base["n4"], name="COMP_I3_MIR_n4", s_TM=-1, s_TY=-1),
    "MM": r(base["n5"], name="COMP_I3_MIR_MM", n_nodes=4, gene_gene="none", s_TM=-1, s_TY=-1),
    "GG+MM": r(base["n5"], name="COMP_I3_MIR_n5", s_TM=-1, s_TY=-1),
    "MM+TT": r(base["n6"], name="COMP_I3neg_MM+TT", n_nodes=5, gene_gene="none", **signs),
}
QUEUE = ["core", "TT", "GG+TT", "GG+MM+TT", "GG", "MM", "GG+MM", "MM+TT"]


CHUNKDIR = os.path.join(HERE, "chunks")


def chunk_path(mod, lo):
    return os.path.join(CHUNKDIR, f"{mod}__{lo:05d}.csv.gz")


def run_chunk_ckpt(args):
    """claim-file protocol so that several processes can share the work; each finished chunk is saved at once."""
    mod, lo, hi = args
    fn = chunk_path(mod, lo)
    if os.path.exists(fn): return mod, lo, hi, "exists", 0.0
    try:
        fd = os.open(fn + ".claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.write(fd, str(os.getpid()).encode()); os.close(fd)
    except FileExistsError:
        return mod, lo, hi, "claimed elsewhere", 0.0
    m, lo_, hi_, df, dt = run_chunk(args)
    df.to_csv(f"{fn}.{os.getpid()}.tmp", index=False, compression="gzip"); os.replace(f"{fn}.{os.getpid()}.tmp", fn)
    return mod, lo, hi, "done", dt


def assemble():
    for mod in QUEUE:
        out = os.path.join(HERE, "perset", f"{FAMILY}__{mod}.csv.gz")
        fs = [chunk_path(mod, lo) for lo in range(0, RM.N_ALL, CHUNK)]
        if os.path.exists(out) or not all(os.path.exists(x) for x in fs): continue
        A = pd.concat([pd.read_csv(x) for x in fs], ignore_index=True).sort_values("param_set")
        assert len(A) == RM.N_ALL and (A.param_set.values == range(RM.N_ALL)).all()
        tmp = f"{out}.{os.getpid()}.tmp"                 # per-process name: two processes may assemble at once
        A.to_csv(tmp, index=False, compression="gzip"); os.replace(tmp, out)
        print(f"  WROTE {out} ({len(A)} sets)", flush=True)


def run_chunk(args):
    mod, lo, hi = args
    RM.hos.SMP = RM._ChunkSampler(lo, hi); RM.hos.N_SET = hi - lo
    t0 = time.time()
    df = RM.hos.run_module((FAMILY, mod, TOP[mod])); df["param_set"] = df["param_set"] + lo
    return mod, lo, hi, df, time.time() - t0


if __name__ == "__main__":
    CHUNK = int(sys.argv[1]) if len(sys.argv) > 1 else 2048
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    os.makedirs(os.path.join(HERE, "perset"), exist_ok=True); os.makedirs(CHUNKDIR, exist_ok=True)
    t0 = time.time(); tasks = []
    for mod in QUEUE:
        t = TOP[mod]
        print(f"{FAMILY} {mod:9s} {t.name:20s} gg={t.gene_gene:6s} mm={t.mir_mir!s:5s} tt={t.tf_tf!s:5s} s_TM={t.s_TM:+d} "
              f"s_TY={t.s_TY:+d} s_T12={t.s_T12:+d} s_T2Y={t.s_T2Y:+d} s_T2M={t.s_T2M:+d} mir_to_TF={t.mir_to_TF}", flush=True)
        for lo in range(0, RM.N_ALL, CHUNK):
            if not os.path.exists(chunk_path(mod, lo)): tasks.append((mod, lo, min(lo + CHUNK, RM.N_ALL)))
    print(f"{len(tasks)} chunk tasks pending, workers={workers}", flush=True)
    with Pool(processes=workers) as pool:
        for mod, lo, hi, st, dt in pool.imap_unordered(run_chunk_ckpt, tasks, chunksize=1):
            if st == "done": print(f"  [{(time.time()-t0)/60:6.1f} min] {mod} sets {lo}-{hi-1} in {dt:.0f}s", flush=True)
            assemble()
    assemble()
    print(f"ALL DONE in {(time.time()-t0)/60:.1f} min", flush=True)
