#!/usr/bin/env python3
"""A1 nested check, certification part (SETTINGS.md, Changes and decisions, 2026-09-28 20:05): the sets that the nested
run (extended model, TF2 -> TF1 absent) flags is_sustained_osc are certified exactly as a1_certify.py does it (1,000 tau
from two initial conditions), and compared set by set with the paper's certification:
  COMP_C2_toggle      results/v3/dynamics_n6_limit_cycle_certified.csv (03f; 83 of 108)
  COMP_I1_negfeedback analyses/six_node_pattern S2, s2_limit_cycle_certification.csv, module GG+MM+TT (31 of 37)
usage: python3 a1_nested_certify.py <analysis root> <nested module list csv> <nested perset folder> <S2 certification csv>
                                    <output csv> [workers]
Exit status 0 only if the flagged sets and every set's certification equal the stored ones."""
import sys, os
sys.dont_write_bytecode = True
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import a1_certify as C                                  # certify(): the extended model, 1,000 tau, two initial conditions


def stored_certification(root, s2cert):
    st = {"COMP_C2_toggle": pd.read_csv(os.path.join(root, "results/v3/dynamics_n6_limit_cycle_certified.csv"))}
    S2 = pd.read_csv(s2cert)
    st["COMP_I1_negfeedback"] = S2[(S2.family == "COMP_I1_negfeedback") & (S2.module == "GG+MM+TT")]
    for k, v in st.items():
        assert v.certified_limit_cycle.dtype == bool, (k, v.certified_limit_cycle.dtype)
        st[k] = v[["param_set", "certified_limit_cycle"]].astype({"param_set": int})
    return st


def main(root, mlist, perset, s2cert, outf, workers=3):
    st = stored_certification(root, s2cert)
    M = pd.read_csv(mlist, dtype=str).fillna("")
    tasks = []
    for row in M.to_dict("records"):
        X = pd.read_csv(os.path.join(perset, f"{row['family']}__{row['module'].replace('|', '')}.csv.gz"), usecols=["param_set", "is_sustained_osc"])
        fl = np.sort(X.loc[X.is_sustained_osc > 0.5, "param_set"].to_numpy(int))
        for k in range(0, len(fl), 64): tasks.append((row, fl[k:k + 64]))
    with Pool(int(workers)) as pool:
        out = pd.concat(list(pool.imap_unordered(C.certify, tasks)), ignore_index=True).sort_values(["family", "param_set"])
    out.to_csv(outf, index=False)
    L, ok = [], True
    for fam, g in out.groupby("family"):
        m = g[["param_set", "certified_limit_cycle"]].merge(st[fam], on="param_set", how="outer", suffixes=("", "_stored"), indicator=True)
        same = bool((m._merge == "both").all())
        agree = same and bool((m.certified_limit_cycle.astype(bool) == m.certified_limit_cycle_stored.astype(bool)).all())
        ok &= agree
        L.append(f"NESTED CERTIFICATION {fam} GG+MM+TT: {int(g.certified_limit_cycle.sum())} of {len(g)} flagged sets certified here vs "
                 f"{int(st[fam].certified_limit_cycle.sum())} of {len(st[fam])} stored; same flagged sets {same}; "
                 f"set-by-set agreement {agree} -> {'PASS' if agree else 'FAIL'}")
        print(L[-1], flush=True)
    open(os.path.splitext(outf)[0] + ".log", "w").write("\n".join(L) + "\n")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(*sys.argv[1:7])
