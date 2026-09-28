#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02c_characterise_3node_N16384.py
================================
Re-run of 02_characterise_3node.py with the ensemble enlarged from 1,024 to
16,384 scrambled-Sobol parameter sets per topology, matching the size used in
the higher-order sweep (03_higher_order_sweep.py) so that Parts B and C of the
dynamical analysis rest on the same sampling density.

Because scipy's scrambled Sobol engine is deterministic given the seed and
random_base2 returns a PREFIX of the sequence, parameter sets 0..1023 of this
run are IDENTICAL to the whole of the previous 1,024-set run.  The script
therefore also serves as a reproducibility check: metrics for those first 1,024
sets must reproduce the earlier file exactly.

Memory is bounded by evaluating the ensemble in chunks of CHUNK columns.

Outputs are written to a staging directory first and only promoted to
results/v3/ once every topology has completed.
"""
import sys, os, time, shutil, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_metrics as M, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
STAGE = "/path/to/scratch/stage_3node"
N_ENS = 16384
CHUNK = 4096
NPROC = 12
os.makedirs(STAGE, exist_ok=True)
t0 = time.time()

reg = D.build_registry()
Pref = D.default_params(1)
Pens = SMP.sample(N_ENS, seed=20260908)
Pall = {k: np.concatenate([np.asarray(Pref[k], float), np.asarray(Pens[k], float)])
        for k in Pref}
NT = N_ENS + 1
S_grid = np.logspace(np.log10(0.02), np.log10(20.0), 41)


def run_one(nm):
    t = reg[nm]
    mets, aux0 = [], None
    for a in range(0, NT, CHUNK):
        b = min(a + CHUNK, NT)
        Pc = {k: v[a:b] for k, v in Pall.items()}
        met, aux = M.characterise(t, Pc, b - a, S_grid=S_grid)
        if a == 0:
            aux0 = aux
        mets.append(pd.DataFrame({k: np.asarray(v, float) for k, v in met.items()}))
    df = pd.concat(mets, ignore_index=True)
    df.insert(0, "param_set", np.arange(NT) - 1)
    df.insert(0, "topology", nm)
    df.to_csv(f"{STAGE}/metrics_{nm}.csv", index=False)

    tc, dr = [], []
    for i, tt in enumerate(aux0["ts"]):
        tc.append(dict(topology=nm, protocol="ON_step", t=tt, protein=aux0["p_on"][i, 0]))
    for i, tt in enumerate(aux0["ts2"]):
        tc.append(dict(topology=nm, protocol="OFF_step", t=tt, protein=aux0["p_off"][i, 0]))
    for i, tt in enumerate(aux0["ts3"]):
        tc.append(dict(topology=nm, protocol="brief_pulse_0.5tau", t=tt, protein=aux0["p_short"][i, 0]))
    for i, S in enumerate(S_grid):
        dr.append(dict(topology=nm, S=S, protein_ss=aux0["Ps"][i, 0], mRNA_ss=aux0["Ys"][i, 0]))
    pd.DataFrame(tc).to_csv(f"{STAGE}/tc_{nm}.csv", index=False)
    pd.DataFrame(dr).to_csv(f"{STAGE}/dr_{nm}.csv", index=False)
    print(f"  {nm:16s} done  ({time.time()-t0:7.1f}s)", flush=True)
    return nm


if __name__ == "__main__":
    names = list(reg)
    with Pool(processes=NPROC) as pool:
        pool.map(run_one, names)

    E = pd.concat([pd.read_csv(f"{STAGE}/metrics_{nm}.csv") for nm in names], ignore_index=True)
    cas = E[E.topology == "CASCADE"].set_index("param_set")
    E = E.merge(cas[["T_half_ON", "T_half_OFF"]].rename(
        columns={"T_half_ON": "T_half_ON_cascade", "T_half_OFF": "T_half_OFF_cascade"}),
        left_on="param_set", right_index=True, how="left")
    E["rel_T_ON"] = E.T_half_ON / E.T_half_ON_cascade
    E["rel_T_OFF"] = E.T_half_OFF / E.T_half_OFF_cascade
    E["sign_sensitive_delay_index"] = np.log2(E.rel_T_ON / E.rel_T_OFF)

    # ---- reproducibility check against the previous 1,024-set run --------------
    old_path = f"{RES}/dynamics_3node_metrics_ensemble.csv"
    if os.path.exists(old_path):
        O = pd.read_csv(old_path)
        key = ["topology", "param_set"]
        cols = [c for c in O.columns if c not in key and c in E.columns]
        Osub = O[O.param_set < 1024].sort_values(key).reset_index(drop=True)
        Esub = E[E.param_set < 1024].sort_values(key).reset_index(drop=True)
        assert (Osub[key].values == Esub[key].values).all(), "key mismatch"
        worst = {}
        for c in cols:
            a = Osub[c].to_numpy(float); b2 = Esub[c].to_numpy(float)
            both = np.isfinite(a) & np.isfinite(b2)
            nanmis = int((np.isfinite(a) != np.isfinite(b2)).sum())
            d = np.abs(a[both] - b2[both]).max() if both.any() else 0.0
            worst[c] = (float(d), nanmis)
        rep = pd.DataFrame([dict(metric=c, max_abs_diff=v[0], n_finiteness_mismatch=v[1])
                            for c, v in worst.items()]).sort_values("max_abs_diff", ascending=False)
        rep.to_csv(f"{RES}/dynamics_3node_reproducibility_vs_N1024.csv", index=False)
        print("\nREPRODUCIBILITY vs previous 1,024-set run (first 1,024 sets):")
        print(rep.head(12).to_string(index=False))

    np.savez_compressed(f"{RES}/dynamics_3node_sobol_params.npz", **Pall)
    E.to_csv(f"{RES}/dynamics_3node_metrics_ensemble.csv", index=False)
    pd.concat([pd.read_csv(f"{STAGE}/tc_{nm}.csv") for nm in names], ignore_index=True) \
      .to_csv(f"{RES}/dynamics_reference_timecourses.csv", index=False)
    pd.concat([pd.read_csv(f"{STAGE}/dr_{nm}.csv") for nm in names], ignore_index=True) \
      .to_csv(f"{RES}/dynamics_reference_doseresponse.csv", index=False)
    print(f"\nDONE in {time.time()-t0:.1f}s   rows={len(E)}")
