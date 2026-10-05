#!/usr/bin/env python3
"""A1 of SETTINGS.md: certification of sustained oscillation for the A1 modules, exactly as analyses/six_node_pattern
S2/s2_certify_limit_cycles.py does it (after scripts/11_dynamics/03f_n6_limit_cycle_certification.py): every set flagged
is_sustained_osc by the sweep is re-integrated for 1,000 tau (DT 0.01, input S = 2.0, target protein recorded every 20
steps) from two initial conditions (all active states 0, and all 1.5); certified = relative peak-to-peak amplitude over
the last 25 % > 2 % from both, the two amplitudes within 5 % and the two periods within 5 % (amp_period() verbatim).
The extended model (dyn_models_bhat.py) and the A1 parameter sample (including Kx21, nx21) are used.
usage: python3 a1_certify.py <analysis root> <module list csv> <perset folder> <output csv> [workers]
"""
import sys, os
sys.dont_write_bytecode = True
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import a1_run as A                                      # the extended model, topologies and parameter sample
D_ = A.D
DT, S_OP, T_LONG, REC = 0.01, 2.0, 1000.0, 20


def amp_period(ts, y, frac=0.25):                     # verbatim from 03f (via S2)
    k = int((1 - frac) * len(ts))
    t = ts[k:]; Y = y[k:]
    pk = Y.max(axis=0) - Y.min(axis=0)
    mean = np.maximum(np.abs(Y.mean(axis=0)), 1e-9)
    rel = pk / mean
    Cc = Y - Y.mean(axis=0)[None, :]
    per = np.full(Y.shape[1], np.nan)
    for j in range(Y.shape[1]):
        s = np.sign(Cc[:, j])
        idx = np.where(np.diff(s) != 0)[0]
        if len(idx) >= 3:
            per[j] = 2.0 * np.median(np.diff(t[idx]))
    return rel, per


def certify(args):
    row, sets = args
    topo = A.topo_of(row)
    P0 = A._ChunkSampler(0, A.N_ALL).sample(A.N_ALL, A.SEED)
    P = {k: np.asarray(v, float)[sets] for k, v in P0.items()}
    N = len(sets); act = topo.active_states(); rhs = D_.make_rhs(topo, P)
    out = {}
    for ic_name, ic in (("zero", 0.0), ("high", 1.5)):
        Z0 = np.zeros((D_.NSTATE, N)); Z0[act] = ic
        Zf, ts, tr = D_.integrate(rhs, Z0, lambda t: S_OP, T_LONG, DT, record_every=REC, active=act, record_states=(D_.IP1,))
        out[ic_name] = amp_period(ts, tr[:, 0, :])
    rows = []
    for j, ps in enumerate(sets):
        a0, a1 = out["zero"][0][j], out["high"][0][j]; p0, p1 = out["zero"][1][j], out["high"][1][j]
        ok_amp = (a0 > 0.02) and (a1 > 0.02)
        ok_conv = (abs(a0 - a1) / max(a0, a1) < 0.05) if (a0 > 0 and a1 > 0) else False
        ok_per = (np.isfinite(p0) and np.isfinite(p1) and abs(p0 - p1) / max(p0, p1) < 0.05)
        rows.append(dict(family=row["family"], module=row["module"], param_set=int(ps), rel_amp_zero_IC=a0, rel_amp_high_IC=a1,
                         period_tau_zero_IC=p0, period_tau_high_IC=p1, amplitude_sustained_1000tau=ok_amp,
                         independent_of_initial_condition=bool(ok_conv and ok_per), certified_limit_cycle=bool(ok_amp and ok_conv and ok_per)))
    return pd.DataFrame(rows)


def main(root, mlist, perset, outf, workers=8):
    # validation first: the stored COMP_C2_toggle six-node certification (83 of 108), with the extended code and TF2 -> TF1 absent
    stored = pd.read_csv(os.path.join(root, "results/v3/dynamics_n6_limit_cycle_certified.csv"))
    vrow = dict(family="COMP_C2_toggle", module="GG+MM+TT@validate", registry="C2_MIR:composite", size_key="n6", s_T21="0")
    for k in ("mir_to_TF", "input_node", "s_TM", "s_TY", "s_T12", "s_T2Y", "s_T2M"): vrow[k] = ""
    v = certify((vrow, np.sort(stored.param_set.to_numpy(int))))
    m = v.merge(stored[["param_set", "certified_limit_cycle"]], on="param_set")
    agree = bool((m.certified_limit_cycle_x == m.certified_limit_cycle_y.astype(bool)).all())
    vline = (f"VALIDATION stored n6 certification: {int(v.certified_limit_cycle.sum())} of {len(v)} certified here vs "
             f"{int(stored.certified_limit_cycle.sum())} stored; set-by-set agreement {agree}")
    print(vline, flush=True)
    LOG = [vline]
    assert agree
    M = pd.read_csv(mlist, dtype=str).fillna("")
    tasks = []
    for row in M.to_dict("records"):
        fn = os.path.join(perset, f"{row['family']}__{row['module'].replace('|', '')}.csv.gz")
        X = pd.read_csv(fn, usecols=["param_set", "is_sustained_osc"])
        fl = np.sort(X.loc[X.is_sustained_osc > 0.5, "param_set"].to_numpy(int))
        print(f"{row['family']:22s} {row['module']:30s} flagged {len(fl)}", flush=True)
        LOG.append(f"{row['family']} {row['module']}: flagged sustained oscillation {len(fl)}")
        for k in range(0, len(fl), 64): tasks.append((row, fl[k:k + 64]))
    res = []
    with Pool(int(workers)) as pool:
        for r in pool.imap_unordered(certify, tasks): res.append(r)
    out = pd.concat(res, ignore_index=True) if res else pd.DataFrame(columns=["family", "module", "param_set", "certified_limit_cycle"])
    out = out.sort_values(["family", "module", "param_set"])
    out.to_csv(outf, index=False)
    if len(out):
        for (f, m), g in out.groupby(["family", "module"]): LOG.append(f"{f} {m}: certified {int(g.certified_limit_cycle.sum())} of {len(g)} flagged")
    open(os.path.splitext(outf)[0] + ".log", "w").write("\n".join(LOG) + "\n")
    print("\n".join(LOG[-5:]), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:6])
