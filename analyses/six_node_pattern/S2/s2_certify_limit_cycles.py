#!/usr/bin/env python3
"""
Certification of sustained oscillation for every S2/S3 module, exactly as
scripts/11_dynamics/03f_n6_limit_cycle_certification.py does it for COMP_C2_toggle n6:
every set flagged is_sustained_osc by the sweep is re-integrated for 1,000 tau (DT 0.01, input S = 2.0,
target protein recorded every 20 steps) from two initial conditions (all active states 0, and all 1.5);
certified = relative peak-to-peak amplitude over the last 25 % > 2 % from both, the two amplitudes within
5 % and the two periods within 5 % (amp_period() copied verbatim from 03f).

Validation first: the stored COMP_C2_toggle n6 certification (results/v3/dynamics_n6_limit_cycle_certified.csv,
83 of 108) is recomputed here and must agree set by set.

usage: s2_certify_limit_cycles.py [validate]    (validate = only the stored n6 check)
Output: s2_limit_cycle_certification.csv (one row per flagged set per module); log to stdout.
"""
import sys, os
sys.dont_write_bytecode = True
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s2_common as C
import s2_run_modules as RM          # topologies (topo_of) and the dyn_models / dyn_sample modules
D_, SMP = RM.D, RM.SMP
N_SET, SEED, DT, S_OP, T_LONG, REC = 16384, 20260908, 0.01, 2.0, 1000.0, 20


def amp_period(ts, y, frac=0.25):                     # verbatim from 03f
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


def certify(topo, sets):
    P0 = SMP.sample(N_SET, seed=SEED)
    P = {k: np.asarray(v, float)[sets] for k, v in P0.items()}
    N = len(sets); act = topo.active_states(); rhs = D_.make_rhs(topo, P)
    out = {}
    for ic_name, ic in (("zero", 0.0), ("high", 1.5)):
        Z0 = np.zeros((D_.NSTATE, N)); Z0[act] = ic
        Zf, ts, tr = D_.integrate(rhs, Z0, lambda t: S_OP, T_LONG, DT, record_every=REC, active=act,
                                  record_states=(D_.IP1,))
        out[ic_name] = amp_period(ts, tr[:, 0, :])
    rows = []
    for j, ps in enumerate(sets):
        a0, a1 = out["zero"][0][j], out["high"][0][j]; p0, p1 = out["zero"][1][j], out["high"][1][j]
        ok_amp = (a0 > 0.02) and (a1 > 0.02)
        ok_conv = (abs(a0 - a1) / max(a0, a1) < 0.05) if (a0 > 0 and a1 > 0) else False
        ok_per = (np.isfinite(p0) and np.isfinite(p1) and abs(p0 - p1) / max(p0, p1) < 0.05)
        rows.append(dict(param_set=int(ps), rel_amp_zero_IC=a0, rel_amp_high_IC=a1, period_tau_zero_IC=p0,
                         period_tau_high_IC=p1, amplitude_sustained_1000tau=ok_amp,
                         independent_of_initial_condition=bool(ok_conv and ok_per),
                         certified_limit_cycle=bool(ok_amp and ok_conv and ok_per)))
    return pd.DataFrame(rows)


def topology(f, mod):
    """S2/S3 topologies from s2_run_modules; the S6 family from S6/s6_run_modules.py."""
    if f == "COMP_I3_allneg_TT":
        import importlib.util
        spec = importlib.util.spec_from_file_location("s6", os.path.join(HERE, "..", "S6", "s6_run_modules.py"))
        s6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s6)
        return s6.TOP[mod]
    return RM.topo_of(f, mod)


if __name__ == "__main__":
    stored = pd.read_csv(f"{C.REV}/results/v3/dynamics_n6_limit_cycle_certified.csv")
    sets = np.sort(stored.param_set.to_numpy(int))
    v = certify(RM.topo_of("COMP_C2_toggle", "GG+MM+TT"), sets)
    m = v.merge(stored[["param_set", "certified_limit_cycle", "n6_rel_amp_zero_IC", "n6_rel_amp_high_IC"]], on="param_set")
    agree = (m.certified_limit_cycle_x == m.certified_limit_cycle_y.astype(bool)).all()
    dmax = max((m.rel_amp_zero_IC - m.n6_rel_amp_zero_IC).abs().max(), (m.rel_amp_high_IC - m.n6_rel_amp_high_IC).abs().max())
    print(f"VALIDATION stored n6: {int(v.certified_limit_cycle.sum())} of {len(v)} certified here vs "
          f"{int(stored.certified_limit_cycle.sum())} stored; set-by-set agreement {agree}; max |amplitude diff| {dmax:.2e}", flush=True)
    assert agree
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        sys.exit(0)
    D = C.load_all()
    res = []
    OUTF = os.path.join(HERE, "s2_limit_cycle_certification.csv")
    prev = pd.read_csv(OUTF) if os.path.exists(OUTF) else None          # incremental: keep modules already certified
    done = set() if prev is None else set(map(tuple, prev.loc[prev.source == "s2_certify_limit_cycles.py", ["family", "module"]].drop_duplicates().values))
    if prev is not None:
        res.append(prev[prev.source == "s2_certify_limit_cycles.py"].drop(columns=["source"]))
    for (f, mod), g in D.groupby(["family", "module"]):
        fl = np.sort(g.loc[g.is_sustained_osc > 0.5, "param_set"].to_numpy(int))
        if len(fl) == 0 or (f, mod) == ("COMP_C2_toggle", "GG+MM+TT") or (f, mod) in done:
            if (f, mod) in done: print(f"{f:22s} {mod:26s} already certified (kept)", flush=True)
            continue
        r = certify(topology(f, mod), fl); r.insert(0, "module", mod); r.insert(0, "family", f)
        res.append(r)
        print(f"{f:22s} {mod:26s} flagged {len(fl):4d}  certified {int(r.certified_limit_cycle.sum()):4d}", flush=True)
    stored_rows = stored.rename(columns={"n6_rel_amp_zero_IC": "rel_amp_zero_IC", "n6_rel_amp_high_IC": "rel_amp_high_IC",
                                         "n6_period_tau": "period_tau_zero_IC"})
    stored_rows = stored_rows[["param_set", "rel_amp_zero_IC", "rel_amp_high_IC", "period_tau_zero_IC",
                               "amplitude_sustained_1000tau", "independent_of_initial_condition", "certified_limit_cycle"]].copy()
    stored_rows.insert(0, "module", "GG+MM+TT"); stored_rows.insert(0, "family", "COMP_C2_toggle")
    stored_rows["source"] = "results/v3/dynamics_n6_limit_cycle_certified.csv (03f, stored)"
    out = pd.concat([stored_rows] + [r.assign(source="s2_certify_limit_cycles.py") for r in res], ignore_index=True)
    out.to_csv(OUTF + ".tmp", index=False); os.replace(OUTF + ".tmp", OUTF)
    print("WROTE s2_limit_cycle_certification.csv", flush=True)
