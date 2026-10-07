#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02g_fcd_recomputed.py
=====================
FOLD-CHANGE-DETECTION TEST WITH THE AMPLITUDE MEASURED AT THE PEAK.

In dyn_metrics.characterise, a circuit is scored as fold-change detecting when
its normalised trajectory is invariant across three absolute input levels AND
the response amplitude exceeds 10% of baseline, with the amplitude measured
at the LAST time point of the observation window.  Fold-change detection in
the incoherent type-1 FFL comes with (near-)perfect adaptation, so a
fold-change-detecting circuit returns to its baseline by the end of the
window, where its amplitude is small; this script measures it at the peak
instead (Goentoro, Shoval, Kirschner
& Alon 2009 Mol Cell, doi:10.1016/j.molcel.2009.11.018).

Here responsiveness is measured as the PEAK excursion of the normalised
trajectory instead, which is the quantity Goentoro et al. actually plot, and the
test is repeated at three fold-changes (F = 2, 3, 5).  Everything else is
unchanged.  Both the old and the new call are reported per parameter set so the
difference is auditable.

Out: results/v3/dynamics_fcd_recomputed.csv          (per parameter set)
     results/v3/dynamics_fcd_recomputed_summary.csv  (per topology and fold)
"""
import sys, os, time, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_metrics as M

RES = "/path/to/revision/results/v3"
BASES = (0.25, 0.5, 1.0)
FOLDS = (2.0, 3.0, 5.0)
T_STEP = 30.0
CHUNK = 4096
NPROC = 6

Z = np.load(f"{RES}/dynamics_3node_sobol_params.npz")
Pall = {k: Z[k] for k in Z.files}
NT = len(Pall["gx"])
reg = D.build_registry()


def fcd_one(topo, P, N):
    out = {}
    for F in FOLDS:
        curves = []
        for S0 in BASES:
            Z0 = M.equilibrate(topo, P, S0, N)
            b0 = Z0[D.IP1].copy()
            _, ts, tr = M.step_response(topo, P, N, Z0, F * S0, T_STEP,
                                        states=(D.IY1, D.IP1))
            curves.append(tr[:, 1, :] / np.maximum(b0[None, :], 1e-9))
        C = np.array(curves)                       # (n_base, n_t, N)
        span = C.max(axis=0) - C.min(axis=0)
        scale = np.maximum(np.abs(C.mean(axis=0)), 1e-9)
        err = (span / scale).max(axis=0)
        mean_curve = C.mean(axis=0)
        peak_exc = np.abs(mean_curve - 1.0).max(axis=0)     # NEW responsiveness
        end_exc = np.abs(mean_curve[-1] - 1.0)              # OLD responsiveness
        out[F] = dict(FCD_error=err, peak_excursion=peak_exc, end_excursion=end_exc,
                      is_FCD_peak=((err < 0.10) & (peak_exc > 0.10)).astype(float),
                      is_FCD_end=((err < 0.10) & (end_exc > 0.10)).astype(float))
    return out


def run_one(nm):
    t0 = time.time()
    topo = reg[nm]
    frames = []
    for a in range(0, NT, CHUNK):
        b = min(a + CHUNK, NT)
        P = {k: v[a:b] for k, v in Pall.items()}
        res = fcd_one(topo, P, b - a)
        for F, d in res.items():
            frames.append(pd.DataFrame(dict(topology=nm, fold=F,
                                            param_set=np.arange(a, b) - 1, **d)))
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(f"/path/to/scratch/"
              f"891c54c3-ee6b-4b4f-8631-1f39abee7f27/scratchpad/fcd_{nm}.csv", index=False)
    print(f"  {nm:16s} ({time.time()-t0:6.1f}s)", flush=True)
    return nm


if __name__ == "__main__":
    t0 = time.time()
    names = list(reg)
    with Pool(processes=NPROC) as pool:
        pool.map(run_one, names)
    A = pd.concat([pd.read_csv(f"/path/to/scratch/"
                               f"891c54c3-ee6b-4b4f-8631-1f39abee7f27/scratchpad/fcd_{n}.csv")
                   for n in names], ignore_index=True)
    A.to_csv(f"{RES}/dynamics_fcd_recomputed.csv", index=False)
    E = A[A.param_set >= 0]
    rows = []
    for (nm, F), g in E.groupby(["topology", "fold"]):
        t = reg[nm]
        rows.append(dict(topology=nm, census_class=t.core_class,
                         n_sign_resolved_cores_in_BRCA_network=t.n_cores_in_network,
                         fold_change=F, n_param_sets=len(g),
                         pct_FCD_peak_criterion=100 * float(g.is_FCD_peak.mean()),
                         pct_FCD_end_criterion_old=100 * float(g.is_FCD_end.mean()),
                         median_FCD_error=float(np.nanmedian(g.FCD_error)),
                         pct_trajectory_invariant_err_lt_0p10=100 * float((g.FCD_error < 0.10).mean()),
                         median_peak_excursion=float(np.nanmedian(g.peak_excursion)),
                         median_end_excursion=float(np.nanmedian(g.end_excursion))))
    S = pd.DataFrame(rows)
    S.to_csv(f"{RES}/dynamics_fcd_recomputed_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print()
    print(S[S.fold_change == 3.0].round(3).to_string(index=False))
    print(f"\nDONE in {(time.time()-t0)/60:.1f} min")
