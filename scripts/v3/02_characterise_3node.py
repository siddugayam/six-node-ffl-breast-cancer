#!/usr/bin/env python3
"""
02_characterise_3node.py  -- Part A + B.
Simulate every 3-node topology (Alon transcriptional references, miRNA-mediated
cores, composite cores, and the regulated-cascade control) at (i) a stated
reference parameter set and (ii) a matched scrambled-Sobol ensemble, and compute
the full dynamical battery.  Compute only -- figures are made by 05_figures.py.

Outputs (results/v3/):
  dynamics_parameter_definitions.csv
  dynamics_parameter_ranges.csv
  dynamics_topologies.csv
  dynamics_reference_timecourses.csv
  dynamics_reference_doseresponse.csv
  dynamics_3node_metrics_ensemble.csv   (per parameter set, per topology)
  dynamics_3node_comparison_table.csv   (summary; the paper table)
"""
import sys, os, time, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_metrics as M, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
N_ENS = 1024
os.makedirs(RES, exist_ok=True)
t0 = time.time()

# ---------------------------------------------------------------- metadata --
pd.DataFrame(D.PARAM_UNITS, columns=["symbol", "quantity", "units", "meaning"]) \
  .to_csv(f"{RES}/dynamics_parameter_definitions.csv", index=False)
SMP.ranges_table().to_csv(f"{RES}/dynamics_parameter_ranges.csv", index=False)

reg = D.build_registry()
pd.DataFrame([dict(topology=t.name, label=t.label, n_nodes=t.n_nodes, gate=t.gate,
                   miRNA_arm_post_transcriptional=t.mir_arm,
                   sign_TF_to_intermediate=t.s_TM, sign_intermediate_to_target=t.s_MY,
                   sign_TF_to_target=t.s_TY, reciprocal_miRNA_represses_TF=t.mir_to_TF,
                   signal_enters_at=t.input_node,
                   n_sign_resolved_cores_in_BRCA_network=t.n_cores_in_network,
                   census_class=t.core_class, note=t.note)
              for t in reg.values()]).to_csv(f"{RES}/dynamics_topologies.csv", index=False)

# --------------------------------------------------------------- ensemble ---
# Column 0 of every parameter array is the REFERENCE parameter set; columns
# 1..N_ENS are the scrambled-Sobol ensemble.  One pass therefore yields both the
# reference time courses (for the figures) and the ensemble metrics.
Pref = D.default_params(1)
Pens = SMP.sample(N_ENS, seed=20260908)
Pall = {k: np.concatenate([np.asarray(Pref[k], float), np.asarray(Pens[k], float)])
        for k in Pref}
NT = N_ENS + 1
np.savez_compressed(f"{RES}/dynamics_3node_sobol_params.npz", **Pall)
pd.DataFrame({"parameter": list(Pref), "reference_value": [float(v[0]) for v in Pref.values()]}) \
  .to_csv(f"{RES}/dynamics_reference_parameters.csv", index=False)

S_grid = np.logspace(np.log10(0.02), np.log10(20.0), 41)


def run_one(nm):
    t = reg[nm]
    met, aux = M.characterise(t, Pall, NT, S_grid=S_grid)
    tc, dr = [], []
    for i, tt in enumerate(aux["ts"]):
        tc.append(dict(topology=nm, protocol="ON_step", t=tt, protein=aux["p_on"][i, 0]))
    for i, tt in enumerate(aux["ts2"]):
        tc.append(dict(topology=nm, protocol="OFF_step", t=tt, protein=aux["p_off"][i, 0]))
    for i, tt in enumerate(aux["ts3"]):
        tc.append(dict(topology=nm, protocol="brief_pulse_0.5tau", t=tt, protein=aux["p_short"][i, 0]))
    for i, S in enumerate(S_grid):
        dr.append(dict(topology=nm, S=S, protein_ss=aux["Ps"][i, 0], mRNA_ss=aux["Ys"][i, 0]))
    df = pd.DataFrame({k: np.asarray(v, float) for k, v in met.items()})
    df.insert(0, "param_set", np.arange(NT) - 1)          # -1 = reference set
    df.insert(0, "topology", nm)
    print(f"  {nm:16s} done  ({time.time()-t0:6.1f}s)", flush=True)
    return tc, dr, df


if __name__ == "__main__":
    with Pool(processes=4) as pool:
        parts = pool.map(run_one, list(reg))
    rows_tc = [r for p_ in parts for r in p_[0]]
    rows_dr = [r for p_ in parts for r in p_[1]]
    all_rows = [p_[2] for p_ in parts]

    pd.DataFrame(rows_tc).to_csv(f"{RES}/dynamics_reference_timecourses.csv", index=False)
    pd.DataFrame(rows_dr).to_csv(f"{RES}/dynamics_reference_doseresponse.csv", index=False)
    E = pd.concat(all_rows, ignore_index=True)

    # ---- relative-to-cascade response times, matched by parameter set -----------
    cas = E[E.topology == "CASCADE"].set_index("param_set")
    E = E.merge(cas[["T_half_ON", "T_half_OFF"]].rename(
            columns={"T_half_ON": "T_half_ON_cascade", "T_half_OFF": "T_half_OFF_cascade"}),
            left_on="param_set", right_index=True, how="left")
    E["rel_T_ON"] = E.T_half_ON / E.T_half_ON_cascade
    E["rel_T_OFF"] = E.T_half_OFF / E.T_half_OFF_cascade
    E["sign_sensitive_delay_index"] = np.log2(E.rel_T_ON / E.rel_T_OFF)
    E.to_csv(f"{RES}/dynamics_3node_metrics_ensemble.csv", index=False)

    # ---- summary table ----------------------------------------------------------
    def q(s, p):
        return np.nanpercentile(s.astype(float), p) if np.isfinite(s.astype(float)).any() else np.nan

    summ = []
    Eens = E[E.param_set >= 0]
    for nm, g in Eens.groupby("topology", sort=False):
        t = reg[nm]
        summ.append(dict(
            topology=nm, label=t.label, gate=t.gate,
            n_cores_in_BRCA_network=t.n_cores_in_network,
            n_param_sets=len(g),
            rel_T_ON_median=q(g.rel_T_ON, 50), rel_T_ON_q25=q(g.rel_T_ON, 25), rel_T_ON_q75=q(g.rel_T_ON, 75),
            rel_T_OFF_median=q(g.rel_T_OFF, 50), rel_T_OFF_q25=q(g.rel_T_OFF, 25), rel_T_OFF_q75=q(g.rel_T_OFF, 75),
            pct_accelerated_ON=100 * np.nanmean(g.rel_T_ON < 1),
            pct_delayed_ON=100 * np.nanmean(g.rel_T_ON > 1),
            pct_accelerated_OFF=100 * np.nanmean(g.rel_T_OFF < 1),
            pct_delayed_OFF=100 * np.nanmean(g.rel_T_OFF > 1),
            SSD_index_median=q(g.sign_sensitive_delay_index, 50),
            pct_pulse=100 * np.nanmean(g.is_pulse),
            overshoot_ratio_median=q(g.overshoot_ratio, 50),
            pulse_width_above_half_peak_median=q(g.pulse_width_above_half_peak, 50),
            adaptation_median=q(g.adaptation, 50),
            pct_adaptation_gt_0p5=100 * np.nanmean(g.adaptation > 0.5),
            FCD_error_median=q(g.FCD_error, 50), pct_FCD=100 * np.nanmean(g.is_FCD),
            FCD_response_amplitude_median=q(g.FCD_response_amplitude, 50),
            noise_transmission_median=q(g.noise_transmission, 50),
            filter_index_median=q(g.filter_index, 50),
            pct_filter_index_gt_0p5=100 * np.nanmean(g.filter_index > 0.5),
            n_eff_protein_median=q(g.n_eff_protein, 50),
            pct_ultrasensitive_neff_gt2=100 * np.nanmean(g.n_eff_protein > 2),
            dyn_range_protein_median=q(g.dyn_range_protein, 50),
        dyn_range_protein_floored_median=q(g.dyn_range_protein_floored, 50),
        ss_min_protein_median=q(g.ss_min_protein, 50),
        ss_max_protein_median=q(g.ss_max_protein, 50),
            pct_nonmonotone_dose=100 * np.nanmean(g.nonmonotone_dose),
            n_turning_points_median=q(g.n_turning_points, 50),
        ))
    pd.DataFrame(summ).to_csv(f"{RES}/dynamics_3node_comparison_table.csv", index=False)
    print(f"\nDONE in {time.time()-t0:.1f}s")
    print(pd.DataFrame(summ)[["topology", "rel_T_ON_median", "rel_T_OFF_median",
                              "pct_pulse", "overshoot_ratio_median_when_pulse",
                              "pct_FCD", "filter_index_median",
                              "n_eff_protein_median",
                              "dyn_range_protein_floored_median"]].to_string(index=False))
