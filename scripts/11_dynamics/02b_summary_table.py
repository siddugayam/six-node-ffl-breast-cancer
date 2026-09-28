#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02b_summary_table.py
Build the paper's 3-node comparison table from the saved per-parameter-set metrics
(dynamics_3node_metrics_ensemble.csv).  Kept separate from 02 so the table can be
rebuilt without re-running the 45-minute simulation.

Pulse amplitude, width and adaptation are summarised ONLY over the parameter sets that
actually produce a pulse; reporting their median over all sets makes them look like 1.0
and 0.0 merely because most sets are not pulsing.
Out: dynamics_3node_comparison_table.csv, dynamics_parameter_definitions.csv
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D

RES = "/path/to/revision/results/v3"
reg = D.build_registry()
E = pd.read_csv(f"{RES}/dynamics_3node_metrics_ensemble.csv")
Eens = E[E.param_set >= 0]

pd.DataFrame(D.PARAM_UNITS, columns=["symbol", "quantity", "units", "meaning"]) \
  .to_csv(f"{RES}/dynamics_parameter_definitions.csv", index=False)


def q(s, p):
    s = np.asarray(s, float)
    return np.nanpercentile(s, p) if np.isfinite(s).any() else np.nan


summ = []
for nm, g in Eens.groupby("topology", sort=False):
    t = reg[nm]
    gp = g[g.is_pulse > 0]
    summ.append(dict(
        topology=nm, label=t.label, gate=t.gate,
        census_class=t.core_class,
        n_sign_resolved_cores_in_BRCA_network=t.n_cores_in_network,
        n_param_sets=len(g),
        rel_T_ON_median=q(g.rel_T_ON, 50), rel_T_ON_q25=q(g.rel_T_ON, 25), rel_T_ON_q75=q(g.rel_T_ON, 75),
        rel_T_OFF_median=q(g.rel_T_OFF, 50), rel_T_OFF_q25=q(g.rel_T_OFF, 25), rel_T_OFF_q75=q(g.rel_T_OFF, 75),
        pct_accelerated_ON=100 * np.nanmean(g.rel_T_ON < 1),
        pct_delayed_ON=100 * np.nanmean(g.rel_T_ON > 1),
        pct_accelerated_OFF=100 * np.nanmean(g.rel_T_OFF < 1),
        pct_delayed_OFF=100 * np.nanmean(g.rel_T_OFF > 1),
        SSD_index_median=q(g.sign_sensitive_delay_index, 50),
        T_half_ON_median_tau=q(g.T_half_ON, 50), T_half_OFF_median_tau=q(g.T_half_OFF, 50),
        pct_pulse=100 * np.nanmean(g.is_pulse),
        n_pulsing_sets=int((g.is_pulse > 0).sum()),
        overshoot_ratio_median_when_pulse=q(gp.overshoot_ratio, 50) if len(gp) else np.nan,
        overshoot_ratio_q90_when_pulse=q(gp.overshoot_ratio, 90) if len(gp) else np.nan,
        pulse_fwhm_median_tau_when_pulse=q(gp.pulse_fwhm, 50) if len(gp) else np.nan,
        pct_of_pulses_with_censored_fwhm=(100 * np.nanmean(gp.pulse_fwhm_censored)
                                          if len(gp) else np.nan),
        t_peak_median_tau_when_pulse=q(gp.t_peak, 50) if len(gp) else np.nan,
        adaptation_median_when_pulse=q(gp.adaptation, 50) if len(gp) else np.nan,
        pct_adaptation_gt_0p5=100 * np.nanmean(g.adaptation > 0.5),
        overshoot_ratio_median_all_sets=q(g.overshoot_ratio, 50),
        adaptation_median_all_sets=q(g.adaptation, 50),
        FCD_error_median=q(g.FCD_error, 50), pct_FCD=100 * np.nanmean(g.is_FCD),
        FCD_response_amplitude_median=q(g.FCD_response_amplitude, 50),
        noise_transmission_median=q(g.noise_transmission, 50),
        filter_index_median=q(g.filter_index, 50),
        pct_filter_index_gt_0p5=100 * np.nanmean(g.filter_index > 0.5),
        n_eff_protein_median=q(g.n_eff_protein, 50),
        n_eff_mRNA_median=q(g.n_eff_mRNA, 50),
        pct_ultrasensitive_neff_gt2=100 * np.nanmean(g.n_eff_protein > 2),
        dyn_range_protein_floored_median=q(g.dyn_range_protein_floored, 50),
        dyn_range_protein_floored_q90=q(g.dyn_range_protein_floored, 90),
        dyn_range_protein_raw_median=q(g.dyn_range_protein, 50),
        ss_min_protein_median=q(g.ss_min_protein, 50),
        ss_max_protein_median=q(g.ss_max_protein, 50),
        pct_nonmonotone_dose=100 * np.nanmean(g.nonmonotone_dose),
        n_turning_points_median=q(g.n_turning_points, 50),
    ))
S = pd.DataFrame(summ)
S.to_csv(f"{RES}/dynamics_3node_comparison_table.csv", index=False)
pd.set_option("display.width", 250)
print(S[["topology", "census_class", "n_sign_resolved_cores_in_BRCA_network",
         "rel_T_ON_median", "rel_T_OFF_median", "pct_pulse",
         "overshoot_ratio_median_when_pulse", "pulse_fwhm_median_tau_when_pulse",
         "adaptation_median_when_pulse", "pct_FCD", "filter_index_median",
         "n_eff_protein_median", "dyn_range_protein_floored_median",
         "pct_nonmonotone_dose"]].to_string(index=False))
print("\nrows:", len(S), "->", f"{RES}/dynamics_3node_comparison_table.csv")
