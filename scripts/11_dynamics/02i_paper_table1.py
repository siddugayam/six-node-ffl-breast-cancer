#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02i_paper_table1.py
A compact, print-ready version of the 3-node comparison table: the columns a
journal table would actually carry, formatted as text, with the full 46-column
version left in dynamics_3node_comparison_table.csv for the supplement.

Out: results/v3/dynamics_TABLE1_paper.csv
"""
import os, numpy as np, pandas as pd

RES = "/path/to/revision/results/v3"
S = pd.read_csv(f"{RES}/dynamics_3node_comparison_table.csv")
O = None
if os.path.exists(f"{RES}/dynamics_3node_oscillation_by_titration.csv"):
    O = pd.read_csv(f"{RES}/dynamics_3node_oscillation_by_titration.csv")
L = None
if os.path.exists(f"{RES}/dynamics_pulse_longwindow_summary.csv"):
    L = pd.read_csv(f"{RES}/dynamics_pulse_longwindow_summary.csv").set_index("topology")

rows = []
for _, r in S.iterrows():
    t = r.topology
    cap = ""
    if O is not None:
        o0 = O[(O.topology == t) & (O.titration.str.startswith("theta = 0"))]
        o1 = O[(O.topology == t) & (O.titration.str.startswith("theta sampled"))]
        if len(o0) and len(o1):
            cap = (f"complex eigenvalues: {float(o0.pct_any_complex_eigenvalue.iloc[0]):.1f}% "
                   f"(catalytic RISC) / {float(o1.pct_any_complex_eigenvalue.iloc[0]):.1f}% "
                   f"(with titration)")
    fwhm = ""
    if L is not None and t in L.index and np.isfinite(L.loc[t, "fwhm_median_tau_uncensored"]):
        fwhm = (f"{L.loc[t,'fwhm_median_tau_uncensored']:.1f} "
                f"(n={int(L.loc[t,'n_uncensored'])} uncensored)")
    rows.append(dict(
        topology=t,
        circuit=r.label,
        census_class=r.census_class,
        sign_resolved_cores_in_BRCA_network=int(r.n_sign_resolved_cores_in_BRCA_network),
        ON_response_vs_cascade=f"{r.rel_T_ON_median:.2f} [{r.rel_T_ON_q25:.2f}-{r.rel_T_ON_q75:.2f}]",
        pct_accelerated_ON=round(float(r.pct_accelerated_ON), 1),
        OFF_response_vs_cascade=f"{r.rel_T_OFF_median:.2f} [{r.rel_T_OFF_q25:.2f}-{r.rel_T_OFF_q75:.2f}]",
        sign_sensitive_delay_index_log2=round(float(r.SSD_index_median), 3),
        pct_parameter_space_pulsing=round(float(r.pct_pulse), 1),
        overshoot_when_pulsing=("-" if not np.isfinite(r.overshoot_ratio_median_when_pulse)
                                else f"{r.overshoot_ratio_median_when_pulse:.2f}"),
        adaptation_when_pulsing=("-" if not np.isfinite(r.adaptation_median_when_pulse)
                                 else f"{r.adaptation_median_when_pulse:.2f}"),
        pulse_FWHM_tau_uncensored=fwhm if fwhm else "-",
        pct_fold_change_detecting=round(float(r.pct_FCD), 2),
        noise_filter_index=round(float(r.filter_index_median), 3),
        steady_state_steepness_n_eff=round(float(r.n_eff_protein_median), 2),
        dynamic_range_floored=round(float(r.dyn_range_protein_floored_median), 1),
        pct_nonmonotone_dose_response=round(float(r.pct_nonmonotone_dose), 1),
        oscillatory_capacity=cap if cap else "-",
    ))
T = pd.DataFrame(rows)
T.to_csv(f"{RES}/dynamics_TABLE1_paper.csv", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 30)
print(T[["topology", "census_class", "sign_resolved_cores_in_BRCA_network",
         "ON_response_vs_cascade", "OFF_response_vs_cascade",
         "sign_sensitive_delay_index_log2", "pct_parameter_space_pulsing",
         "overshoot_when_pulsing", "adaptation_when_pulsing",
         "pct_fold_change_detecting", "noise_filter_index",
         "steady_state_steepness_n_eff", "dynamic_range_floored"]].to_string(index=False))
print(f"\nwrote {len(T)} rows, {T.shape[1]} columns -> dynamics_TABLE1_paper.csv")
