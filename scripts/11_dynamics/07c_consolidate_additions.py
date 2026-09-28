#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07c_consolidate_additions.py
Appends the results produced after the first pass to dynamics_key_results.csv, each
with the file it came from, so that every number quoted in the manuscript can be
traced to one line of one CSV.  Run AFTER 07_consolidate.py.
"""
import os, numpy as np, pandas as pd

RES = "/path/to/revision/results/v3"
K = f"{RES}/dynamics_key_results.csv"
rows = []


def add(section, statement, value, units, source):
    rows.append(dict(section=section, statement=statement, value=value,
                     units=units, source_file=source))


def ex(f):
    return os.path.exists(f"{RES}/{f}")


# ---- numerical validation -----------------------------------------------------
if ex("dynamics_solver_crosscheck_summary.csv"):
    S = pd.read_csv(f"{RES}/dynamics_solver_crosscheck_summary.csv").set_index("metric")
    f = "dynamics_solver_crosscheck_summary.csv"
    for m in S.index:
        add("validation", f"{m}: median |fixed-step RK4 - adaptive LSODA| over 320 "
            f"topology x parameter-set combinations", float(S.loc[m, "median_abs_diff"]),
            str(S.loc[m, "unit"]), f)
        add("validation", f"{m}: worst-case |RK4 - LSODA|", float(S.loc[m, "max_abs_diff"]),
            str(S.loc[m, "unit"]), f)
        add("validation", f"{m}: Pearson r between the two solvers",
            round(float(S.loc[m, "pearson_r"]), 6), "r", f)

if ex("dynamics_3node_reproducibility_vs_N1024.csv"):
    R = pd.read_csv(f"{RES}/dynamics_3node_reproducibility_vs_N1024.csv")
    add("validation", "largest disagreement between the 16,384-set run and the earlier "
        "1,024-set run on their 1,024 shared parameter sets (Sobol prefix property)",
        float(R.max_abs_diff.max()), "metric units",
        "dynamics_3node_reproducibility_vs_N1024.csv")
    add("validation", "metrics with any finite/NaN mismatch between the two runs",
        int((R.n_finiteness_mismatch > 0).sum()), "metrics",
        "dynamics_3node_reproducibility_vs_N1024.csv")

# ---- textbook validation ------------------------------------------------------
if ex("dynamics_textbook_validation.csv"):
    V = pd.read_csv(f"{RES}/dynamics_textbook_validation.csv")
    add("validation", "published FFL results the model was required to reproduce",
        len(V), "criteria", "dynamics_textbook_validation.csv")
    add("validation", "of those, PASS", int((V.verdict == "PASS").sum()), "criteria",
        "dynamics_textbook_validation.csv")
    for _, r in V.iterrows():
        add("validation", f"[{r.verdict}] {r.published_result[:110]}", r.observed, "",
            "dynamics_textbook_validation.csv")

# ---- paired tests -------------------------------------------------------------
if ex("dynamics_paired_tests_vs_cascade.csv"):
    P = pd.read_csv(f"{RES}/dynamics_paired_tests_vs_cascade.csv")
    f = "dynamics_paired_tests_vs_cascade.csv"
    for met in ("T_half_ON", "T_half_OFF"):
        g = P[P.metric == met]
        for _, r in g.iterrows():
            add("partB", f"{r.topology}: {met} relative to the matched cascade "
                f"(median, 95% bootstrap CI; paired Wilcoxon)",
                f"{r.median_effect:.3f} [{r.ci_lo:.3f}, {r.ci_hi:.3f}], "
                f"p_Holm={r.p_holm_within_metric:.3g}", "ratio", f)

# ---- pulse, long window -------------------------------------------------------
if ex("dynamics_pulse_longwindow_summary.csv"):
    L = pd.read_csv(f"{RES}/dynamics_pulse_longwindow_summary.csv")
    f = "dynamics_pulse_longwindow_summary.csv"
    for _, r in L.iterrows():
        if r.n_flagged_pulse_60tau > 0:
            add("partB", f"{r.topology}: median pulse FWHM among UNCENSORED pulses over a "
                f"400-tau window (n={int(r.n_uncensored)})",
                round(float(r.fwhm_median_tau_uncensored), 2), "tau", f)
            add("partB", f"{r.topology}: % of flagged pulses still censored at 400 tau",
                round(float(r.pct_of_flagged_censored_at_400tau), 1), "%", f)
            add("partB", f"{r.topology}: % of non-flagged sets that do pulse over 400 tau "
                f"(false-negative rate of the 60-tau screen)",
                round(float(r.pct_nonflagged_that_pulse_at_400tau), 2), "%", f)

# ---- corrected FCD ------------------------------------------------------------
if ex("dynamics_fcd_recomputed_summary.csv"):
    S = pd.read_csv(f"{RES}/dynamics_fcd_recomputed_summary.csv")
    f = "dynamics_fcd_recomputed_summary.csv"
    for _, r in S[S.fold_change == 3.0].iterrows():
        add("partB", f"{r.topology}: % of parameter space with fold-change detection at a "
            f"3-fold step (corrected peak-excursion criterion)",
            round(float(r.pct_FCD_peak_criterion), 2), "%", f)
    add("partB", "highest % of parameter space with fold-change detection in any topology "
        "(corrected criterion, any fold)", round(float(S.pct_FCD_peak_criterion.max()), 2),
        "%", f)

# ---- higher order -------------------------------------------------------------
if ex("dynamics_n6_limit_cycle_certified.csv"):
    C = pd.read_csv(f"{RES}/dynamics_n6_limit_cycle_certified.csv")
    f = "dynamics_n6_limit_cycle_certified.csv"
    add("partC", "6-node composite sets flagged as sustained oscillators", len(C), "sets", f)
    add("partC", "of those, CERTIFIED limit cycles (amplitude held for 1,000 tau and "
        "identical from two initial conditions)", int(C.certified_limit_cycle.sum()), "sets", f)
    add("partC", "certified limit cycles as a fraction of the 16,384-set design",
        round(100 * float(C.certified_limit_cycle.sum()) / 16384, 3), "%", f)
    add("partC", "certified oscillators whose embedded 3-node core also oscillates",
        int((C.certified_limit_cycle & C.core_3node_oscillates).sum()), "sets", f)
    add("partC", "largest |Im(lambda)| in the 3-node cores of the oscillating parameter sets",
        float(C.n3_lead_eig_im.max()), "1/tau", f)
    cc = C[C.certified_limit_cycle]
    if len(cc):
        add("partC", "median period of the certified limit cycles",
            round(float(np.nanmedian(cc.n6_period_tau)), 2), "tau", f)
        add("partC", "median period of the certified limit cycles, in hours "
            "(tau = 12.98 h, median mammalian mRNA t1/2 9 h)",
            round(float(np.nanmedian(cc.n6_period_hours)), 1), "h", f)

if ex("dynamics_higher_order_empirical_support.csv"):
    E = pd.read_csv(f"{RES}/dynamics_higher_order_empirical_support.csv")
    f = "dynamics_higher_order_empirical_support.csv"
    for _, r in E.iterrows():
        add("partC_grounding", str(r.quantity), r.value, "", f)

if ex("dynamics_higher_order_convergence_audit.csv"):
    A = pd.read_csv(f"{RES}/dynamics_higher_order_convergence_audit.csv")
    f = "dynamics_higher_order_convergence_audit.csv"
    for _, r in A.iterrows():
        add("partC", f"{r.family}/{r.module}: % of sets with at least one input level not at "
            f"a certified steady state after the 20-tau continuation step used in the sweep",
            round(float(r.pct_sets_with_any_uncertified_input_level), 2), "%", f)
        add("partC", f"{r.family}/{r.module}: % of ultrasensitivity calls that change when "
            f"the dose-response is recomputed with certified settling",
            round(float(r.pct_ultrasensitivity_call_changed), 2), "%", f)

if ex("dynamics_higher_order_fractions_compI1.csv"):
    F = pd.read_csv(f"{RES}/dynamics_higher_order_fractions_compI1.csv")
    f = "dynamics_higher_order_fractions_compI1.csv"
    for _, r in F.iterrows():
        for b in ("pct_is_bistable", "pct_is_sustained_osc", "pct_is_ultrasensitive",
                  "pct_is_memory", "pct_is_noise_rejecting"):
            add("partC", f"composite I1 (negative-feedback) family, {r.module}: {b}",
                round(float(r[b]), 3), "%", f)

# ---- stochastic ---------------------------------------------------------------
if ex("dynamics_stochastic_noise_summary.csv"):
    S = pd.read_csv(f"{RES}/dynamics_stochastic_noise_summary.csv")
    f = "dynamics_stochastic_noise_summary.csv"
    for _, r in S.iterrows():
        if r.topology == "CASCADE":
            continue
        if "median_CV_ratio_to_cascade" in S.columns and np.isfinite(r.get("median_CV_ratio_to_cascade", np.nan)):
            add("partB_noise", f"{r.topology} at Omega={r.omega:.0f}: median CV of the target "
                f"protein relative to the matched cascade (<1 = quieter)",
                round(float(r.median_CV_ratio_to_cascade), 4), "ratio", f)
            add("partB_noise", f"{r.topology} at Omega={r.omega:.0f}: % of parameter sets "
                f"quieter than the cascade",
                round(float(r.pct_sets_quieter_than_cascade), 1), "%", f)

A = pd.DataFrame(rows)
if os.path.exists(K):
    old = pd.read_csv(K)
    A = pd.concat([old, A], ignore_index=True)
A.to_csv(K, index=False)
print(f"wrote {len(A)} rows -> dynamics_key_results.csv ({len(rows)} added here)")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 95)
print(pd.DataFrame(rows).head(30).to_string(index=False))
