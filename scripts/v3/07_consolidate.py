#!/usr/bin/env python3
"""
07_consolidate.py -- pull every headline number of the dynamical analysis into one
table with its provenance (which script and which file produced it).
Out: results/v3/dynamics_key_results.csv
"""
import os, numpy as np, pandas as pd
RES = "/path/to/revision/results/v3"
rows = []
def add(section, statement, value, units, source):
    rows.append(dict(section=section, statement=statement, value=value,
                     units=units, source_file=source))

# ---- empirical grounding -----------------------------------------------------
G = pd.read_csv(f"{RES}/dynamics_empirical_grounding.csv").set_index("ffl_class")
add("grounding", "3-node FFL cores in the BRCA network", int(G.loc['ALL','n_cores']), "cores",
    "dynamics_empirical_grounding.csv")
add("grounding", "cores with all three signs resolved", int(G.loc['ALL','n_resolved']), "cores",
    "dynamics_empirical_grounding.csv")
add("grounding", "coherent fraction of resolved cores", round(float(G.loc['ALL','pct_coherent']),1), "%",
    "dynamics_empirical_grounding.csv")
add("grounding", "cores of Alon type C1 (all-activating) anywhere in the network",
    int(G.loc['ALL','n_C1']), "cores", "dynamics_empirical_grounding.csv")
add("grounding", "I1 cores (commonest resolved type)", int(G.loc['ALL','n_I1']), "cores",
    "dynamics_empirical_grounding.csv")
add("grounding", "miRNA-FFL coherent fraction", round(float(G.loc['miRNA-FFL','pct_coherent']),1), "%",
    "dynamics_empirical_grounding.csv")
add("grounding", "composite-FFL coherent fraction", round(float(G.loc['Composite-FFL','pct_coherent']),1), "%",
    "dynamics_empirical_grounding.csv")
SP = pd.read_csv(f"{RES}/dynamics_composite_feedback_split.csv")
for _, r in SP.iterrows():
    add("grounding", f"composite cores with {r.subclass}", int(r.n_resolved),
        f"resolved cores ({r.pct_of_resolved_composite:.1f}% of resolved composite)",
        "dynamics_composite_feedback_split.csv")

# ---- integrator validation ----------------------------------------------------
if os.path.exists(f"{RES}/dynamics_integrator_validation.csv"):
    V = pd.read_csv(f"{RES}/dynamics_integrator_validation.csv")
    v = V[V.dt == 0.01]
    add("validation", "max abs deviation of the vectorised RK4 (dt=0.01) from scipy LSODA "
        "(rtol 1e-9) across all topologies", float(v.max_abs_err.max()), "concentration units",
        "dynamics_integrator_validation.csv")

# ---- 3-node comparison --------------------------------------------------------
if os.path.exists(f"{RES}/dynamics_3node_comparison_table.csv"):
    S = pd.read_csv(f"{RES}/dynamics_3node_comparison_table.csv").set_index("topology")
    add("partB", "parameter sets per topology (scrambled Sobol)", int(S.n_param_sets.iloc[0]),
        "sets", "dynamics_3node_comparison_table.csv")
    for t in S.index:
        add("partB", f"{t}: median ON-step response time relative to the cascade",
            round(float(S.loc[t,'rel_T_ON_median']),3), "ratio", "dynamics_3node_comparison_table.csv")
        add("partB", f"{t}: median OFF-step response time relative to the cascade",
            round(float(S.loc[t,'rel_T_OFF_median']),3), "ratio", "dynamics_3node_comparison_table.csv")
        add("partB", f"{t}: % of parameter space showing a pulse",
            round(float(S.loc[t,'pct_pulse']),1), "%", "dynamics_3node_comparison_table.csv")
        add("partB", f"{t}: median noise-filter index (1 = brief input fully rejected)",
            round(float(S.loc[t,'filter_index_median']),3), "index", "dynamics_3node_comparison_table.csv")
        add("partB", f"{t}: % of parameter space with fold-change detection",
            round(float(S.loc[t,'pct_FCD']),1), "%", "dynamics_3node_comparison_table.csv")
        add("partB", f"{t}: median effective Hill coefficient of the steady-state dose-response",
            round(float(S.loc[t,'n_eff_protein_median']),2), "n_eff", "dynamics_3node_comparison_table.csv")

# ---- higher order --------------------------------------------------------------
if os.path.exists(f"{RES}/dynamics_higher_order_fractions.csv"):
    F = pd.read_csv(f"{RES}/dynamics_higher_order_fractions.csv")
    add("partC", "parameter sets per module (scrambled Sobol)", int(F.n_param_sets.iloc[0]),
        "sets", "dynamics_higher_order_fractions.csv")
    for _, r in F.iterrows():
        for b in [c for c in F.columns if c.startswith("pct_is_")]:
            add("partC", f"{r.family} / {r.module}: {b.replace('pct_is_','')}",
                round(float(r[b]),2), "% of parameter sets", "dynamics_higher_order_fractions.csv")
    GLf = pd.read_csv(f"{RES}/dynamics_higher_order_gain_loss.csv")
    gb = GLf[(GLf.metric_type == "binary")]
    for _, r in gb.iterrows():
        add("partC_gain_loss", f"{r.family} / {r.module} vs its own 3-node core: {r.behaviour}",
            f"gain {r.pct_gain:.2f}% , loss {r.pct_loss:.2f}% (McNemar p={r.mcnemar_p:.3g})",
            "% of parameter sets", "dynamics_higher_order_gain_loss.csv")

# ---- miR-29 -------------------------------------------------------------------
if os.path.exists(f"{RES}/dynamics_mir29_fit.csv"):
    FI = pd.read_csv(f"{RES}/dynamics_mir29_fit.csv").set_index("parameter")
    for k in FI.index:
        add("partD", k, round(float(FI.loc[k,'estimate']),4)
            if np.isfinite(FI.loc[k,'estimate']) else np.nan,
            str(FI.loc[k,'units']), "dynamics_mir29_fit.csv")
    P = pd.read_csv(f"{RES}/dynamics_mir29_predictions.csv")
    key = P[(P.cell_state == "activated myofibroblast") & (P.protein_half_life_h == 46.0)]
    for _, r in key.iterrows():
        add("partD_prediction",
            f"miR-29 raised {r.miR29_fold_over_cohort_median:.2f}x over the cohort median "
            f"(activated myofibroblast, COL1A1 mRNA t1/2 24 h, protein t1/2 46 h)",
            f"COL1A1 mRNA falls {r.steady_state_COL1A1_mRNA_pct_drop:.1f}%; "
            f"mRNA t1/2 {r.t_half_COL1A1_mRNA_hours:.1f} h; protein t1/2 "
            f"{r.t_half_COL1A1_protein_hours:.1f} h; protein 90% at {r.t90_COL1A1_protein_hours:.0f} h",
            "prediction", "dynamics_mir29_predictions.csv")

# ---- model well-posedness -------------------------------------------------------
if os.path.exists(f"{RES}/dynamics_model_wellposedness_check.csv"):
    W = pd.read_csv(f"{RES}/dynamics_model_wellposedness_check.csv")
    for var, g in W.groupby("mutual_stabilisation"):
        gg = g[g.module.isin(["n4", "n5", "n6"])]
        add("model_check",
            f"gene-gene edge modelled as {var} mutual stabilisation: parameter sets whose "
            f"4/5/6-node module had NOT reached a steady state after 250 tau",
            round(float(gg.pct_rel_residual_gt_1e4.mean()), 1), "% (mean over modules)",
            "dynamics_model_wellposedness_check.csv")
        add("model_check",
            f"gene-gene edge modelled as {var} mutual stabilisation: largest target-protein "
            f"level reached in any 4/5/6-node module",
            round(float(gg.p1_max.max()), 1), "scaled concentration",
            "dynamics_model_wellposedness_check.csv")

# ---- oscillation capacity -------------------------------------------------------
if os.path.exists(f"{RES}/dynamics_3node_oscillation_by_titration.csv"):
    O = pd.read_csv(f"{RES}/dynamics_3node_oscillation_by_titration.csv")
    for _, r in O.iterrows():
        add("partB_oscillation",
            f"{r.topology} ({r.titration}): parameter sets with a complex eigenvalue at the "
            f"operating point (capacity to ring or oscillate)",
            round(float(r.pct_any_complex_eigenvalue), 2), "% of parameter sets",
            "dynamics_3node_oscillation_by_titration.csv")
        add("partB_oscillation",
            f"{r.topology} ({r.titration}): parameter sets at an unstable focus "
            f"(sustained oscillation)",
            round(float(r.pct_hopf_unstable_focus), 2), "% of parameter sets",
            "dynamics_3node_oscillation_by_titration.csv")

# ---- paired miRNA vs transcriptional --------------------------------------------
if os.path.exists(f"{RES}/dynamics_mir_vs_txn_paired.csv"):
    PA = pd.read_csv(f"{RES}/dynamics_mir_vs_txn_paired.csv")
    for _, r in PA.iterrows():
        add("partB_paired", f"{r.comparison}: {r.metric}",
            f"median {r.median_A:.3f} vs {r.median_B:.3f}; paired difference "
            f"{r.median_paired_difference:+.3f} ({r.pct_A_greater:.1f}% of sets higher; "
            f"Wilcoxon p={r.wilcoxon_p:.3g})", "paired over parameter sets",
            "dynamics_mir_vs_txn_paired.csv")

# ---- analytic bounds -------------------------------------------------------------
if os.path.exists(f"{RES}/dynamics_analytic_bounds.csv"):
    AB = pd.read_csv(f"{RES}/dynamics_analytic_bounds.csv")
    for _, r in AB.iterrows():
        add("analytic_bounds", r.statement,
            (f"predicted {r.median_predicted:.4g}; observed {r.median_observed:.4g}; "
             f"median relative error {r.median_relative_error:.2e}")
            if np.isfinite(r.median_observed) else f"predicted {r.median_predicted:.4g}",
            "closed form vs simulation", "dynamics_analytic_bounds.csv")

# ---- references ----------------------------------------------------------------
R = pd.read_csv(f"{RES}/dynamics_references_verified.csv")
add("references", "citations resolved on CrossRef", f"{(R.status=='FOUND').sum()}/{len(R)}",
    "references", "dynamics_references_verified.csv")

out = pd.DataFrame(rows)
out.to_csv(f"{RES}/dynamics_key_results.csv", index=False)
print(f"wrote {len(out)} rows -> dynamics_key_results.csv")
print(out[out.section.isin(["grounding","validation","references"])].to_string(index=False))
