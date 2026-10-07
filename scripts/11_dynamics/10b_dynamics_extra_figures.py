#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10b_dynamics_extra_figures.py
Figures for the later dynamics analyses (07c_consolidate_additions.py).

Fig_dyn9   the one qualitative gain of function from higher order: a certified
           limit cycle in the 6-node composite module whose own 3-node core is
           dead flat at the same parameters
Fig_dyn10  intrinsic-noise buffering measured with chemical Langevin equations
           (the Osella et al. 2011 prediction, tested)
Fig_dyn11  paired effect sizes versus the matched cascade control, and the
           peak-amplitude fold-change-detection test (02g)
Every panel is skipped, with a message, if its input file is missing.
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_style as ST
ST.apply()
import matplotlib.pyplot as plt

RES = "/path/to/revision/results/v3"
FIG = "/path/to/revision/figures/v3"
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/{name}.{ext}")
    plt.close(fig)
    print("  wrote", name)


def have(*fs):
    miss = [f for f in fs if not os.path.exists(f"{RES}/{f}")]
    if miss:
        print("  SKIP - missing:", ", ".join(miss))
        return False
    return True


# ======================================================== FIGURE 9 ============
if have("dynamics_n6_limit_cycle_timecourses.csv", "dynamics_n6_limit_cycle_certified.csv",
        "dynamics_higher_order_fractions.csv"):
    TC = pd.read_csv(f"{RES}/dynamics_n6_limit_cycle_timecourses.csv")
    CE = pd.read_csv(f"{RES}/dynamics_n6_limit_cycle_certified.csv")
    F = pd.read_csv(f"{RES}/dynamics_higher_order_fractions.csv")
    fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.9))

    ax = axes[0]
    for mod, col, lab in (("n3", ST.MUTED, "3-node core"), ("n6", ST.VIOLET, "6-node module")):
        g = TC[TC.module == mod].sort_values("t_hours")
        if len(g):
            mu = float(g.target_protein.mean())
            ax.plot(g.t_hours - g.t_hours.min(), g.target_protein / max(mu, 1e-12),
                    color=col, lw=1.7)
            ax.text(0.98, 0.94 if mod == "n6" else 0.10,
                    f"{lab} (mean {mu:.3g})", transform=ax.transAxes,
                    ha="right", va="center", color=col, fontsize=6.8, fontweight="bold")
    ax.axhline(1.0, color=ST.AXIS, lw=0.7, ls=(0, (4, 3)))
    ax.set_xlabel("time (h), last 60 $\\tau$ of a 1,000 $\\tau$ run")
    ax.set_ylabel("target protein / its own mean")
    ax.set_title("same parameter vector, two module sizes", loc="left", color=ST.INK, fontsize=8)
    ST.panel_label(ax, "A", dx=-0.30)

    ax = axes[1]
    mods = ["n3", "n4", "n5", "n6"]
    f2 = f"{RES}/dynamics_higher_order_fractions_compI1.csv"
    if os.path.exists(f2):
        F = pd.concat([F, pd.read_csv(f2)], ignore_index=True)
    fams = list(F.family.unique())
    cmap = {"COMP_C2_toggle": ST.VIOLET, "I1_miRNA_FFL": ST.AQUA,
            "COMP_I1_negfeedback": ST.ORANGE}
    w = 0.8 / len(fams)
    for k, fam in enumerate(fams):
        sub = F[F.family == fam].set_index("module")
        v = [sub.loc[m, "pct_is_sustained_osc"] if m in sub.index else 0 for m in mods]
        ax.bar(np.arange(len(mods)) + (k - (len(fams) - 1) / 2) * w, v, width=w * 0.9,
               color=cmap.get(fam, ST.MUTED), label=fam.replace("_", " "))
    ax.set_xticks(np.arange(len(mods))); ax.set_xticklabels([m[1] + "-node" for m in mods])
    ax.set_ylabel("% of 16,384 parameter sets\nwith sustained oscillation")
    ax.legend(fontsize=6.6, loc="upper left")
    ax.set_title("sustained oscillation appears only at n = 6", loc="left",
                 color=ST.INK, fontsize=8)
    ST.panel_label(ax, "B", dx=-0.30)

    ax = axes[2]
    steps = ["flagged\nby screen", "amplitude\nheld to\n1,000 $\\tau$",
             "same from\n2 initial\nconditions", "CERTIFIED\nlimit cycle",
             "3-node core\nalso\noscillates"]
    vals = [len(CE), int(CE.amplitude_sustained_1000tau.sum()),
            int(CE.independent_of_initial_condition.sum()),
            int(CE.certified_limit_cycle.sum()),
            int((CE.certified_limit_cycle & CE.core_3node_oscillates).sum())]
    cols = [ST.SEQ[4], ST.SEQ[6], ST.SEQ[8], ST.VIOLET, ST.MUTED]
    ax.bar(np.arange(len(vals)), vals, color=cols, width=0.68)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.5, str(v), ha="center", va="bottom", fontsize=7.5,
                fontweight="bold", color=ST.INK2)
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(steps, fontsize=5.8, rotation=24, ha="right", rotation_mode="anchor")
    ax.set_ylabel("parameter sets")
    ax.set_ylim(0, max(vals) * 1.25)
    ax.set_title("certification of the 6-node oscillators", loc="left", color=ST.INK, fontsize=8)
    ST.panel_label(ax, "C", dx=-0.30)
    fig.tight_layout(w_pad=2.0)
    save(fig, "Fig_dyn9_order_adds_oscillation")

# ======================================================== FIGURE 10 ===========
if have("dynamics_stochastic_noise.csv", "dynamics_stochastic_noise_summary.csv"):
    A = pd.read_csv(f"{RES}/dynamics_stochastic_noise.csv")
    S = pd.read_csv(f"{RES}/dynamics_stochastic_noise_summary.csv")
    om = 500.0
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
    ax = axes[0]
    style = {"CASCADE": (ST.MUTED, "cascade"), "I1_TXN_AND": (ST.BLUE, "I1 transcriptional"),
             "I1_MIR_AND": (ST.AQUA, "I1 miRNA"), "COMP_I1_AND": (ST.VIOLET, "composite I1")}
    sub = A[A.omega == om]
    for nm, (c, lab) in style.items():
        g = sub[(sub.topology == nm) & (sub.mean_protein > 1e-3) & (sub.cv_protein > 0)]
        if not len(g):
            continue
        ax.scatter(g.mean_protein, g.cv_protein, s=2.0, alpha=0.20, color=c, linewidths=0)
        x = np.log(g.mean_protein.to_numpy()); y = np.log(g.cv_protein.to_numpy())
        b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(np.exp(xs), np.exp(np.polyval(b, xs)), color=c, lw=1.8, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("mean target protein (scaled)")
    ax.set_ylabel("CV of target protein")
    ax.set_title(f"intrinsic noise, chemical Langevin, $\\Omega$ = {om:.0f}",
                 loc="left", color=ST.INK, fontsize=8)
    ax.legend(fontsize=6.6, loc="lower left")
    ST.panel_label(ax, "A", dx=-0.22)

    ax = axes[1]
    order = [t for t in ["I1_TXN_AND", "I1_MIR_AND", "C2_MIR_AND", "COMP_I1_AND", "COMP_C2_AND"]
             if t in set(S.topology)]
    MM = None
    mmf = f"{RES}/dynamics_stochastic_noise_matched_mean.csv"
    if os.path.exists(mmf):
        MM = pd.read_csv(mmf)
    w = 0.38
    raw = [float(S[(S.topology == t) & (S.omega == om)]["median_CV_ratio_to_cascade"].iloc[0])
           for t in order]
    ax.bar(np.arange(len(order)) - w / 2, raw, width=w * 0.9, color=ST.SEQ[4],
           label="as simulated")
    if MM is not None:
        mm = [float(MM[(MM.topology == t) & (MM.omega == om)]
                    ["median_CV_ratio_at_matched_mean"].iloc[0]) for t in order]
        ax.bar(np.arange(len(order)) + w / 2, mm, width=w * 0.9, color=ST.VIOLET,
               label="at matched mean output")
    ax.axhline(1.0, color=ST.INK, lw=0.9)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([t.replace("_AND", "").replace("_", " ") for t in order],
                       rotation=20, ha="right", fontsize=6.6)
    ax.set_ylabel("median CV / cascade CV\n(<1 = quieter than a cascade)")
    ax.legend(fontsize=6.6)
    ax.set_title(f"paired at identical parameters, $\\Omega$ = {om:.0f}",
                 loc="left", color=ST.INK, fontsize=8)
    ST.panel_label(ax, "B", dx=-0.24)
    fig.tight_layout(w_pad=2.0)
    save(fig, "Fig_dyn10_stochastic_noise")

# ======================================================== FIGURE 11 ===========
if have("dynamics_paired_tests_vs_cascade.csv"):
    P = pd.read_csv(f"{RES}/dynamics_paired_tests_vs_cascade.csv")
    fcd_ok = os.path.exists(f"{RES}/dynamics_fcd_recomputed_summary.csv")
    ncol = 2 if fcd_ok else 1
    fig, axes = plt.subplots(1, ncol, figsize=(7.6 if fcd_ok else 4.2, 4.6))
    axes = np.atleast_1d(axes)

    ax = axes[0]
    sub = P[P.metric.isin(["T_half_ON", "T_half_OFF"])]
    tops = [t for t in sub.topology.unique()]
    y = np.arange(len(tops))
    for met, col, off in (("T_half_ON", ST.ORANGE, -0.18), ("T_half_OFF", ST.BLUE, 0.18)):
        g = sub[sub.metric == met].set_index("topology").reindex(tops)
        ax.errorbar(g.median_effect, y + off,
                    xerr=[g.median_effect - g.ci_lo, g.ci_hi - g.median_effect],
                    fmt="o", ms=3.2, lw=1.0, color=col, capsize=1.6,
                    label="ON step" if met == "T_half_ON" else "OFF step")
    ax.axvline(1.0, color=ST.INK, lw=0.9)
    ax.set_yticks(y); ax.set_yticklabels(tops, fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlabel("response half-time / cascade half-time\n(<1 accelerated, >1 delayed)")
    ax.legend(fontsize=6.8, loc="lower right")
    ax.set_title("paired against the cascade at identical parameters",
                 loc="left", color=ST.INK, fontsize=8)
    ST.panel_label(ax, "A", dx=-0.55)

    if fcd_ok:
        S = pd.read_csv(f"{RES}/dynamics_fcd_recomputed_summary.csv")
        S3 = S[S.fold_change == 3.0].set_index("topology")
        tops2 = [t for t in tops if t in S3.index]
        y2 = np.arange(len(tops2))
        ax = axes[1]
        ax.barh(y2 - 0.19, S3.loc[tops2, "pct_FCD_end_criterion_old"], height=0.34,
                color=ST.MUTED, label="end-point criterion (old)")
        ax.barh(y2 + 0.19, S3.loc[tops2, "pct_FCD_peak_criterion"], height=0.34,
                color=ST.GREEN, label="peak criterion (corrected)")
        ax.set_yticks(y2); ax.set_yticklabels(tops2, fontsize=6.4)
        ax.invert_yaxis()
        ax.set_xlabel("% of parameter sets showing fold-change detection\n(3-fold input step)")
        ax.legend(fontsize=6.8, loc="lower right")
        ax.set_title("fold-change detection, corrected responsiveness test",
                     loc="left", color=ST.INK, fontsize=8)
        ST.panel_label(ax, "B", dx=-0.50)
    fig.tight_layout(w_pad=2.6)
    save(fig, "Fig_dyn11_paired_and_fcd")

print("done")
