#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""05_figures.py -- publication figures for the dynamical analysis.

Fig_dyn1  reference time courses: delay / acceleration / pulse / filtering / dose-response
Fig_dyn2  comparison heat map across all 25 modelled 3-node topologies
Fig_dyn3  response-time distributions over the parameter ensemble
Fig_dyn4  higher-order (n = 3..6): fraction of parameter space showing each behaviour
Fig_dyn5  higher-order: matched gain/loss versus the embedded 3-node core
Fig_dyn6  the fitted miR-29 / COL1A1 circuit and its predictions
Fig_dyn7  what only a COMPOSITE FFL can do: hysteresis and ringing
Fig_dyn8  the titration threshold: a post-transcriptional arm that no promoter can copy
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_style as ST
ST.apply()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RES = "/path/to/revision/results/v3"
FIG = "/path/to/revision/figures/final"
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/{name}.{ext}", dpi=600)
    plt.close(fig)
    print("  wrote", name)


def dlabel(ax, x, y, s, color, ha="left", va="center", dx=0.0, dy=0.0, fs=7.5):
    ax.text(x + dx, y + dy, s, color=color, fontsize=fs, ha=ha, va=va,
            fontweight="bold", clip_on=False)


# ============================================================ FIGURE 1 ========
TC = pd.read_csv(f"{RES}/dynamics_reference_timecourses.csv")
DR = pd.read_csv(f"{RES}/dynamics_reference_doseresponse.csv")


def curve(topo, proto):
    d = TC[(TC.topology == topo) & (TC.protocol == proto)].sort_values("t")
    return d.t.values, d.protein.values


fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
ax = axes[0, 0]
setA = [("CASCADE", ST.MUTED, "cascade (control)"),
        ("C1_TXN_AND", ST.BLUE, "C1 transcriptional (AND)"),
        ("I1_TXN_AND", ST.ORANGE, "I1 transcriptional (AND)"),
        ("I1_MIR_AND", ST.AQUA, "I1 miRNA-mediated (AND)")]
for nm, col, lab in setA:
    t, p = curve(nm, "ON_step")
    ax.plot(t, p / p[-1], color=col, lw=1.8, label=lab)
ax.axhline(1.0, color=ST.AXIS, lw=0.7, ls=(0, (4, 3)))
ax.axhline(0.5, color=ST.GRID, lw=0.7)
ax.set_xlim(0, 20); ax.set_ylim(0, 3.2)
ax.set_xlabel("time (target-mRNA lifetimes, $\\tau$)")
ax.set_ylabel("target protein / its own steady state")
ax.set_title("ON step: delay, acceleration and pulse", loc="left", color=ST.INK)
ax.legend(loc="upper right", handlelength=1.6)
ST.panel_label(ax, "A")

ax = axes[0, 1]
setB = [("CASCADE", ST.MUTED, "cascade"),
        ("C1_TXN_AND", ST.BLUE, "C1 AND"),
        ("C1_TXN_OR", ST.ORANGE, "C1 OR (SUM)"),
        ("C2_MIR_AND", ST.AQUA, "C2 miRNA (AND)")]
for nm, col, lab in setB:
    t, p = curve(nm, "OFF_step")
    p0, pe = p[0], p[-1]
    ax.plot(t, (p - pe) / max(p0 - pe, 1e-9), color=col, lw=1.8, label=lab)
ax.axhline(0.5, color=ST.GRID, lw=0.7)
ax.set_xlim(0, 20); ax.set_ylim(-0.05, 1.05)
ax.set_xlabel("time (target-mRNA lifetimes, $\\tau$)")
ax.set_ylabel("fraction of the OFF step remaining")
ax.set_title("OFF step: the delay is sign-sensitive", loc="left", color=ST.INK)
ax.legend(loc="upper right", handlelength=1.6)
ST.panel_label(ax, "B")

ax = axes[1, 0]
setC = [("CASCADE", ST.MUTED, "cascade"),
        ("C1_TXN_AND", ST.BLUE, "C1 AND"),
        ("I1_MIR_AND", ST.AQUA, "I1 miRNA (AND)"),
        ("COMP_C2_AND", ST.VIOLET, "composite C2 (AND)")]
for nm, col, lab in setC:
    t, ps = curve(nm, "brief_pulse_0.5tau")
    _, pl = curve(nm, "ON_step")
    base = ps[0]
    denom = np.abs(pl - base).max()
    ax.plot(t, np.abs(ps - base) / max(denom, 1e-9), color=col, lw=1.8, label=lab)
ax.axvspan(0, 0.5, color=ST.GRID, alpha=0.8, lw=0)
ax.set_xlim(0, 12); ax.set_ylim(0, 1.02)
ax.set_xlabel("time (target-mRNA lifetimes, $\\tau$)")
ax.set_ylabel("output / full persistent response")
ax.set_title("brief (0.5 $\\tau$) spurious input is filtered", loc="left", color=ST.INK)
ax.legend(loc="upper right", handlelength=1.6)
ST.panel_label(ax, "C")

ax = axes[1, 1]
setD = [("CASCADE", ST.MUTED, "cascade"),
        ("I1_MIR_AND", ST.AQUA, "I1 miRNA"),
        ("C2_MIR_AND", ST.BLUE, "C2 miRNA"),
        ("COMP_C2_AND", ST.VIOLET, "composite C2")]
for nm, col, lab in setD:
    d = DR[DR.topology == nm].sort_values("S")
    v = d.protein_ss.values
    ax.plot(d.S.values, v / max(v.max(), 1e-12), color=col, lw=1.8, label=lab)
ax.set_xscale("log")
ax.set_xlabel("input signal $S$ (dimensionless)")
ax.set_ylabel("steady-state target protein (scaled)")
ax.set_title("steady-state input-output", loc="left", color=ST.INK)
ax.legend(loc="best", handlelength=1.6)
ST.panel_label(ax, "D")
fig.tight_layout(w_pad=2.4, h_pad=2.2)
save(fig, "Fig_dyn1_circuit_dynamics")

# ============================================================ FIGURE 2 ========
S = pd.read_csv(f"{RES}/dynamics_3node_comparison_table.csv")
order = ["CASCADE", "C1_TXN_AND", "C1_TXN_OR", "I1_TXN_AND", "I1_TXN_OR",
         "C2_MIR_AND", "C2_MIR_OR", "C3_MIR_AND", "C3_MIR_OR",
         "C4_MIR_AND", "C4_MIR_OR", "I1_MIR_AND", "I1_MIR_OR",
         "I3_MIR_AND", "I3_MIR_OR", "I3_MIRLED_AND", "I3_MIRLED_OR",
         "COMP_C2_AND", "COMP_C2_OR", "COMP_C4_AND", "COMP_C4_OR",
         "COMP_I1_AND", "COMP_I1_OR", "COMP_I3_AND", "COMP_I3_OR"]
S = S.set_index("topology").reindex([o for o in order if o in set(S.topology)]).reset_index()
ncount = S["n_sign_resolved_cores_in_BRCA_network"].values
ylab = [f"{t}   ({int(c)})" if c > 0 else f"{t}   (–)" for t, c in zip(S.topology, ncount)]

fig, axes = plt.subplots(1, 2, figsize=(7.6, 6.6),
                         gridspec_kw=dict(width_ratios=[1.0, 1.45]))
ax = axes[0]
cols = ["rel_T_ON_median", "rel_T_OFF_median"]
Mrel = np.log2(S[cols].values.astype(float))
vmax = np.nanmax(np.abs(Mrel))
im = ax.imshow(Mrel, cmap=ST.div_cmap(), vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks([0, 1]); ax.set_xticklabels(["ON step", "OFF step"])
ax.set_yticks(range(len(S))); ax.set_yticklabels(ylab, fontsize=6.4)
ax.grid(False)
for i in range(len(S)):
    for j in range(2):
        v = S[cols[j]].values[i]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                color=ST.INK if abs(Mrel[i, j]) < 0.6 * vmax else "white")
ax.set_title("response time relative to the cascade\n(<1 accelerated, >1 delayed)",
             loc="left", color=ST.INK, fontsize=8.5)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.set_label("$\\log_2$ (relative response time)", fontsize=7)
cb.ax.tick_params(labelsize=6.5)
ST.panel_label(ax, "A", dx=-0.62)

ax = axes[1]
fcols = ["pct_pulse", "pct_adaptation_gt_0p5", "pct_FCD",
         "pct_filter_index_gt_0p5", "pct_ultrasensitive_neff_gt2",
         "pct_nonmonotone_dose"]
flab = ["pulse", "adaptation >0.5", "fold-change detection",
        "noise filter index >0.5", "ultrasensitive $n_{eff}$>2",
        "non-monotone dose-response"]
Mf = S[fcols].values.astype(float)
im2 = ax.imshow(Mf, cmap=ST.seq_cmap(), vmin=0, vmax=100, aspect="auto")
ax.set_xticks(range(len(fcols)))
ax.set_xticklabels(flab, fontsize=6.6, rotation=32, ha="right",
                   rotation_mode="anchor")
ax.set_yticks(range(len(S))); ax.set_yticklabels([])
ax.grid(False)
for i in range(len(S)):
    for j in range(len(fcols)):
        ax.text(j, i, f"{Mf[i,j]:.0f}", ha="center", va="center", fontsize=5.8,
                color=ST.INK if Mf[i, j] < 55 else "white")
ax.set_title("% of 1,024 sampled parameter sets showing each behaviour",
             loc="left", color=ST.INK, fontsize=8.5)
cb2 = fig.colorbar(im2, ax=ax, fraction=0.04, pad=0.02)
cb2.set_label("% of parameter space", fontsize=7); cb2.ax.tick_params(labelsize=6.5)
ST.panel_label(ax, "B", dx=-0.06)
fig.text(0.005, 0.005, "row labels: topology (number of sign-resolved cores of that "
                       "architecture in the BRCA network; – = reference circuit, not observed)",
         fontsize=6.2, color=ST.MUTED)
fig.tight_layout(w_pad=1.0, rect=(0, 0.02, 1, 1))
save(fig, "Fig_dyn2_comparison_heatmap")

# ============================================================ FIGURE 3 ========
E = pd.read_csv(f"{RES}/dynamics_3node_metrics_ensemble.csv")
E = E[E.param_set >= 0]
keep = [o for o in order if o in set(E.topology)]
fig, axes = plt.subplots(1, 2, figsize=(7.4, 6.0), sharey=True)
for ax, col, ttl in zip(axes, ["rel_T_ON", "rel_T_OFF"], ["ON step", "OFF step"]):
    data = [np.log2(E.loc[E.topology == t, col].replace([np.inf, -np.inf], np.nan).dropna())
            for t in keep]
    bp = ax.boxplot(data, vert=False, widths=0.62, patch_artist=True, showfliers=False,
                    medianprops=dict(color=ST.INK, lw=1.2),
                    whiskerprops=dict(color=ST.AXIS, lw=0.8),
                    capprops=dict(color=ST.AXIS, lw=0.8))
    for patch, t in zip(bp["boxes"], keep):
        c = (ST.MUTED if t == "CASCADE" else
             ST.BLUE if "TXN" in t else
             ST.VIOLET if t.startswith("COMP") else ST.AQUA)
        patch.set_facecolor(c); patch.set_alpha(0.55)
        patch.set_edgecolor(c); patch.set_linewidth(1.0)
    ax.axvline(0, color=ST.INK, lw=0.9, ls=(0, (4, 3)))
    ax.set_yticks(range(1, len(keep) + 1)); ax.set_yticklabels(keep, fontsize=6.4)
    ax.set_xlabel("$\\log_2$ (response time / cascade)")
    ax.set_title(ttl, loc="left", color=ST.INK)
axes[0].invert_yaxis()
handles = [Line2D([], [], color=c, lw=6, alpha=0.55, label=l) for c, l in
           [(ST.MUTED, "cascade control"), (ST.BLUE, "transcriptional FFL (Alon)"),
            (ST.AQUA, "miRNA-mediated FFL"), (ST.VIOLET, "composite FFL")]]
axes[1].legend(handles=handles, loc="lower right", fontsize=7)
ST.panel_label(axes[0], "A", dx=-0.42); ST.panel_label(axes[1], "b", dx=-0.06)
fig.tight_layout(w_pad=1.2)
save(fig, "Fig_dyn3_response_times")

# ============================================================ FIGURE 4 ========
if os.path.exists(f"{RES}/dynamics_higher_order_fractions.csv"):
    F = pd.read_csv(f"{RES}/dynamics_higher_order_fractions.csv")
    GL = pd.read_csv(f"{RES}/dynamics_higher_order_gain_loss.csv")
    nsets = int(F.n_param_sets.iloc[0])
    mods = ["n3", "n4", "n5", "n6"]
    beh = ["pct_is_bistable", "pct_is_sustained_osc", "pct_is_damped_osc",
           "pct_is_ultrasensitive", "pct_is_memory", "pct_is_noise_rejecting",
           "pct_is_pulse"]
    blab = ["bistability\n(certified)", "sustained\noscillation", "damped\noscillation",
            "ultrasensitivity\n$n_{eff}$>2", "memory\n$T_{1/2}^{OFF}$>5$\\tau$",
            "noise\nrejection", "pulse"]
    fams = list(F.family.unique())
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4))
    ramp = [ST.SEQ[2], ST.SEQ[5], ST.SEQ[8], ST.SEQ[11]]
    for ax, fam in zip(axes, fams):
        sub = F[F.family == fam].set_index("module")
        x = np.arange(len(beh)); w = 0.20
        for k, mod in enumerate(mods):
            if mod not in sub.index:
                continue
            vals = [sub.loc[mod, b] for b in beh]
            ax.bar(x + (k - 1.5) * w, vals, width=w * 0.9, color=ramp[k],
                   edgecolor=ST.SURFACE, linewidth=0.8, label=f"{mod[1]}-node")
            for xi, v in zip(x + (k - 1.5) * w, vals):
                if v > 0.4:
                    ax.text(xi, v + 0.8, f"{v:.0f}", ha="center", va="bottom",
                            fontsize=5.6, color=ST.INK2)
        ax.set_xticks(x); ax.set_xticklabels(blab, fontsize=7)
        ax.set_ylabel(f"% of {nsets:,} parameter sets")
        unres = " | ".join(f"{m}: {sub.loc[m, 'pct_bistability_unresolved']:.1f}%"
                           for m in mods if m in sub.index)
        ax.set_title(f"{fam.replace('_', ' ')}      "
                     f"(steady state not reached within 500$\\tau$ — {unres})",
                     loc="left", color=ST.INK, fontsize=8)
        ax.legend(ncol=4, loc="upper right", fontsize=7, handlelength=1.2)
    ST.panel_label(axes[0], "A"); ST.panel_label(axes[1], "B")
    fig.tight_layout(h_pad=2.0)
    save(fig, "Fig_dyn4_higher_order_fractions")

    g = GL[(GL.metric_type == "binary") & (GL.module.isin(["n4", "n5", "n6"]))].copy()
    fig, axes = plt.subplots(1, len(fams), figsize=(7.6, 4.4), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, fam in zip(axes, fams):
        sub = g[g.family == fam]
        y = np.arange(len(beh)); w = 0.24
        for k, mod in enumerate(["n4", "n5", "n6"]):
            gg = sub[sub.module == mod].set_index("behaviour")
            gain = [gg.loc[b[4:], "pct_gain"] if b[4:] in gg.index else 0 for b in beh]
            loss = [-gg.loc[b[4:], "pct_loss"] if b[4:] in gg.index else 0 for b in beh]
            ax.barh(y + (k - 1) * w, gain, height=w * 0.9, color=ST.DIV_HI,
                    alpha=0.45 + 0.2 * k, edgecolor=ST.SURFACE, linewidth=0.6,
                    label=f"{mod[1]}-node" if ax is axes[0] else None)
            ax.barh(y + (k - 1) * w, loss, height=w * 0.9, color=ST.DIV_LO,
                    alpha=0.45 + 0.2 * k, edgecolor=ST.SURFACE, linewidth=0.6)
        ax.axvline(0, color=ST.INK, lw=0.9)
        ax.set_yticks(y); ax.set_yticklabels(blab, fontsize=7)
        ax.set_xlabel("% of parameter sets\n$\\leftarrow$ lost vs core   |   gained vs core $\\rightarrow$")
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=7, loc="lower left", title="module", title_fontsize=7)
    ST.panel_label(axes[0], "a", dx=-0.45)
    if len(axes) > 1:
        ST.panel_label(axes[1], "b", dx=-0.06)
    fig.tight_layout(w_pad=1.4)
    # Figure 6 of the paper is drawn from the same table by 62_fig6_gain_loss.R (print-size version).
    save(fig, "Fig_dyn5_higher_order_gain_loss")

# ============================================================ FIGURE 6 ========
if os.path.exists(f"{RES}/dynamics_mir29_fit.csv"):
    FIT = pd.read_csv(f"{RES}/dynamics_mir29_fit.csv").set_index("parameter")
    T = pd.read_csv(f"{RES}/dynamics_mir29_tcga_measurements.csv")
    DRm = pd.read_csv(f"{RES}/dynamics_mir29_doseresponse.csv")
    PRED = pd.read_csv(f"{RES}/dynamics_mir29_predictions.csv")
    med_m = float(np.median(T.miR29_family)); med_x = float(np.median(T.NFKB1))
    lam = FIT.loc["lambda", "estimate"]; K = FIT.loc["K", "estimate"]; n = FIT.loc["n", "estimate"]
    c0 = FIT.loc["c0", "estimate"]; a = FIT.loc["a_NFKB1", "estimate"]
    fh = FIT.loc["miR29_fold_to_halve_COL1A1_mRNA", "estimate"]
    fh_lo = FIT.loc["miR29_fold_to_halve_COL1A1_mRNA", "ci_lo"]
    fh_hi = FIT.loc["miR29_fold_to_halve_COL1A1_mRNA", "ci_hi"]

    fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.0))
    ax = axes[0]
    x = T.miR29_family.values - med_m
    yres = T.COL1A1.values - a * (T.NFKB1.values - med_x)
    ax.scatter(2 ** x, yres, s=5, color=ST.BLUE, alpha=0.22, linewidths=0, rasterized=True)
    gx = np.logspace(np.log10(0.2), np.log10(5), 200)
    ax.plot(gx, c0 - np.log2(1 + lam * gx ** n / (K ** n + gx ** n)), color=ST.ORANGE, lw=2.0)
    dlabel(ax, 3.2, c0 - np.log2(1 + lam * 3.2 ** n / (K ** n + 3.2 ** n)) + 0.7,
           "fitted circuit", ST.ORANGE, ha="center")
    ax.set_xscale("log")
    ax.set_xlabel("miR-29 family (fold of cohort median)")
    ax.set_ylabel("$\\log_2$ COL1A1, NFKB1-adjusted")
    ax.set_title(f"TCGA-BRCA, n = {len(T):,}", loc="left", color=ST.INK)
    ST.panel_label(ax, "A", dx=-0.32)

    ax = axes[1]
    ax.plot(DRm.miR29_fold_over_median, DRm.COL1A1_mRNA_rel, color=ST.BLUE, lw=2.0)
    ax.plot(DRm.miR29_fold_over_median, DRm.COL1A1_protein_rel, color=ST.ORANGE, lw=2.0)
    dlabel(ax, 16, np.interp(16, DRm.miR29_fold_over_median, DRm.COL1A1_mRNA_rel),
           "mRNA (TCGA)", ST.BLUE, dy=-0.055)
    dlabel(ax, 16, np.interp(16, DRm.miR29_fold_over_median, DRm.COL1A1_protein_rel),
           "protein (CPTAC)", ST.ORANGE, dy=0.055)
    ax.axhline(0.5, color=ST.AXIS, lw=0.8, ls=(0, (4, 3)))
    if np.isfinite(fh):
        ax.plot([fh], [0.5], marker="o", ms=6, color=ST.INK, zorder=5)
        ax.annotate(f"{fh:.2f}$\\times$ halves COL1A1 mRNA\n(95% CI {fh_lo:.2f}–{fh_hi:.2f})",
                    xy=(fh, 0.5), xytext=(fh * 1.6, 0.76), fontsize=6.6, color=ST.INK,
                    arrowprops=dict(arrowstyle="-", color=ST.MUTED, lw=0.8))
    ax.set_xscale("log"); ax.set_xlim(0.1, 100); ax.set_ylim(0, 1.15)
    ax.set_xlabel("miR-29 family (fold of cohort median)")
    ax.set_ylabel("COL1A1 relative to cohort median")
    ax.set_title("fitted steady-state dose-response", loc="left", color=ST.INK)
    ST.panel_label(ax, "B", dx=-0.32)

    ax = axes[2]
    sub = PRED[(PRED.protein_half_life_h == 46.0) &
               (np.isclose(PRED.miR29_fold_over_cohort_median, fh, rtol=1e-3))]
    cols = {24.0: ST.BLUE, 9.0: ST.ORANGE, 1.5: ST.AQUA}
    for _, r in sub.iterrows():
        F0 = 1 + lam / (K ** n + 1); F1 = 1 + lam * fh ** n / (K ** n + fh ** n)
        tau = r.tau_hours_per_time_unit; gp = r.mRNA_half_life_h / 46.0
        tt = np.linspace(0, 240 / tau, 800)
        y0, yinf = 1 / F0, 1 / F1
        A = (y0 - yinf) * gp / (gp - F1); B = (y0 - yinf) - A
        p = yinf + A * np.exp(-F1 * tt) + B * np.exp(-gp * tt)
        c = cols.get(r.mRNA_half_life_h, ST.MUTED)
        ax.plot(tt * tau, p / y0, color=c, lw=1.8,
                label=f"mRNA $t_{{1/2}}$ {r.mRNA_half_life_h:g} h")
    ax.axhline(0.5, color=ST.AXIS, lw=0.8, ls=(0, (4, 3)))
    ax.set_xlim(0, 240); ax.set_ylim(0.4, 1.03)
    ax.set_xlabel("hours after a %.2f-fold miR-29 step" % fh)
    ax.set_ylabel("COL1A1 protein (relative)")
    ax.legend(fontsize=6.4, loc="upper right")
    ax.set_title("predicted kinetics (protein $t_{1/2}$ 46 h)",
                 loc="left", color=ST.INK, fontsize=8)
    ST.panel_label(ax, "C", dx=-0.32)
    fig.tight_layout(w_pad=2.2)
    save(fig, "Fig_dyn6_mir29_circuit")

# ============================================================ FIGURE 7 ========
hy = f"{RES}/dynamics_composite_hysteresis_curves.csv"
ri = f"{RES}/dynamics_composite_ringing_timecourses.csv"
osc = f"{RES}/dynamics_3node_oscillation_by_titration.csv"
if os.path.exists(hy) and os.path.exists(ri) and os.path.exists(osc):
    H = pd.read_csv(hy); R = pd.read_csv(ri); O = pd.read_csv(osc)
    fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.1))

    ax = axes[0]
    style = {"COMP_C2_AND": (ST.VIOLET, "composite C2 (TF -| miR -| TF)"),
             "C2_MIR_AND": (ST.AQUA, "same core, no reciprocal arm"),
             "C1_TXN_AND": (ST.BLUE, "Alon C1 (transcriptional)")}
    for nm, (c, lab) in style.items():
        d = H[H.topology == nm].sort_values("S")
        if not len(d):
            continue
        ycol_lo = "protein_low_IC" if "protein_low_IC" in d.columns else "protein_up_sweep"
        ycol_hi = "protein_high_IC" if "protein_high_IC" in d.columns else "protein_down_sweep"
        ax.plot(d.S, d[ycol_lo], color=c, lw=1.9, label=lab)
        ax.plot(d.S, d[ycol_hi], color=c, lw=1.9, ls=(0, (3, 2)))
    ax.set_xscale("log")
    ax.set_xlabel("input signal $S$ (dimensionless)")
    ax.set_ylabel("steady-state target protein")
    ax.set_title("two certified steady states at one input", loc="left",
                 color=ST.INK, fontsize=8)
    ax.text(0.30, 0.34, "solid: settled from a low initial state\ndashed: settled from a high one",
            transform=ax.transAxes, fontsize=6.0, color=ST.INK2, va="top")
    ax.legend(loc="upper left", fontsize=6.0, handlelength=1.4,
              bbox_to_anchor=(0.20, 0.62))
    ST.panel_label(ax, "A", dx=-0.32)

    ax = axes[1]
    style2 = {"COMP_I1_AND": (ST.VIOLET, "composite I1 (TF -> miR -| TF)"),
              "I1_MIR_AND": (ST.AQUA, "same core, no reciprocal arm"),
              "I1_TXN_AND": (ST.BLUE, "Alon I1"),
              "CASCADE": (ST.MUTED, "cascade")}
    for nm, (c, lab) in style2.items():
        d = R[R.topology == nm].sort_values("t")
        if not len(d):
            continue
        ax.plot(d.t, d.protein_over_steady_state, color=c, lw=1.8, label=lab)
    ax.axhline(1.0, color=ST.AXIS, lw=0.7, ls=(0, (4, 3)))
    ax.set_xlim(0, 40)
    ax.set_xlabel("time (target-mRNA lifetimes, $\\tau$)")
    ax.set_ylabel("target protein / steady state")
    ax.set_title("ringing: same parameters, one extra edge", loc="left",
                 color=ST.INK, fontsize=8)
    ax.legend(loc="lower right", fontsize=6.3, handlelength=1.4)
    ST.panel_label(ax, "B", dx=-0.32)

    ax = axes[2]
    sel = ["C1_TXN_AND", "I1_TXN_AND", "I1_MIR_AND", "C2_MIR_AND",
           "COMP_I1_AND", "COMP_C2_AND"]
    O2 = O[O.topology.isin(sel)].copy()
    piv = O2.pivot_table(index="topology", columns="titration",
                         values="pct_any_complex_eigenvalue")
    piv = piv.reindex([s for s in sel if s in piv.index])
    y = np.arange(len(piv)); w = 0.36
    for k, cname in enumerate(piv.columns):
        ax.barh(y + (k - 0.5) * w, piv[cname].values, height=w * 0.9,
                color=(ST.ORANGE if "0 (c" in cname else ST.VIOLET),
                edgecolor=ST.SURFACE, linewidth=0.7,
                label=("catalytic $\\theta=0$" if "0 (c" in cname
                       else "stoichiometric $\\theta>0$"))
    ax.set_yticks(y); ax.set_yticklabels(piv.index, fontsize=6.6)
    ax.invert_yaxis()
    ax.set_xlabel("% of parameter sets with a complex\neigenvalue at the operating point")
    ax.set_title("capacity to oscillate at all", loc="left", color=ST.INK, fontsize=8)
    ax.legend(fontsize=6.4, loc="lower right")
    ST.panel_label(ax, "C", dx=-0.36)
    fig.tight_layout(w_pad=2.2)
    save(fig, "Fig_dyn7_composite_only")

# ============================================================ FIGURE 8 ========
tt = f"{RES}/dynamics_titration_threshold.csv"
if os.path.exists(tt):
    TT = pd.read_csv(tt)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax = axes[0]
    ramp = {0.0: ST.MUTED, 1.0: ST.SEQ[4], 5.0: ST.SEQ[7], 20.0: ST.SEQ[9], 100.0: ST.SEQ[12]}
    for th, g in TT[TT.arm == "miRNA_titration"].groupby("theta"):
        g = g.sort_values("beta_target")
        ax.plot(g.beta_target, g.target_ss, color=ramp.get(th, ST.INK), lw=1.7)
        dlabel(ax, g.beta_target.values[-1], g.target_ss.values[-1],
               f"$\\theta$={th:g}", ramp.get(th, ST.INK), ha="left", dx=0.4, fs=6.4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("target transcription rate $\\beta_y$")
    ax.set_ylabel("steady-state target")
    ax.set_title("stoichiometric titration makes a threshold", loc="left",
                 color=ST.INK, fontsize=8.5)
    ST.panel_label(ax, "A", dx=-0.30)

    ax = axes[1]
    for th, g in TT[TT.arm == "miRNA_titration"].groupby("theta"):
        g = g.sort_values("beta_target")
        ax.plot(g.beta_target, g.local_loglog_slope, color=ramp.get(th, ST.INK), lw=1.7)
    for rl, g in TT[TT.arm == "transcriptional_repressor"].groupby("repressor_level"):
        g = g.sort_values("beta_target")
        ax.plot(g.beta_target, g.local_loglog_slope, color=ST.ORANGE, lw=1.4,
                ls=(0, (3, 2)))
    dlabel(ax, 20, 1.0, "transcriptional repressor:\nslope 1 at every level",
           ST.ORANGE, ha="right", dy=-0.35, fs=6.4)
    ax.axhline(1.0, color=ST.AXIS, lw=0.8, ls=(0, (4, 3)))
    ax.set_xscale("log")
    ax.set_xlabel("target transcription rate $\\beta_y$")
    ax.set_ylabel("local $d\\log y / d\\log \\beta_y$")
    ax.set_title("slope > 1 is a threshold no promoter can copy", loc="left",
                 color=ST.INK, fontsize=8.5)
    ST.panel_label(ax, "B", dx=-0.30)
    fig.tight_layout(w_pad=2.0)
    save(fig, "Fig_dyn8_titration_threshold")

print("figures done ->", FIG)
