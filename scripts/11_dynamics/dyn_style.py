#!/usr/bin/env python3
"""Shared figure style. Palette validated with the dataviz validator:
   node scripts/validate_palette.js "#2a78d6,#eb6834,#1baf7a,#4a3aa7" --mode light
   -> all checks PASS (aqua carries a contrast WARN, relieved by direct labels)."""
import matplotlib as mpl
import matplotlib.pyplot as plt

SURFACE   = "#fcfcfb"
INK       = "#0b0b0b"
INK2      = "#52514e"
MUTED     = "#898781"
GRID      = "#e1e0d9"
AXIS      = "#c3c2b7"

# categorical slots, fixed order (never cycled)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = C

# sequential blue ramp, light -> dark
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
       "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
# diverging blue <-> red, neutral gray midpoint
DIV_LO, DIV_MID, DIV_HI = "#2a78d6", "#f0efec", "#e34948"


def apply():
    mpl.rcParams.update({
        "font.family": ["Liberation Sans", "DejaVu Sans", "sans-serif"],
        "font.size": 8,
        "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
        "axes.labelcolor": INK2, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.grid": True, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.6, "lines.markersize": 4,
        "legend.frameon": False,
        "figure.dpi": 130, "savefig.dpi": 400, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel_label(ax, s, dx=-0.16, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="top", ha="left", color=INK)


def seq_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("seqblue", SEQ)


def div_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("divbr", [DIV_LO, DIV_MID, DIV_HI])
