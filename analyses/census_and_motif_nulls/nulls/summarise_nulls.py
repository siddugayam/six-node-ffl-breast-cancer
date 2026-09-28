#!/usr/bin/env python3
"""Summarise the null runs with the formulas of scripts/04_motif_significance/v2_analyse_null.py (Z, fold,
p = (1 + #{random >= obs}) / (R + 1)).  Writes nulls/motif_nulls_all_graphs.csv."""
import os, csv
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
NAMES = {"A": "NULL-A (class+degree, reciprocity destroyed)",
         "B": "NULL-B (as A but reciprocal TF<->miRNA pairs frozen)",
         "C": "NULL-C (B + curveball on bipartite miRNA->target & TF->target layers)"}
GRAPHS = {"dep": "deposited network, 6,859 arcs (incl. 30 legacy miRNA-miRNA)",
          "dep_nolegacy": "deposited network minus the 30 legacy miRNA-miRNA edges, 6,829 arcs",
          "pub": "census graph (Table 3), 9,254 arcs",
          "pub_nolegacy": "census graph minus 28 legacy arcs, 9,226 arcs",
          "pub_nostring": "census graph minus 833 STRING arcs, 8,421 arcs",
          "pub_nolegacy_nostring": "census graph minus legacy and STRING arcs, 8,393 arcs"}
rows = []
for g in GRAPHS:
    for model in "ABC":
        dat = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, f"null_{model}_{g}.tsv"))]
        hdr = dat[0]; obs = dict(zip(hdr[1:], map(int, dat[1][1:])))
        arr = np.array([[int(x) for x in r[1:]] for r in dat[2:]], dtype=float); R = arr.shape[0]
        d = {c: arr[:, i] for i, c in enumerate(hdr[1:])}
        d["typed"] = d["comp"] + d["tf"] + d["mir"]; obs["typed"] = obs["comp"] + obs["tf"] + obs["mir"]
        for metric, label in (("comp", "Composite-FFL (once per reciprocal pair)"), ("tf", "TF-FFL"),
                              ("mir", "miRNA-FFL"), ("typed", "Typed total (Comp+TF+miRNA)"),
                              ("total", "All 3-node FFL modules"), ("recip", "reciprocal TF<->miRNA pairs")):
            v = d[metric]; o = obs[metric]; m, s = v.mean(), v.std(ddof=1)
            rows.append(dict(graph=g, graph_desc=GRAPHS[g], null=model, null_desc=NAMES[model], n_randomisations=R,
                             metric=label, observed=o, random_mean=round(m, 2), random_sd=round(s, 2),
                             Z=round((o - m) / s, 3) if s > 0 else float("nan"), fold=round(o / m, 4) if m > 0 else float("nan"),
                             excess_pct=round(100 * (o / m - 1), 1) if m > 0 else float("nan"),
                             p_enrich=(1 + int((v >= o).sum())) / (R + 1), p_deplete=(1 + int((v <= o).sum())) / (R + 1),
                             direction=("enriched" if o > m else "depleted")))
with open(os.path.join(D, "motif_nulls_all_graphs.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    if r["metric"] in ("All 3-node FFL modules", "Composite-FFL (once per reciprocal pair)", "reciprocal TF<->miRNA pairs"):
        print(f"{r['graph']:22s} {r['null']} {r['metric'][:30]:30s} obs {r['observed']:>6} rand {r['random_mean']:>9.1f} ± {r['random_sd']:>6.1f}  "
              f"fold {r['fold']:.3f} ({r['excess_pct']:+.1f} %)  Z {r['Z']:+.2f}  p_enr {r['p_enrich']:.4f}  p_dep {r['p_deplete']:.4f}")
