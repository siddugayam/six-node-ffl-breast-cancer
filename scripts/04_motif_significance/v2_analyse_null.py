#!/usr/bin/env python3
import sys, os, csv, math
import numpy as np

RES = "/path/to/revision/results/v2"
rows = []
NAMES = {"A": "NULL-A (class+degree, reciprocity destroyed)",
         "B": "NULL-B (as A but reciprocal TF<->miRNA pairs frozen)",
         "C": "NULL-C (B + curveball on bipartite miRNA->target & TF->target layers)"}

for graph in ("pub", "dep"):
    for model in ("A", "B", "C"):
        path = os.path.join(RES, "null_%s_%s.tsv" % (model, graph))
        if not os.path.exists(path):
            continue
        dat = [l.rstrip("\n").split("\t") for l in open(path)]
        hdr = dat[0]
        obs = dict(zip(hdr[1:], map(int, dat[1][1:])))
        arr = np.array([[int(x) for x in r[1:]] for r in dat[2:]], dtype=float)
        R = arr.shape[0]
        cols = hdr[1:]
        d = {c: arr[:, i] for i, c in enumerate(cols)}
        d["typed"] = d["comp"] + d["tf"] + d["mir"]
        obs["typed"] = obs["comp"] + obs["tf"] + obs["mir"]
        d["comp_x2"] = 2 * d["comp"]
        obs["comp_x2"] = 2 * obs["comp"]
        for metric, label in (("comp", "Composite-FFL (once per reciprocal pair)"),
                              ("comp_x2", "Composite-FFL (twice, once per arc)"),
                              ("tf", "TF-FFL"), ("mir", "miRNA-FFL"),
                              ("typed", "Typed total (Comp+TF+miRNA)"),
                              ("total", "All 3-node FFL modules"),
                              ("recip", "reciprocal TF<->miRNA pairs")):
            v = d[metric]; o = obs[metric]
            m, s = v.mean(), v.std(ddof=1)
            z = (o - m) / s if s > 0 else float("nan")
            fold = o / m if m > 0 else float("nan")
            p_hi = (1 + int((v >= o).sum())) / (R + 1)
            p_lo = (1 + int((v <= o).sum())) / (R + 1)
            rows.append(dict(graph=graph, null=model, null_desc=NAMES[model], n_randomisations=R,
                             metric=label, observed=o, random_mean=round(m, 2),
                             random_sd=round(s, 2), Z=round(z, 3), fold=round(fold, 4),
                             p_enrich=p_hi, p_deplete=p_lo,
                             direction=("enriched" if o > m else "depleted")))

out = os.path.join(RES, "verify_motif.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows: w.writerow(r)
print("wrote", out, len(rows), "rows")
for r in rows:
    if r["graph"] == "pub":
        print("%-3s %-42s obs=%8d rand=%10.1f+-%7.1f Z=%+8.2f fold=%6.3f p_enr=%.4f p_dep=%.4f  R=%d"
              % (r["null"], r["metric"], r["observed"], r["random_mean"], r["random_sd"],
                 r["Z"], r["fold"], r["p_enrich"], r["p_deplete"], r["n_randomisations"]))
