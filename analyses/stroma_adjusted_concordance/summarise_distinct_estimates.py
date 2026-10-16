#!/usr/bin/env python3
"""The 45 stromal estimates contain three runs that reproduce other estimates exactly (MCPcounter_Fibroblasts_recomputed =
MCPcounter_Fibroblasts, NNLS_Wu641_full_CAFs = NNLS_Wu_CAFs, QPROG_Wu641_full_CAFs = QPROG_Wu_CAFs). This script checks that
identity in concordance_adjusted_by_estimate.csv and summarises the 42 distinct estimates per edge set, as concordance_summary.csv
does for all 45 (median and range of the mean p, the estimates with mean p < 0.05, and the estimates with concordance above the
matched null). Output: concordance_summary_42_distinct.csv. Run from this folder."""
import csv, statistics as st
from collections import defaultdict

DUP = {"MCPcounter_Fibroblasts_recomputed": "MCPcounter_Fibroblasts", "NNLS_Wu641_full_CAFs": "NNLS_Wu_CAFs",
       "QPROG_Wu641_full_CAFs": "QPROG_Wu_CAFs"}
rows = [r for r in csv.DictReader(open("concordance_adjusted_by_estimate.csv", encoding="utf8")) if r["analysis"] == "adjusted"]
by = defaultdict(dict)
for r in rows:
    by[(r["edges"], r["set"])][r["estimate"]] = r
for k, d in by.items():                                    # identity check
    for dup, orig in DUP.items():
        a, b = d[dup], d[orig]
        assert all(abs(float(a[c]) - float(b[c])) < 1e-12 for c in ("conc", "null", "p_mean")) and a["n_tumours"] == b["n_tumours"], (k, dup)
out = []
for (edges, s), d in sorted(by.items()):
    ests = {e: r for e, r in d.items() if e not in DUP}
    assert len(ests) == 42
    p = [float(r["p_mean"]) for r in ests.values()]
    ex = [float(r["conc"]) - float(r["null"]) for r in ests.values()]
    out.append({"edges": edges, "set": s, "n_estimates": len(ests), "p_median": st.median(p), "p_min": min(p), "p_max": max(p),
                "n_p_below_0.05": sum(x < 0.05 for x in p), "n_above_null": sum(x > 0 for x in ex),
                "excess_median": st.median(ex), "excess_min": min(ex), "excess_max": max(ex)})
w = csv.DictWriter(open("concordance_summary_42_distinct.csv", "w", encoding="utf8", newline=""), fieldnames=list(out[0]))
w.writeheader(); w.writerows(out)
for o in out:
    print(o["edges"], o["set"], "p<0.05:", o["n_p_below_0.05"], "above null:", o["n_above_null"], "median p %.4g" % o["p_median"], "range %.3g-%.3g" % (o["p_min"], o["p_max"]))
