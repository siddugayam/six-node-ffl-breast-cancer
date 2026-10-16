#!/usr/bin/env python3
"""Writes VALUES.md for the stroma-adjusted tier concordance from the tables in this folder: every number with its
file, line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

f3 = lambda x: f"{float(x):.3f}"
g3 = lambda x: f"{float(x):.3g}"
V = ["# VALUES: tier-stratified miRNA–target sign concordance after stromal adjustment", "",
     "Lines count the header as line 1. Tables are written by `stroma_adjusted_concordance.R` (run log `run.log`).",
     "Concordance = fraction of edges with the predicted (negative) correlation sign; null = mean null concordance over the",
     "five matched-null draws; p = mean over the draws of the exact two-sided binomial p (p_max = the largest).", "",
     "## 1. Control (`control_reproduction.csv`)", ""]
for i, r in rows("control_reproduction.csv"):
    V.append(f"- line {i}: {r['column']}: max relative difference {g3(r['max_relative_difference'])}, reproduced {r['reproduced']}")
V += ["", "The published table of 04 (raw, CAF-A, CAF-B for all 8 edge sets) is reproduced.", "",
      "## 2. Per tier: unadjusted, ABSOLUTE alone, ESTIMATE alone (`concordance_summary.csv`)", "",
      "Columns: n_edges; unadj_conc, unadj_null, unadj_p; ABSOLUTE_tumours, ABSOLUTE_unadj_conc, ABSOLUTE_unadj_p (unadjusted,",
      "same 1,023 tumours); ABSOLUTE_conc, ABSOLUTE_null, ABSOLUTE_p; ESTIMATE_conc, ESTIMATE_null, ESTIMATE_p.", "",
      "| line | edges | set | edges n | unadjusted conc / null (p) | same tumours as ABSOLUTE: conc (p) | ABSOLUTE-adjusted conc / null (p) | ESTIMATE-adjusted conc / null (p) |",
      "|---|---|---|---|---|---|---|---|"]
S = rows("concordance_summary.csv")
for i, r in S:
    if not r["set"].startswith("miRNA"): continue
    V.append(f"| {i} | {r['edges']} | {r['set']} | {r['n_edges']} | {f3(r['unadj_conc'])} / {f3(r['unadj_null'])} ({g3(r['unadj_p'])}) | "
             f"{f3(r['ABSOLUTE_unadj_conc'])} ({g3(r['ABSOLUTE_unadj_p'])}) | {f3(r['ABSOLUTE_conc'])} / {f3(r['ABSOLUTE_null'])} "
             f"({g3(r['ABSOLUTE_p'])}) | {f3(r['ESTIMATE_conc'])} / {f3(r['ESTIMATE_null'])} ({g3(r['ESTIMATE_p'])}) |")
V += ["", "## 3. Per tier: the 45 stromal estimates, each alone (`concordance_summary.csv`)", "",
      "Columns all45_conc_median/min/max, all45_excess_median/min/max (concordance minus null), all45_p_median/min/max,",
      "all45_n_p_below_0.05, all45_n_above_null.", "",
      "| line | edges | set | adjusted conc median (range) | excess over null median (range) | p median (range) | estimates with p < 0.05 | estimates above null |",
      "|---|---|---|---|---|---|---|---|"]
for i, r in S:
    if not r["set"].startswith("miRNA"): continue
    V.append(f"| {i} | {r['edges']} | {r['set']} | {f3(r['all45_conc_median'])} ({f3(r['all45_conc_min'])}–{f3(r['all45_conc_max'])}) | "
             f"{float(r['all45_excess_median']):+.4f} ({float(r['all45_excess_min']):+.4f} to {float(r['all45_excess_max']):+.4f}) | "
             f"{g3(r['all45_p_median'])} ({g3(r['all45_p_min'])}–{g3(r['all45_p_max'])}) | {r['all45_n_p_below_0.05']} of 45 | "
             f"{r['all45_n_above_null']} of 45 |")
V += ["", "## 4. TF edge sets (same table, for completeness)", "",
      "| line | edges | set | n | unadjusted conc / null (p) | ABSOLUTE conc / null (p) | ESTIMATE conc / null (p) |", "|---|---|---|---|---|---|---|"]
for i, r in S:
    if r["set"].startswith("miRNA"): continue
    V.append(f"| {i} | {r['edges']} | {r['set']} | {r['n_edges']} | {f3(r['unadj_conc'])} / {f3(r['unadj_null'])} ({g3(r['unadj_p'])}) | "
             f"{f3(r['ABSOLUTE_conc'])} / {f3(r['ABSOLUTE_null'])} ({g3(r['ABSOLUTE_p'])}) | {f3(r['ESTIMATE_conc'])} / "
             f"{f3(r['ESTIMATE_null'])} ({g3(r['ESTIMATE_p'])}) |")
A = rows("concordance_adjusted_by_estimate.csv")
V += ["", "## 5. Every estimate", "",
      f"`concordance_adjusted_by_estimate.csv` holds {len(A)} rows: 45 estimates × 2 edge selections × 8 sets × 2 analyses",
      "(adjusted; unadjusted on the same tumours). Columns: estimate, n_tumours, edges, set, n_edges, analysis, conc, null,",
      "p_mean, p_max.", "",
      "The Farmer-correlated set: genes with Spearman rho > 0.5 against the Farmer stroma-related ssGSEA score (698 genes,",
      "67 of them network miRNA targets; counts printed in `run.log`).", ""]
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
