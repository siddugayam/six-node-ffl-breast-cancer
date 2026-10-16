#!/usr/bin/env python3
"""Writes VALUES.md for the tier- and sign-aware null models from the tables in this folder: every number with its
file, line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

V = ["# VALUES: tier- and sign-aware null models for the three-node FFL test", "",
     "Lines count the header as line 1. Tables are written by `summarise_tier_sign_nulls.py` from the replicate files in",
     "`runs/` (written by `motif_null_tier_sign.c`, run by `run_tier_sign_nulls.sh`). 1,000 randomised networks per run.", "",
     "## 1. Controls (`control_reproduction.csv`, column identical)", ""]
for i, r in rows("control_reproduction.csv"):
    V.append(f"- line {i}: {r['graph']}, {r['typing']}, {r['null']}: identical to `{r['published_file']}` = {r['identical']} "
             f"({r['replicate_lines']} lines)")
V += ["", "## 2. Arcs per tier and sign (`arc_attribute_counts.csv`, column arcs)", "",
      "| line | graph | typing | attribute | value | arc class | arcs |", "|---|---|---|---|---|---|---|"]
for i, r in rows("arc_attribute_counts.csv"):
    V.append(f"| {i} | {r['graph']} | {r['typing']} | {r['attribute']} | {r['value']} | {r['arc_class']} | {r['arcs']} |")
V += ["", "Tiers are counted over all miRNA -> {Gene, TF} arcs. The reciprocal arcs keep their tier, but they are frozen and",
      "never swapped.", "",
      "## 3. Every run (`null_model_summary.csv`)", "",
      "Columns observed, expected (mean of 1,000), sd, fold, excess_pct, p_upper, p_lower.", "",
      "| line | graph | typing | null | metric | observed | expected (SD) | fold | excess % | p upper | p lower |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("null_model_summary.csv"):
    V.append(f"| {i} | {r['graph']} | {r['typing']} | {r['null']} | {r['metric']} | {r['observed']} | "
             f"{float(r['expected']):.1f} ({float(r['sd']):.1f}) | {float(r['fold']):.3f} | {float(r['excess_pct']):+.1f} | "
             f"{float(r['p_upper']):.4f} | {float(r['p_lower']):.4f} |")
V += ["", "## 4. Waterfall of the residual excess (`waterfall.csv`, columns excess_pct and change_pp)", "",
      "| line | graph | typing | metric | step | excess % | change (points) |", "|---|---|---|---|---|---|---|"]
for i, r in rows("waterfall.csv"):
    ch = "" if r["change_pp"] in ("", "nan") else f"{float(r['change_pp']):+.1f}"
    V.append(f"| {i} | {r['graph']} | {r['typing']} | {r['metric']} | {r['step']} | {float(r['excess_pct']):+.1f} | {ch} |")
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
