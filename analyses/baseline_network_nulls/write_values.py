#!/usr/bin/env python3
"""Writes VALUES.md for the baseline network analysis from the tables in this folder: every number with its file, line
(header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

V = ["# VALUES: a baseline network without FFL pre-selection", "",
     "Lines count the header as line 1. Networks: `build_baseline_network.R` (log `runs/build_baseline_network.log`); tests:",
     "`summarise_baseline_nulls.py` from the replicate files in `runs/`.", "",
     "The curated lists from which the network was assembled (1,057 genes, 421 miRNAs) are not available. The substitute keeps the 587 nodes",
     "of the analysed network and takes every edge among them from the primary tables.", "",
     "## 1. Networks (`network_arc_counts.csv`)", ""]
A = rows("network_arc_counts.csv"); hdr = list(A[0][1].keys())
V += ["| line | " + " | ".join(hdr) + " |", "|" + "---|" * (len(hdr) + 1)]
for i, r in A:
    V.append(f"| {i} | " + " | ".join(r[h] for h in hdr) + " |")
V += ["", "Checks (`runs/build_baseline_network.log`):"]
for i, l in enumerate(open("runs/build_baseline_network.log"), 1):
    if "equal to the project's" in l:
        V.append(f"- line {i}: {l.strip()}")
V += ["", "## 2. Overlap with the analysed network (`analysed_network_overlap.csv`)", "",
      "| line | network | edge type | analysed arcs | in the primary tables | baseline arcs |", "|---|---|---|---|---|---|"]
for i, r in rows("analysed_network_overlap.csv"):
    V.append(f"| {i} | {r['network']} | {r['edge_type']} | {r['analysed_arcs']} | {r['in_primary_tables']} | {r.get('baseline_arcs', '')} |")
V += ["", "## 3. Three-node FFL motif test (`null_model_summary.csv`)", "",
      "| line | network | typing | null model | metric | observed | expected (SD) | fold | excess % | p upper | p lower |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("null_model_summary.csv"):
    V.append(f"| {i} | {r['network']} | {r['typing']} | {r['null_model']} | {r['metric']} | {r['observed']} | "
             f"{float(r['expected']):.1f} ({float(r['sd']):.1f}) | {float(r['fold']):.3f} | {float(r['excess_pct']):+.1f} | "
             f"{float(r['p_upper']):.4f} | {float(r['p_lower']):.4f} |")
V += ["", "## 4. Decision (`decision.csv`; persists = p_upper < 0.05 under both the reciprocity-preserving and the",
      "bipartite-layer null; reverses = p_lower < 0.05 under both; vanishes otherwise)", "",
      "| line | network | typing | metric | reciprocity-preserving: excess % (p upper) | bipartite-layer: excess % (p upper) | decision |",
      "|---|---|---|---|---|---|---|"]
for i, r in rows("decision.csv"):
    V.append(f"| {i} | {r['network']} | {r['typing']} | {r['metric']} | {float(r['excess_reciprocity_preserving_pct']):+.1f} "
             f"({float(r['p_upper_reciprocity_preserving']):.4f}) | {float(r['excess_bipartite_layer_pct']):+.1f} "
             f"({float(r['p_upper_bipartite_layer']):.4f}) | {r['decision']} |")
V.append("")
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
