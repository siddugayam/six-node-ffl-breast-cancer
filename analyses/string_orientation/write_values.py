#!/usr/bin/env python3
"""Writes VALUES.md for the STRING orientation analysis from the tables in this folder: every number with its file,
line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

def n0(x): return f"{float(x):,.0f}" if "." not in str(x) or float(x) > 1000 else f"{float(x):,.1f}"

V = ["# VALUES: orientation of the STRING associations in the census graph", "",
     "Lines count the header as line 1. Tables are written by `string_orientation_graphs.py` (graphs) and",
     "`summarise_string_orientation.py` (counts from `counts/`, run by `run_string_orientation.sh`). The census graph",
     "built here reproduces `census_noexemplar_fine.txt` arc for arc; the published counts at n = 3–6 and the published",
     "n = 7 samples seed by seed are reproduced (the script stops otherwise).", "",
     "## (a) Random orientation of the 941 associations (`string_bootstrap_summary.csv`)", "",
     "Columns replicates, published_count, median, p2_5, p97_5, mean, sd, cv, published_percentile (where the published",
     "count lies in the bootstrap distribution), recovery_median / p2_5 / p97_5 (share of the published modules found",
     "again), jaccard_median / p2_5 / p97_5.", "",
     "| line | n | replicates | published | median | 95 % interval | mean | SD | CV | published percentile | recovery median (95 %) | Jaccard median (95 %) |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("string_bootstrap_summary.csv"):
    V.append(f"| {i} | {r['n']} | {r['replicates']} | {int(r['published_count']):,} | {float(r['median']):,.0f} | "
             f"{float(r['p2_5']):,.0f}–{float(r['p97_5']):,.0f} | {float(r['mean']):,.1f} | {float(r['sd']):,.1f} | "
             f"{float(r['cv']):.4f} | {float(r['published_percentile']):.1f} | {float(r['recovery_median']):.4f} "
             f"({float(r['recovery_p2_5']):.4f}–{float(r['recovery_p97_5']):.4f}) | {float(r['jaccard_median']):.4f} "
             f"({float(r['jaccard_p2_5']):.4f}–{float(r['jaccard_p97_5']):.4f}) |")
G = rows("string_bootstrap_graphs.csv")
sa = [int(r["string_arcs_added"]) for _, r in G]
V += ["", f"STRING arcs entering the graph per replicate (`string_bootstrap_graphs.csv`, column string_arcs_added): "
      f"{min(sa)}–{max(sa)} (the listed orientation gives 833). Per-replicate counts: `string_bootstrap_counts.csv`.", "",
      "## (b) Gene–gene layer variants (`gene_gene_layer_variants.csv`, column modules; sd for n = 7)", "",
      "| line | graph | n | method | modules | SD |", "|---|---|---|---|---|---|"]
for i, r in rows("gene_gene_layer_variants.csv"):
    V.append(f"| {i} | {r['description']} | {r['n']} | {r['method']} | {n0(r['modules'])} | {r['sd']} |")
V += ["", "## (c) Modules that traverse no STRING arc (`orientation_robust_core.csv`)", "",
      "no_string_arc: no STRING arc in the module's induced subgraph. no_string_pair: no two of its vertices form",
      "a STRING association, so the module is unchanged under every orientation. n = 7: mean over the 8 RAND-ESU seeds.", "",
      "| line | n | method | modules | no STRING arc (fraction) | no STRING pair (fraction) |", "|---|---|---|---|---|---|"]
for i, r in rows("orientation_robust_core.csv"):
    V.append(f"| {i} | {r['n']} | {r['method']} | {n0(r['modules'])} | {n0(r['no_string_arc'])} "
             f"({float(r['fraction_no_string_arc']):.4f}) | {n0(r['no_string_pair'])} ({float(r['fraction_no_string_pair']):.4f}) |")
V.append("")
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
