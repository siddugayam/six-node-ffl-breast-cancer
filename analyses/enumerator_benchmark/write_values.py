#!/usr/bin/env python3
"""Writes VALUES.md for the benchmark of the topological-definition enumerator from the tables in this folder: every number with its file,
line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

S = {r["measure"]: (i, r) for i, r in rows("benchmark_summary.csv")}
def g(x):
    x = float(x)
    return f"{int(x):,}" if x == int(x) else f"{x:,.1f}" if abs(x) >= 100 else f"{x:.4g}"

def cell(m, fmt=g):
    i, r = S[m]
    return f"{fmt(r['median'])} ({fmt(r['min'])}–{fmt(r['max'])}) [line {i}]"

def reps(m):
    return S[m][1]["replicates"]

V = ["# VALUES: benchmark of the topological-definition enumerator on planted modules", "",
     "Lines count the header as line 1. Tables are written by `benchmark.py`.",
     "`benchmark_summary.csv` gives, per measure (column measure), the number of replicates and the median, min, max",
     "and mean over replicates. Below, each value is the median (min–max) over replicates, with its line in that file.",
     "Per-replicate values are in `benchmark_per_replicate.csv` (one row per replicate, columns as the measures).", "",
     "## 1. the topological definition topologies (`planted_topologies.csv`, column D1_D4_topologies_up_to_isomorphism)", ""]
for i, r in rows("planted_topologies.csv"):
    V.append(f"- line {i}: n = {r['n']}: {int(r['D1_D4_topologies_up_to_isomorphism']):,} topologies")
V += ["", "## 2. Graphs", "",
      f"- background edges (the contracted analysed network, rewired): {cell('background_edges')}",
      f"- noise arcs: {cell('noise_arcs')}",
      f"- arcs of the final graph: {cell('arcs')}",
      f"- replicates: {reps('arcs')}", "",
      "## 3. The topological-definition enumerator", "",
      "| n | planted | intact | reported | recall of intact planted | destroyed planted reported | precision (sample of 1,000) |",
      "|---|---|---|---|---|---|---|"]
for n in (3, 4, 5, 6):
    V.append(f"| {n} | {cell(f'n{n}_planted')} | {cell(f'n{n}_intact')} | {cell(f'n{n}_reported')} | {cell(f'n{n}_recall_intact')} | "
             f"{cell(f'n{n}_destroyed_reported')} | {cell(f'n{n}_precision_sample')} |")
V += ["", "Exhaustive cross-check against the independent Python enumeration:", ""]
for n in (3, 4):
    V.append(f"- n = {n} ({reps(f'n{n}_exhaustive_truth')} replicates): modules {cell(f'n{n}_exhaustive_truth')}; recall "
             f"{cell(f'n{n}_exhaustive_recall')}; precision {cell(f'n{n}_exhaustive_precision')}")
V += ["", "## 4. Relaxations", "",
      "ratio = modules reported under the relaxed test / modules reported under the topological definition; contains = the relaxed set contains",
      "every topological-definition module (1 = true in every replicate); share planted = planted modules reported / all reported.", "",
      "| n | test | reported | ratio | contains the topological definition set | recall of intact planted | share planted |",
      "|---|---|---|---|---|---|---|"]
for n in (3, 4, 5, 6):
    V.append(f"| {n} | the topological definition | {cell(f'n{n}_reported')} | 1 | – | {cell(f'n{n}_recall_intact')} | {cell(f'n{n}_share_planted')} |")
    for mode, lab in (("d3relax", "D3 relaxed"), ("d4relax", "D4 relaxed")):
        V.append(f"| {n} | {lab} | {cell(f'n{n}_{mode}_reported')} | {cell(f'n{n}_{mode}_ratio')} | "
                 f"{cell(f'n{n}_{mode}_contains_strict')} | {cell(f'n{n}_{mode}_recall_intact')} | {cell(f'n{n}_{mode}_share_planted')} |")
same = all(float(S[f"n{n}_d3relax_ratio"][1][k]) == 1 for n in (3, 4, 5, 6) for k in ("min", "max"))
V += ["", f"D3 relaxed, ratio 1 in every replicate at every n: {same}. When true, the D3-relaxed test reports the same modules as the topological definition. topological-definition modules are",
      "acyclic, and in an acyclic subgraph with one source and one sink every vertex lies on a source–sink path, so dropping",
      "that clause of D3 removes nothing.", ""]
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
