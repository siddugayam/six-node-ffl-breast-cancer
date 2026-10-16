#!/usr/bin/env python3
"""Writes VALUES.md for the reciprocal-handling census from the tables in this folder: every number with its file,
line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

def num(x, f="{:,.0f}"):
    try: return f.format(float(x))
    except ValueError: return x

V = ["# VALUES: reciprocal TF<->miRNA handling in the topological-definition census", "",
     "Lines count the header as line 1. Tables are written by `summarise_reciprocal_variants.py`,",
     "`six_node_pattern_variants.py` and `make_scc_file.py` from the counts in `counts/` (run by",
     "`run_reciprocal_variants.sh`). The baseline is the published census (contracted census graph, 8,003 arcs).", "",
     "## 1. Strongly connected components of the census graph (`scc_structure.csv`)", "",
     "| line | graph | component (smallest vertex) | size | TF | miRNA | Gene | note |", "|---|---|---|---|---|---|---|---|"]
for i, r in rows("scc_structure.csv"):
    V.append(f"| {i} | {r['graph']} | {r['component_smallest_vertex']} | {r['size']} | {r['TF']} | {r['miRNA']} | {r['Gene']} | {r['note']} |")
V += ["", "## 2. Counts, overlap with the baseline, lost and gained (`census_reciprocal_variants.csv`)", "",
      "Columns baseline_count, variant_count, shared, jaccard, lost_fraction (of the baseline), gained_fraction (of the",
      "variant), variant_sd and baseline_sd (n = 7, SD of the 8 seed estimates).", "",
      "| line | n | variant | method | baseline | variant | shared | Jaccard | lost | gained | SD (variant / baseline) |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("census_reciprocal_variants.csv"):
    j = "" if r["jaccard"] in ("", "nan") else f"{float(r['jaccard']):.4f}"
    g = "" if r["gained_fraction"] in ("", "nan") else f"{float(r['gained_fraction']):.4f}"
    sd = f"{r['variant_sd']} / {r['baseline_sd']}" if r["variant_sd"] else ""
    V.append(f"| {i} | {r['n']} | {r['variant']} | {r['method']} | {num(r['baseline_count'])} | {num(r['variant_count'])} | "
             f"{num(r['shared'])} | {j} | {float(r['lost_fraction']):.4f} | {g} | {sd} |")
V += ["", "Notes on the rows:", ""]
for i, r in rows("census_reciprocal_variants.csv"):
    if r["note"]: V.append(f"- line {i}: {r['note']}")
V += ["", "## 3. Baseline modules lost, by first failing condition (`lost_under_module_condensation_first_condition.csv`)", "",
      "| line | n | baseline | lost under M | of which D3 fails first | D4 fails first | lost under U (all fail D2) | lost under S: two nodes share an SCC | lost under S: image not a module |",
      "|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("lost_under_module_condensation_first_condition.csv"):
    V.append(f"| {i} | {r['n']} | {num(r['baseline_modules'])} | {num(r['lost_under_M'])} | {num(r['first_failing_D3'])} | "
             f"{num(r['first_failing_D4'])} | {num(r['lost_under_U_all_by_D2'])} | {num(r['lost_under_S_share_an_SCC'])} | "
             f"{num(r['lost_under_S_image_not_a_module'])} |")
V += ["", "## 4. Composite three-node FFLs under module-level condensation (`composite_three_node_under_module_condensation.csv`)", ""]
for i, r in rows("composite_three_node_under_module_condensation.csv"):
    V.append(f"- line {i}: composite three-node FFLs {r['composite_three_node_FFLs']}; valid under M {r['valid_under_M']}; "
             f"fail D4 {r['fail_D4_under_M']}; fail D3 {r['fail_D3_under_M']}.")
V += ["", "## 5. The six-node composite pattern (`six_node_pattern_variants.csv`, columns instances and node_sets)", "",
      "| line | treatment | instances | node sets | note |", "|---|---|---|---|---|"]
for i, r in rows("six_node_pattern_variants.csv"):
    V.append(f"| {i} | {r['treatment']} | {r['instances']} | {r['node_sets']} | {r['note']} |")
V += ["", "## 6. Independent check at small n (`check_small_n_output.txt`)", "",
      "A separate Python implementation (undirected ESU over all weakly connected subgraphs; its own D2–D4 test) gives the",
      "same counts as the C enumerator:"]
for line in open("check_small_n_output.txt"):
    if line.startswith("n="): V.append(f"- `{line.strip()}`")
V.append("")
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
