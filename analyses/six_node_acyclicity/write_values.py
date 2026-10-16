#!/usr/bin/env python3
"""Writes VALUES.md for the six-node pattern and acyclicity analysis from the tables in this folder: every number with
its file, line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

V = ["# VALUES: the six-node pattern and acyclicity", "",
     "Lines count the header as line 1. Part (a) is written by `pattern_first_failing.py`; part (b) by `summarise_d2prime.py`",
     "from the counts of `six_node_d2prime_nulls.c` (run by `run_d2prime_nulls.sh`).", "",
     "## (a) The 2,498 node sets of the six-node composite pattern against the topological definition (`pattern_first_failing_condition.csv`)", "",
     "Graph: the census graph, contracted as in the census. Columns first_failing_condition, cycle_class,",
     "with_D2_waived_literal (D3–D4 on the module as it is), with_D2_waived_condensed (D3–D4 with its cycles condensed),",
     "node_sets.", "",
     "| line | first failing | cycles present | D2 waived: literal | D2 waived: condensed | node sets |", "|---|---|---|---|---|---|"]
tot = 0
for i, r in rows("pattern_first_failing_condition.csv"):
    tot += int(r["node_sets"])
    V.append(f"| {i} | {r['first_failing_condition']} | {r['cycle_class']} | {r['with_D2_waived_literal']} | "
             f"{r['with_D2_waived_condensed']} | {r['node_sets']} |")
V += ["", f"Total node sets: {tot}. Every set holds both a co-transcribed miRNA pair and a reciprocal TF<->TF pair; the",
      "composite pattern requires the TF1<->TF2 link and the 10 kb miRNA pair.", "",
      "## (b) D2′ six-node modules (composite core plus a second TF)", "",
      "### Observed (`d2prime_counts_observed.csv`, column value)", ""]
for i, r in rows("d2prime_counts_observed.csv"):
    V.append(f"- line {i}: {r['quantity']}: {int(r['value']):,}")
V += ["", "### Against NULL-C and NULL-L (`d2prime_null_summary.csv`)", "",
      "Columns observed, null_mean, null_sd, fold, p_upper, p_lower, replicates.", "",
      "| line | null | count | observed | null mean (SD) | fold | p upper | p lower | replicates |", "|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("d2prime_null_summary.csv"):
    V.append(f"| {i} | {r['null']} | {r['count']} | {int(float(r['observed'])):,} | {float(r['null_mean']):,.1f} "
             f"({float(r['null_sd']):,.1f}) | {float(r['fold']):.3f} | {float(r['p_upper']):.4f} | {float(r['p_lower']):.4f} | {r['replicates']} |")
V += ["", "### Decision (`d2prime_decision.csv`, column decision; rule: p_upper < 0.05 under both NULL-C and NULL-L)", ""]
for i, r in rows("d2prime_decision.csv"):
    V.append(f"- line {i}: {r['count']}: {r['decision']} (NULL-C fold {float(r['NULL_C_fold']):.3f}, p_upper "
             f"{float(r['NULL_C_p_upper']):.4f}; NULL-L fold {float(r['NULL_L_fold']):.3f}, p_upper {float(r['NULL_L_p_upper']):.4f}); "
             f"rule applies: {r['rule_applies']}")
V += ["", "### Replicate check (`check_replicates_against_s1.txt`)", ""]
V += [f"- {l.strip()}" for l in open("check_replicates_against_s1.txt") if "check" in l or "identical" in l]
V += ["", "The number of replicates (1,000 per model) follows the rule in SETTINGS.md (section \"Run decisions\").", ""]
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
