#!/usr/bin/env python3
"""Writes VALUES.md for the combined concordance tests from the tables in this folder: every number with its file,
line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

g = lambda x: f"{float(x):.3g}"
KEY_ADJ = {"unadjusted": "unadjusted", "CAF_A": "CAF-A", "CAF_B": "CAF-B",
           "estimate:ABSOLUTE_nontumour": "ABSOLUTE purity", "estimate:ESTIMATE_StromalScore": "ESTIMATE stromal score"}
KEY_SET = ["TF_target_Activation", "miRNA_target_all", "miRNA_target_strong", "miRNA_target_weak", "miRNA_target_predicted_only"]
R = rows("combined_tests_by_adjustment.csv")
V = ["# VALUES: calibrated combination of the stromal-adjusted concordance tests", "",
     "Lines count the header as line 1. Tables are written by `combined_concordance_tests.R` (run log `run.log`).", "",
     "- **Matched pairs per edge per draw:** 10, five draws, seeds 501–505 (`combined_tests_by_adjustment.csv`, column",
     "  pairs_per_edge_per_draw).",
     "- **The 58.1 % and 56.0 % of the TF -> target activation test** come from the CAF-A stromal score: the ESTIMATE stromal",
     "  signature without collagens and network nodes, 123 genes, mean z. They are the CAF_A row below.",
     "- **Columns:**",
     "  - conc, null_pooled;",
     "  - p_mean_two_sided (a), p_twice_mean_two_sided (2 × a, valid under any dependence);",
     "  - p_pooled_two_sided (b);",
     "  - p_cauchy_two_sided, p_hmp_vovk_wang_two_sided (c);",
     "  - p_pooled_one_sided (d), p_mean_one_sided, p_cauchy_one_sided.",
     "- **One-sided:** the alternative is concordance greater than the null, as in the unadjusted test of Fig. S6.", "",
     "## 1. Controls (`control_reproduction.csv`)", "",
     "The published table of 04 is reproduced. The unadjusted, CAF-A and CAF-B rows equal its conc and p columns; the script",
     "asserts this.", "",
     "## 2. Key tests (`combined_tests_by_adjustment.csv`)", "",
     "| line | adjustment | edges | set | n | conc | null (pooled) | (a) mean | 2 × mean | (b) pooled | (c) Cauchy | (c) HM Vovk–Wang | (d) pooled one-sided | mean one-sided | Cauchy one-sided |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in R:
    if r["adjustment"] in KEY_ADJ and r["set"] in KEY_SET:
        V.append(f"| {i} | {KEY_ADJ[r['adjustment']]} | {r['edges']} | {r['set']} | {r['n_edges']} | {float(r['conc']):.4f} | "
                 f"{float(r['null_pooled']):.4f} | {g(r['p_mean_two_sided'])} | {g(r['p_twice_mean_two_sided'])} | "
                 f"{g(r['p_pooled_two_sided'])} | {g(r['p_cauchy_two_sided'])} | {g(r['p_hmp_vovk_wang_two_sided'])} | "
                 f"{g(r['p_pooled_one_sided'])} | {g(r['p_mean_one_sided'])} | {g(r['p_cauchy_one_sided'])} |")
V += ["", "## 3. The 42 distinct estimates (`combined_tests_summary_42_distinct.csv`)", "",
      "The layout is that of `concordance_summary_42_distinct.csv` (columns edges … excess_max, from the arithmetic-mean protocol), with the combined-test",
      "columns added: <method>_p_median, _p_min, _p_max, _n_below_0.05. The reference row for all edges, miRNA_target_weak",
      "is reproduced to 15 significant digits.", "",
      "| line | edges | set | p median (a) | n p < 0.05: (a) | 2 × mean | (b) pooled | (c) Cauchy | (c) HM Vovk–Wang | (d) one-sided pooled | above null |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("combined_tests_summary_42_distinct.csv"):
    V.append(f"| {i} | {r['edges']} | {r['set']} | {g(r['p_median'])} | {r['n_p_below_0.05']} | {r['twice_mean_n_below_0.05']} | "
             f"{r['pooled_two_sided_n_below_0.05']} | {r['cauchy_two_sided_n_below_0.05']} | "
             f"{r['harmonic_mean_vovk_wang_n_below_0.05']} | {r['pooled_one_sided_n_below_0.05']} | {r['n_above_null']} |")
C = rows("conclusion_changes.csv")
two = ["twice_mean", "pooled_two_sided", "cauchy_two_sided", "harmonic_mean_vovk_wang"]
V += ["", "## 4. Does any conclusion change? (`conclusion_changes.csv`, columns sig_<method> and changes_vs_mean)", "",
      "Significance is read at p < 0.05. Relative to the arithmetic mean (a), for the key tests:", ""]
for i, r in C:
    if r["adjustment"] in KEY_ADJ and r["set"] in KEY_SET and r["changes_vs_mean"]:
        V.append(f"- line {i}: {KEY_ADJ[r['adjustment']]}, {r['edges']}, {r['set']}: (a) significant = "
                 f"{r['sig_mean_two_sided']}; differs for: {r['changes_vs_mean']}")
n2 = sum(1 for i, r in C if any(m in r["changes_vs_mean"].split(";") for m in two))
nall = len(C)
V += ["", f"Over all {nall} rows (45 adjustments: unadjusted, CAF-A, CAF-B and the 42 estimates; × 8 sets × 2 edge selections): the call changes",
      f"under at least one two-sided method in {n2} rows. Under the one-sided tests it also changes wherever the concordance",
      "lies below the null, because there the one-sided p is large.", ""]
V += ["Per method, over the same rows (gained = significant only under the method; lost = significant only under (a)):", "",
      "| method | rows that change | gained | lost |", "|---|---|---|---|"]
for m in ["twice_mean", "pooled_two_sided", "cauchy_two_sided", "harmonic_mean_vovk_wang", "pooled_one_sided",
          "mean_one_sided", "cauchy_one_sided"]:
    ch = [r for i, r in C if r["sig_" + m] != r["sig_mean_two_sided"]]
    gn = sum(1 for r in ch if r["sig_" + m] == "TRUE")
    V.append(f"| {m} | {len(ch)} | {gn} | {len(ch) - gn} |")
V.append("")
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
