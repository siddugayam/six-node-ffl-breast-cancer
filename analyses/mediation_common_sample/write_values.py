#!/usr/bin/env python3
"""Writes VALUES.md for the common-sample mediation from the tables in this folder: every number with its file, line
(header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

f = lambda x: f"{float(x):+.3f}"
V = ["# VALUES: stromal mediation of the twelve paths (seven in Table 3) on a common sample of tumours", "",
     "Lines count the header as line 1. Tables are written by `mediation_common_sample.R` (run log `run.log`).", "",
     "## 1. Control (`control_reproduction.csv`)", ""]
for i, r in rows("control_reproduction.csv"):
    V.append(f"- line {i}: {r['models_compared']} of {r['models_expected']} models (12 paths × 42 estimates, each on its own "
             f"tumours) reproduce the published rows: n equal {r['n_equal']}; largest relative difference in total "
             f"{float(r['max_rel_total']):.1e}, ACME {float(r['max_rel_ACME']):.1e}, ADE {float(r['max_rel_ADE']):.1e}, "
             f"CI bounds {float(r['max_rel_CI']):.1e}; reproduced {r['reproduced']}.")
S = rows("mediation_common_sample_summary.csv")
V += ["", "## 2. Per path: common sample against the published runs (`mediation_common_sample_summary.csv`)", "",
      f"Common sample: {S[0][1]['n_common']} tumours with both assays and all 42 distinct estimates (column n_common). The",
      "total effect is the same for every estimate (checked by the script).", "",
      "| line | path | in Table 3 (published inconsistent there) | total (common) | ACME median (range) | ADE median (range) | proportion mediated median (range) | inconsistent of 42 | of which ADE CI excludes 0 | Sobel p < 0.05 | published: n range | published total range | published inconsistent | published inconsistent, ADE CI excl. 0 |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in S:
    t3 = f"yes ({int(float(r['withexemplar_inconsistent']))})" if r['in_withexemplar'] == 'TRUE' else "no"
    V.append(f"| {i} | {r['x']} -> {r['y']} | {t3} | {f(r['total_common'])} | {f(r['ACME_median'])} ({f(r['ACME_min'])} to {f(r['ACME_max'])}) | "
             f"{f(r['ADE_median'])} ({f(r['ADE_min'])} to {f(r['ADE_max'])}) | {float(r['propmed_median']):.3f} "
             f"({float(r['propmed_min']):.3f} to {float(r['propmed_max']):.3f}) | {r['n_inconsistent']} | "
             f"{r['n_inconsistent_ADE_CI_excl0']} | {r['n_sobel_p05']} | {r['published_n_range']} | {r['published_total_range']} | "
             f"{r['published_n_inconsistent']} | {r['published_n_inconsistent_ADE_CI_excl0']} |")
V += ["", "Columns: total_common; ACME_median, ACME_min, ACME_max; ADE_*; propmed_*; n_inconsistent;",
      "n_inconsistent_ADE_CI_excl0; n_sobel_p05; and published_* (the same 42 estimates in",
      "results/systems_and_dynamics/deconv_mediation_extended.csv).", "",
      "The proportion mediated is ACME / total. It is unstable where the total is near 0.", "",
      "The script checks that the published inconsistent counts of Table 3's seven rows are 23, 32, 38, 22, 6, 6 and 3",
      "(columns published_n_inconsistent and withexemplar_inconsistent).", "",
      "## 3. Every path × estimate (`mediation_common_sample_per_estimate.csv`)", "",
      f"{len(rows('mediation_common_sample_per_estimate.csv'))} rows, with columns x, y, mediator, n, total, ACME (acme_lo, acme_hi, acme_p),",
      "ADE (ade_lo, ade_hi, ade_p), prop_mediated (prop_lo, prop_hi), sobel_p, inconsistent and ADE_CI_excludes_zero.", "",
      "The three exact repeats left out are MCPcounter_Fibroblasts_recomputed (= MCPcounter_Fibroblasts), NNLS_Wu641_full_CAFs",
      "(= NNLS_Wu_CAFs) and QPROG_Wu641_full_CAFs (= QPROG_Wu_CAFs); equality is checked in `run.log`.", ""]
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
