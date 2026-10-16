#!/usr/bin/env python3
"""Writes VALUES.md for the multivariable survival analysis from the tables in this folder: every number with its file,
line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        r = list(csv.DictReader(fh))
    return [(i + 2, d) for i, d in enumerate(r)]

f3 = lambda x: f"{float(x):.3f}"
g3 = lambda x: f"{float(x):.3g}"
V = ["# VALUES: multivariable survival for the miRNAs of the survival meta-analysis", "",
     "Lines count the header as line 1. Every table is written by `mirna_multivariable_cox.R` (run log `run.log`).",
     "HR is per SD of the miRNA (z-scored within cohort, within platform for GSE37405).", "",
     "## 1. Control: Table S8 reproduced (`table_s8_reproduction.csv`)", "",
     "The unadjusted primary-endpoint pool equals Table S8 for all 11 miRNAs (HR, CI, p, I² and q within relative",
     "tolerance 1e-9; column `identical_within_1e9`).", "",
     "| line | miRNA | k | HR (RE) | 95 % CI | p | identical |", "|---|---|---|---|---|---|---|"]
for i, r in rows("table_s8_reproduction.csv"):
    V.append(f"| {i} | {r['miRNA']} | {r['k']} | {f3(r['HR'])} | {f3(r['lo'])}–{f3(r['hi'])} | {g3(r['p'])} | {r['identical_within_1e9']} |")

V += ["", "## 2. Covariates available per cohort (`cohort_covariates.csv`; columns n_samples, n_available, levels)", "",
      "| line | cohort | covariate | samples | available | levels |", "|---|---|---|---|---|---|"]
for i, r in rows("cohort_covariates.csv"):
    V.append(f"| {i} | {r['cohort']} | {r['covariate']} | {r['n_samples']} | {r['n_available']} | {r['levels']} |")
V += ["", "METABRIC: NOT POSSIBLE (no miRNA profiles in its public release). No cohort other than TCGA-BRCA has stage.", ""]

P = rows("pooled_adjusted.csv")
V += ["## 3. Pooled adjusted HR, primary endpoint (`pooled_adjusted.csv`)", "",
      "Columns HR_FE, lo_FE, hi_FE, p_FE, q_BH_FE (fixed effect) and HR_RE, lo_RE, hi_RE, p_RE, q_BH_RE (DerSimonian–Laird),",
      "I2, k. BH is over the 11 miRNAs, separately for FE and RE.", "",
      "| line | analysis | miRNA | k | HR FE (95 % CI) | p FE | q FE | HR RE (95 % CI) | p RE | q RE | I² % |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in P:
    V.append(f"| {i} | {r['analysis']} | {r['miRNA']} | {r['k']} | {f3(r['HR_FE'])} ({f3(r['lo_FE'])}–{f3(r['hi_FE'])}) | "
             f"{g3(r['p_FE'])} | {g3(r['q_BH_FE'])} | {f3(r['HR_RE'])} ({f3(r['lo_RE'])}–{f3(r['hi_RE'])}) | {g3(r['p_RE'])} | "
             f"{g3(r['q_BH_RE'])} | {float(r['I2']):.1f} |")
adj = [(i, r) for i, r in P if r["analysis"].startswith("primary endpoint, adjusted")]
sig = [r["miRNA"] for i, r in adj if float(r["p_RE"]) < 0.05 and float(r["q_BH_RE"]) < 0.05]
V += ["", f"Adjusted, random effects, p < 0.05 and q < 0.05: {', '.join(sig) if sig else 'none'} "
      f"(lines {', '.join(str(i) for i, r in adj if r['miRNA'] in sig)}).", ""]

for title, fn in (("## 4. Harmonised endpoint: overall survival only (`pooled_os_harmonised.csv`)", "pooled_os_harmonised.csv"),
                  ("## 5. Adjusted pool without the stromal score (`pooled_adjusted_without_stromal.csv`)",
                   "pooled_adjusted_without_stromal.csv")):
    V += [title, "", "| line | analysis | miRNA | cohorts | HR FE (95 % CI) | p FE | HR RE (95 % CI) | p RE | q RE | I² % |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in rows(fn):
        V.append(f"| {i} | {r['analysis']} | {r['miRNA']} | {r['cohorts']} | {f3(r['HR_FE'])} ({f3(r['lo_FE'])}–{f3(r['hi_FE'])}) | "
                 f"{g3(r['p_FE'])} | {f3(r['HR_RE'])} ({f3(r['lo_RE'])}–{f3(r['hi_RE'])}) | {g3(r['p_RE'])} | {g3(r['q_BH_RE'])} | "
                 f"{float(r['I2']):.1f} |")
    V.append("")
V += ["The OS-only analysis keeps TCGA-BRCA (OS), GSE19783 (OS; the deposited time variable is the disease-free time) and",
      "GSE37405 (OS). It drops GSE22216 (DRFS), GSE59829 (DMFS) and GSE78870 (TTP).", ""]

C = rows("percohort_cox.csv")
V += ["## 6. Per-cohort fits for hsa-miR-29a (`percohort_cox.csv`)", "",
      "Columns n, nevent, HR, lo, hi, p, schoenfeld_global_p, schoenfeld_mirna_p, covariate_set.", "",
      "| line | cohort | endpoint | model | covariates | n | events | HR (95 % CI) | p | Schoenfeld global p | Schoenfeld miRNA p |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in C:
    if r["miRNA"] == "miR-29a" and not (r["endpoint_set"] == "OS" and r["cohort"] in ("TCGA_BRCA", "GSE19783")):
        V.append(f"| {i} | {r['cohort']} | {r['endpoint']} ({r['endpoint_set']}) | {r['model']} | {r['covariate_set']} | {r['n']} | "
                 f"{r['nevent']} | {f3(r['HR'])} ({f3(r['lo'])}–{f3(r['hi'])}) | {g3(r['p'])} | {g3(r['schoenfeld_global_p'])} | "
                 f"{g3(r['schoenfeld_mirna_p'])} |")
V += ["", "For TCGA-BRCA and GSE19783 the primary endpoint is OS, so their OS rows repeat the primary rows.", ""]
fits = [r for i, r in C if not (r["endpoint_set"] == "OS" and r["cohort"] in ("TCGA_BRCA", "GSE19783"))]
nm = sum(1 for r in fits if r["schoenfeld_mirna_p"] not in ("", "NA") and float(r["schoenfeld_mirna_p"]) < 0.05)
ng = sum(1 for r in fits if r["schoenfeld_global_p"] not in ("", "NA") and float(r["schoenfeld_global_p"]) < 0.05)
bad = [f"{r['cohort']} {r['miRNA']} {r['endpoint']} {r['model']} ({g3(r['schoenfeld_mirna_p'])})" for r in fits
       if r["schoenfeld_mirna_p"] not in ("", "NA") and float(r["schoenfeld_mirna_p"]) < 0.05]
V += ["## 7. Proportional hazards (Schoenfeld, `percohort_cox.csv` columns schoenfeld_mirna_p and schoenfeld_global_p)", "",
      f"- Distinct fits (the duplicated OS rows of TCGA-BRCA and GSE19783 counted once): {len(fits)}.",
      f"- miRNA term p < 0.05: {nm}. Global p < 0.05: {ng}.",
      "- The global test includes the clinical covariates, which carry most of the departures.",
      "- miRNA-term departures: " + "; ".join(bad) + ".", ""]
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
