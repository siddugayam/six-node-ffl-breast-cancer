#!/usr/bin/env python3
"""Writes VALUES.md for the oestrogen-receptor and PAM50 adjustment from the tables in this folder: every number with its
file, line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

f3 = lambda x: f"{float(x):.3f}"
g3 = lambda x: f"{float(x):.3g}"
V = ["# VALUES: oestrogen-receptor and PAM50 adjustment", "",
     "Lines count the header as line 1. Part 1 tables are written by `circuit_survival_er_pam50.R` (log `run_circuits.log`),",
     "part 2 tables by `mir29a_er_pam50.R` (log `run_mir29a.log`).", "",
     "## Part 1. Six-node circuits, TCGA-BRCA, overall survival", "",
     "### Control (`control_reproduction.csv`)", ""]
for i, r in rows("control_reproduction.csv"):
    V.append(f"- line {i}: {r['rows_compared']} of {r['rows_published']} published rows compared; n and events equal "
             f"{r['n_equal']}; largest relative differences HR {float(r['max_rel_HR']):.1e}, p {float(r['max_rel_p_six']):.1e}, "
             f"LRT p {float(r['max_rel_LRT']):.1e}, q {float(r['max_rel_q']):.1e}; reproduced {r['reproduced']}.")
V += ["", "### Summary per fit and family (`circuit_summary.csv`)", "",
      "Fits: published = age + stage on the published sample; same_tumours = age + stage on the tumours with ER and PAM50;",
      "er_pam50 = age + stage + ER + PAM50 on those tumours. LRT = core vs core + extension; BH within fit and family.", "",
      "| line | fit | family | circuits | n | events | LRT p < 0.05 | q < 0.05 | HR p < 0.05 | the 16: LRT p < 0.05 | the 16: q < 0.05 | the 16: HR p < 0.05 | median ΔC |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("circuit_summary.csv"):
    V.append(f"| {i} | {r['fit']} | {r['family']} | {r['circuits']} | {r['n']} | {r['events']} | {r['n_LRT_p05']} | {r['n_q_lt_005']} | "
             f"{r['n_HR_p05']} | {r['sixteen_LRT_p05']} | {r['sixteen_q_lt_005']} | {r['sixteen_HR_p05']} | {float(r['median_deltaC']):+.5f} |")
V += ["", "In the six-node model architecture family, the columns for the 16 count none, because the 16 are all six-node",
      "composite pattern circuits. Median ΔC is the median over circuits of C_six − C_core (column median_deltaC).", "",
      "### Circuits with q < 0.05 in the same-tumour or ER + PAM50 fits (`circuits_er_pam50.csv`)", "",
      "| line | fit | family | members | core | extension | HR (95 % CI) | p | LRT p | q |", "|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("circuits_er_pam50.csv"):
    if r["fit"] != "published" and float(r["q_BH"]) < 0.05:
        V.append(f"| {i} | {r['fit']} | {r['family']} | {r['members']} | {r['core']} | {r['extension']} | {f3(r['HR_six'])} "
                 f"({f3(r['HR_lo'])}–{f3(r['HR_hi'])}) | {g3(r['p_six'])} | {g3(r['LRT_p_extension_adds'])} | {g3(r['q_BH'])} |")
six = {r["members"] for _, r in rows("sixteen_patterns_before_after.csv")}
hits = [r["members"] for _, r in rows("circuits_er_pam50.csv") if r["fit"] != "published" and float(r["q_BH"]) < 0.05]
V += ["", f"Of these {len(hits)}, {sum(m in six for m in hits)} are among the 16 patterns.", "",
      "### The 16 patterns before and after (`sixteen_patterns_before_after.csv`)", "",
      "HR of the six-node score per SD (95 % CI), its p, and the LRT p (q) of each fit.", "",
      "| line | members | published: HR (CI), p; LRT p (q) | same tumours: HR (CI), p; LRT p (q) | ER + PAM50: HR (CI), p; LRT p (q) |",
      "|---|---|---|---|---|"]
for i, r in rows("sixteen_patterns_before_after.csv"):
    cell = lambda f: (f"{f3(r['HR_six_' + f])} ({f3(r['HR_lo_' + f])}–{f3(r['HR_hi_' + f])}), {g3(r['p_six_' + f])}; "
                      f"{g3(r['LRT_p_extension_adds_' + f])} ({g3(r['q_BH_' + f])})")
    V.append(f"| {i} | {r['members']} | {cell('published')} | {cell('same_tumours')} | {cell('er_pam50')} |")
V += ["", "Every circuit and fit: `circuits_er_pam50.csv` (columns fit, covariates, family, members, core, extension, n, events, HR_six, HR_lo, HR_hi,",
      "p_six, C_core, C_six, LRT_p_extension_adds, q_BH). METABRIC: NOT POSSIBLE (no miRNA assay; the circuits cannot be scored).", "",
      "## Part 2. hsa-miR-29a and the other miRNAs, adjusted for ER and/or PAM50", "",
      "### Pooled (`mir29a_er_pam50_pooled.csv`)", "",
      "| line | analysis | miRNA | k | n | events | HR FE (95 % CI) | p FE | HR RE (95 % CI) | p RE | q RE | I² % |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("mir29a_er_pam50_pooled.csv"):
    V.append(f"| {i} | {r['analysis']} | {r['miRNA']} | {r['k']} | {r['n']} | {r['nevent']} | {f3(r['HR_FE'])} ({f3(r['lo_FE'])}–{f3(r['hi_FE'])}) | "
             f"{g3(r['p_FE'])} | {f3(r['HR_RE'])} ({f3(r['lo_RE'])}–{f3(r['hi_RE'])}) | {g3(r['p_RE'])} | {g3(r['q_BH_RE'])} | {float(r['I2']):.1f} |")
V += ["", "### Per cohort, hsa-miR-29a (`mir29a_er_pam50_percohort.csv`)", "",
      "| line | cohort | endpoint | model | covariates | n | events | HR (95 % CI) | p | Schoenfeld miRNA p |", "|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("mir29a_er_pam50_percohort.csv"):
    if r["miRNA"] == "miR-29a":
        V.append(f"| {i} | {r['cohort']} | {r['endpoint']} | {r['model']} | {r['covariates']} | {r['n']} | {r['nevent']} | "
                 f"{f3(r['HR'])} ({f3(r['lo'])}–{f3(r['hi'])}) | {g3(r['p'])} | {g3(r['schoenfeld_mirna_p'])} |")
V += ["", "GSE37405 is ER-positive by design and has no subtype annotation; it is not in this pool.", ""]
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
