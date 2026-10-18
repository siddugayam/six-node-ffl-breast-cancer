#!/usr/bin/env python3
"""Writes VALUES.md for the dynamics sensitivity analysis from the tables in this folder: every number with its file,
line (header = line 1) and column.
usage (from this folder): python3 write_values.py"""
import csv

def rows(fn):
    with open(fn, newline="") as fh:
        return [(i + 2, d) for i, d in enumerate(csv.DictReader(fh))]

g3 = lambda x: f"{float(x):.3g}"
V = ["# VALUES: dynamics of the higher-order modules (Table 4, Fig. 5)", "",
     "Lines count the header as line 1. Written by `behaviour_intervals_and_prcc.py` (a, b, d), `run_excursions.py` and",
     "`summarise_excursions.py` (c). n = 16,384 Sobol points per module.", "",
     "## (a) Exact 95 % Clopper–Pearson intervals (`behaviour_frequencies_clopper_pearson.csv`)", "",
     "Columns k, n, pct, cp95_lo_pct, cp95_hi_pct. Flag 'uncertified (as stored)' is the per-point flag of the sweep",
     "(Table 4); 'certified limit cycle' rows use the 1,000-tau certification.", "",
     "| line | family | module | behaviour | k | % | 95 % CI (%) | flag |", "|---|---|---|---|---|---|---|---|"]
for i, r in rows("behaviour_frequencies_clopper_pearson.csv"):
    V.append(f"| {i} | {r['family']} | {r['module']} | {r['behaviour']} | {r['k']} | {float(r['pct']):.3f} | "
             f"{float(r['cp95_lo_pct']):.3f}–{float(r['cp95_hi_pct']):.3f} | {r['flag']} |")
V += ["", "## (b) PRCC, top 5 per module and behaviour (`prcc_top5.csv`; all parameters in `prcc_binary.csv`)", "",
      "Columns k (points with the behaviour), rank, parameter, prcc, p, q_bh. p values below about 1e-300 print as 0.", "",
      "| line | family | module | behaviour | k | rank | parameter | PRCC | p | q (BH) |", "|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("prcc_top5.csv"):
    V.append(f"| {i} | {r['family']} | {r['module']} | {r['behaviour']} | {r['k']} | {r['rank']} | {r['parameter']} | "
             f"{float(r['prcc']):+.3f} | {g3(r['p'])} | {g3(r['q_bh'])} |")
ne = [(i, r) for i, r in rows("prcc_binary.csv") if r["note"]]
V += ["", "Not estimable (fewer than 10 points with or without the behaviour; `prcc_binary.csv` column note):", ""]
V += [f"- line {i}: {r['family']} {r['module']} {r['behaviour']}: {r['note']}" for i, r in ne]
V += ["", "Parameters that act in each module: `parameters_acting_by_module.csv` (column n_acting):", ""]
V += [f"- line {i}: {r['family']} {r['module']}: {r['n_acting']} acting; not acting: {r['not_acting']}"
      for i, r in rows("parameters_acting_by_module.csv")]
V += ["", "The continuous-proxy PRCCs (secondary) are in `prcc_continuous.csv` (columns proxy, n_points, parameter, prcc, p,",
      "q_bh, rank_abs).", "",
      "## (c) One-at-a-time excursions (`excursion_frequencies.csv`)", "",
      "Columns baseline_k, baseline_pct, k, pct, cp95_lo_pct, cp95_hi_pct, gained, lost, mcnemar_p, status. Rows for memory,",
      "pulse and sustained oscillation are listed here; every behaviour is in the file.", "",
      "| line | family | module | parameter | factor | behaviour | baseline % | excursion % (95 % CI) | gained / lost | McNemar p | status |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("excursion_frequencies.csv"):
    if r["behaviour"] not in ("is_memory", "is_pulse", "is_sustained_osc"): continue
    if r["status"] != "run":
        V.append(f"| {i} | {r['family']} | {r['module']} | {r['parameter']} | x{r['factor']} | {r['behaviour']} | "
                 f"{float(r['baseline_pct']):.2f} | | | | NOT RUN |"); continue
    V.append(f"| {i} | {r['family']} | {r['module']} | {r['parameter']} | x{r['factor']} | {r['behaviour']} | "
             f"{float(r['baseline_pct']):.2f} | {float(r['pct']):.2f} ({float(r['cp95_lo_pct']):.2f}–{float(r['cp95_hi_pct']):.2f}) | "
             f"{r['gained']} / {r['lost']} | {g3(r['mcnemar_p'])} | run |")
V += ["", "### Range over the excursions, per module and behaviour (`excursion_range_by_module.csv`, descriptive)", "",
      "Columns baseline_pct, min_pct (min_at), max_pct (max_at), median_abs_change_pp, max_abs_change_pp, appears (0 at",
      "baseline, above 0 in an excursion), vanishes (above 0 at baseline, 0 in an excursion).", "",
      "| line | family | module | behaviour | excursions | baseline % | lowest % (excursion) | highest % (excursion) | median abs change (pp) | max abs change (pp) | appears | vanishes |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in rows("excursion_range_by_module.csv"):
    V.append(f"| {i} | {r['family']} | {r['module']} | {r['behaviour']} | {r['excursions']} | {float(r['baseline_pct']):.2f} | "
             f"{float(r['min_pct']):.2f} ({r['min_at']}) | {float(r['max_pct']):.2f} ({r['max_at']}) | {float(r['median_abs_change_pp']):.2f} | "
             f"{float(r['max_abs_change_pp']):.2f} | {r['appears']} | {r['vanishes']} |")
V += ["", "## (d) Parameter table (`parameter_table.csv`)", "",
      "| line | parameter | meaning | low | high | distribution | acts in |", "|---|---|---|---|---|---|---|"]
for i, r in rows("parameter_table.csv"):
    V.append(f"| {i} | {r['parameter']} | {r['meaning']} | {r['low']} | {r['high']} | {r['distribution']} | {r['acts_in_modules']} |")
V += ["", "Sampler (every row): " + rows("parameter_table.csv")[0][1]["sampler"] + ".", ""]
open("VALUES.md", "w").write("\n".join(V) + "\n")
print(len(V), "lines")
