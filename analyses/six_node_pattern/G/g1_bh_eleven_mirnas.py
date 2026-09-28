#!/usr/bin/env python3
"""G1: BH q over the eleven miRNAs of the survival meta-analysis (the ten prioritised miRNAs plus
hsa-miR-29a), after removing the miR-130a rows of results/v5/tables/TableS7_mirna_survival_meta.csv.
Input: the RE_pooled row (DerSimonian-Laird, primary endpoint per cohort) of each miRNA, column p.
The stored q_BH column of that table was computed over twelve miRNAs (miR-130a included)."""
import csv, os
REV = '/path/to/revision'
R = [r for r in csv.DictReader(open(f'{REV}/results/v5/tables/TableS7_mirna_survival_meta.csv')) if r['row_type'] == 'RE_pooled']
assert len(R) == 12
R11 = [r for r in R if r['miRNA'] != 'miR-130a']; m = len(R11); assert m == 11
o = sorted(R11, key=lambda r: float(r['p']))
q = [min(1.0, float(r['p']) * m / (i + 1)) for i, r in enumerate(o)]
for i in range(m - 2, -1, -1): q[i] = min(q[i], q[i + 1])            # step-up monotonicity
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'g1_bh_q_eleven_mirnas.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['miRNA', 'pooled_HR', 'pooled_p', 'rank', 'q_BH_11_miRNAs', 'q_BH_stored_12_miRNAs'])
    for i, r in enumerate(o):
        w.writerow([r['miRNA'], r['HR'], r['p'], i + 1, repr(q[i]), r['q_BH']])
        print(f"{r['miRNA']:9s} p={float(r['p']):.4g}  q(11)={q[i]:.4g}  stored q(12)={float(r['q_BH']):.4g}")
