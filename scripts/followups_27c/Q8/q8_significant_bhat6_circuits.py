#!/usr/bin/env python3
"""INBOX_2026-09-27c, Q8: the BHAT6 composite circuits with q < 0.05 for overall survival in the S8 primary analysis
(GSVA score, Cox adjusted for age and stage; LRT of the six-node extension over its three-node core, BH over the
2,383 circuits).  Read from INBOX_2026-09-27b/S8/s8_circuits.csv (one row per circuit; cohort-level n and events only,
no patient-level data).  Writes q8_significant_bhat6_os_circuits.csv and q8_node_sharing.txt."""
import csv, os, collections
H = os.path.dirname(os.path.abspath(__file__))
SRC = '/path/to/revision/INBOX_2026-09-27b/S8/s8_circuits.csv'
R = [r for r in csv.DictReader(open(SRC)) if r['method'] == 'GSVA_adjusted_age_stage' and r['endpoint'] == 'OS' and r['family'] == 'BHAT6_composite']
S = sorted((r for r in R if float(r['q_BH']) < 0.05), key=lambda r: float(r['LRT_p_extension_adds']))
cols = ['members', 'core', 'extension', 'n', 'events', 'HR_six', 'p_six', 'C_core', 'C_six', 'deltaC', 'LRT_p_extension_adds', 'q_BH']
with open(f'{H}/q8_significant_bhat6_os_circuits.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['rank', 'source_line_in_s8_circuits'] + cols)
    lines = {id(r): i for i, r in enumerate(csv.DictReader(open(SRC)), start=2)}
    allrows = list(csv.DictReader(open(SRC)))
    for k, r in enumerate(S, 1):
        ln = next(i for i, x in enumerate(allrows, start=2) if x == r)
        w.writerow([k, ln] + [r[c] for c in cols])
out = []
P = lambda *a: (out.append(' '.join(str(x) for x in a)), print(*a))
P(f'Q8 BHAT6 composite circuits, OS, GSVA adjusted for age and stage: {len(R)} tested; q < 0.05: {len(S)}; deltaC range '
  f'{min(float(r["deltaC"]) for r in S):+.4f} to {max(float(r["deltaC"]) for r in S):+.4f}; q range {min(float(r["q_BH"]) for r in S):.3g}-{max(float(r["q_BH"]) for r in S):.3g}')
node = collections.Counter(n for r in S for n in r['members'].split(';'))
P('Q8 nodes by the number of the ' + str(len(S)) + ' circuits containing them: ' + '; '.join(f'{n} {c}' for n, c in node.most_common()))
cores = collections.Counter(r['core'] for r in S); ext = collections.Counter(r['extension'] for r in S)
P(f'Q8 distinct three-node cores {len(cores)}: ' + '; '.join(f'{k} ({v})' for k, v in cores.most_common()))
P(f'Q8 distinct extensions {len(ext)}: ' + '; '.join(f'{k} ({v})' for k, v in ext.most_common()))
P(f'Q8 nodes in more than half of the circuits: {[n for n, c in node.items() if c > len(S) / 2]}; in every circuit: {[n for n, c in node.items() if c == len(S)]}')
P(f'Q8 cohort: n {sorted({r["n"] for r in S})}, events {sorted({r["events"] for r in S})}')
core4 = {'STAT3', 'HIF1A', 'F3', 'SERPINE1'}
P(f'Q8 circuits containing all of STAT3, HIF1A, F3 and SERPINE1: {sum(1 for r in S if core4 <= set(r["members"].split(";")))}; containing a miR-17~92 or '
  f'miR-106a~363 member (miR-17, -18a, -18b, -19a, -19b, -20a, -20b, -92a, -106a, -363): {sum(1 for r in S if set(r["members"].split(";")) & {"hsa-miR-17", "hsa-miR-18a", "hsa-miR-18b", "hsa-miR-19a", "hsa-miR-19b", "hsa-miR-20a", "hsa-miR-20b", "hsa-miR-92a", "hsa-miR-106a", "hsa-miR-363"})}; '
  f'deltaC < 0: {sum(1 for r in S if float(r["deltaC"]) < 0)}')
none4 = [(k, r) for k, r in enumerate(S, 1) if not core4 & set(r['members'].split(';'))]
P(f'Q8 circuits with HIF1A or STAT3: {sum(1 for r in S if {"HIF1A", "STAT3"} & set(r["members"].split(";")))}; with none of STAT3, HIF1A, F3, SERPINE1: {len(none4)}: ' + '; '.join(f'rank {k} {r["members"]}' for k, r in none4))
open(f'{H}/q8_node_sharing.txt', 'w').write('\n'.join(out) + '\n')
