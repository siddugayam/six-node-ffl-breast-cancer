#!/usr/bin/env python3
"""Rebuild Table S6 (direction-aware hub table), for all 587 nodes ordered by undirected current-flow betweenness,
from the Note S5 information-flow and controllability outputs (../v3/out_<mode>_hash0/).
No script in the project writes results/v5/tables/TableS6_directed_hub_table.csv, so the rank rules are
inferred and checked: for the full network the rebuild must reproduce the stored table exactly.
Ranks: descending, ties share the lowest rank (the tied out-degree-zero nodes all get rank 377).
The single-matching 'is_driver' column is replaced by the exact matching-invariant role
(../v3/exact_control_roles_<mode>.csv); the stored column is reproduced only for the check."""
import csv, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.join(HERE, '..', 'v3')
def rank_min(vals):
    order = sorted(vals, key=lambda k: -vals[k]); r = {}; prev = None
    for i, k in enumerate(order, 1):
        if prev is None or vals[k] != vals[prev]: cur = i
        r[k] = cur; prev = k
    return r
def build(mode):
    F = {r['node']: r for r in csv.DictReader(open(f'{V3}/out_{mode}_hash0/systems_information_flow.csv'))}
    C = {r['node']: r for r in csv.DictReader(open(f'{V3}/out_{mode}_hash0/systems_controllability_nodes.csv'))}
    X = {r['node']: r for r in csv.DictReader(open(f'{V3}/exact_control_roles_{mode}.csv'))}
    f = lambda c: {n: float(F[n][c]) for n in F}
    rd, rc, ri, rw = rank_min(f('degree_total')), rank_min(f('current_flow_betweenness')), rank_min(f('information_centrality')), rank_min(f('rw_throughput_betweenness'))
    top = sorted(F, key=lambda n: (rc[n], n))          # all 587 nodes, as the Table S6 caption states
    rows = []
    for n in top:
        rows.append(dict(node=n, type=F[n]['type'], degree_total=int(float(F[n]['degree_total'])), ffl_participation=int(float(F[n]['ffl_participation'])),
                         rank_degree=rd[n], rank_cfb=rc[n], rank_infocent=ri[n], rank_directed=rw[n], rank_shift=rw[n] - rc[n],
                         is_sink=str(int(float(F[n]['degree_out'])) == 0).upper(), is_driver=C[n]['is_driver_one_config'].upper(),
                         control_role_exact=X[n]['role_exact'], deletion_class=C[n]['deletion_class']))
    return rows
cols_stored = ['node','type','degree_total','ffl_participation','rank_degree','rank_cfb','rank_infocent','rank_directed','rank_shift','is_sink','is_driver','deletion_class']
full = build('full')[:40]                        # the stored table holds the top 40 only
stored = list(csv.DictReader(open('/path/to/revision/results/v5/tables/TableS6_directed_hub_table.csv')))
diff = [(i, c, str(a[c]), b[c]) for i, (a, b) in enumerate(zip(full, stored)) for c in cols_stored if str(a[c]) != b[c]]
print('rows', len(full), len(stored), ' differing cells vs stored (full network):', len(diff))
for d in diff[:20]: print('  ', d)
out = build('nolegacy')
cols = ['node','type','degree_total','ffl_participation','rank_degree','rank_cfb','rank_infocent','rank_directed','rank_shift','is_sink','control_role_exact','deletion_class']
with open(os.path.join(HERE, 'TableS6_directed_hub_table_nolegacy.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore'); w.writeheader(); w.writerows(out)
both = [r['node'] for r in out if r['rank_cfb'] <= 10 and r['rank_directed'] <= 10]
out40 = out[:40]
both_full = [r['node'] for r in full if r['rank_cfb'] <= 10 and r['rank_directed'] <= 10]
print('top-10 on both undirected current-flow and directed throughput: full', both_full, ' nolegacy', both)
print('rows written (nolegacy):', len(out), '; sinks among top 40: full', sum(r['is_sink'] == 'TRUE' for r in full), ' nolegacy', sum(r['is_sink'] == 'TRUE' for r in out40))
