#!/usr/bin/env python3
"""G2: what the degree-matched driver resampling in exact_control_roles_nolegacy.json measures.
Reads analyses/analysed_network_reruns/v3/exact_control_roles_nolegacy.csv and repeats the candidate-pool
rule of exact_control_roles.py (non-hubs with |degree difference| <= max(2, smallest difference)) to
show where the null mean of 0.784 comes from.  Read-only."""
import csv, os, numpy as np
REV = '/path/to/revision'
R = list(csv.DictReader(open(f'{REV}/analyses/analysed_network_reruns/v3/exact_control_roles_nolegacy.csv')))
deg = np.array([float(r['degree_total']) for r in R]); hub = np.array([r['is_FFL_hub'] == 'True' for r in R])
can = np.array([r['role_exact'] != 'redundant' for r in R]); nm = [r['node'] for r in R]
hub_i = np.where(hub)[0]; non_i = np.where(~hub)[0]
L = [f'nodes {len(R)}; can be a driver (critical or intermittent) {can.sum()}; redundant {(~can).sum()}',
     f'FFL hubs {hub.sum()}: can be a driver {can[hub].sum()} ({100*can[hub].mean():.1f}%); non-hubs {(~hub).sum()}: {can[~hub].sum()} ({100*can[~hub].mean():.1f}%)',
     'redundant nodes: ' + ', '.join(f"{r['node']} ({r['type']}, degree {r['degree_total']}, hub {r['is_FFL_hub']})" for r in R if r['role_exact'] == 'redundant')]
exp = []; single = []
rows = []
for i in hub_i:
    d = np.abs(deg[non_i] - deg[i]); c = non_i[d <= max(2.0, d.min())]
    exp.append(can[c].mean())
    rows.append([nm[i], int(deg[i]), len(c), int((~can[c]).sum()), ';'.join(nm[j] for j in c if not can[j])])
    if len(c) == 1: single.append((nm[i], int(deg[i]), nm[c[0]], int(deg[c[0]]), 'redundant' if not can[c[0]] else 'can be driver'))
L.append(f'analytic expectation of the degree-matched null mean = {np.mean(exp):.4f} (simulated 0.7842, 10,000 resamples)')
L.append(f'hubs whose matched pool is a single node: {len(single)} -> ' + '; '.join(f'{a} (deg {b}) -> {c} (deg {d}, {e})' for a, b, c, d, e in single))
out = os.path.dirname(os.path.abspath(__file__))
with open(f'{out}/g2_matched_pools.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['hub', 'degree', 'matched_pool_size', 'redundant_in_pool', 'redundant_nodes_in_pool']); w.writerows(rows)
open(f'{out}/g2_degree_matched_driver_test.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
