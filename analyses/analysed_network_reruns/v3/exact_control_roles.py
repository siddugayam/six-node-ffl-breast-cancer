#!/usr/bin/env python3
"""Exact (matching-invariant) structural-controllability roles, replacing the single-configuration
driver flag behind Note S5's hub comparison.  That flag depends on which of many maximum matchings
the solver returns, so the stored run, a re-run and Note S5 all disagree (see VALUES.md).

For the bipartite out-copy/in-copy graph with a maximum matching M, the in-copy of node v is
unmatched in SOME maximum matching iff it is reachable from an M-exposed in-copy by an even-length
M-alternating path (Dulmage-Mendelsohn).  It is unmatched in EVERY maximum matching iff v has
in-degree 0.  Hence, exactly and independently of M:
  critical      driver in every maximum matching   (in-degree 0)
  intermittent  driver in some, not all
  redundant     driver in none
The hub test of scripts/v3/12_controllability_degree_control.py is repeated on 'can be a driver'
(critical or intermittent): logistic regression on log(degree+1) and node type, and degree-matched
resampling (same design and seed).  Deletion classes are taken from 03_controllability.py's output
(they are deterministic).
usage: python3 exact_control_roles.py <full|nolegacy>   (reads out_<mode>_hash0/ from run_v3.py)"""
import sys, os, csv, json, collections
import numpy as np
import networkx as nx
from networkx.algorithms.bipartite import hopcroft_karp_matching
import statsmodels.api as sm
sys.dont_write_bytecode = True
V3 = '/path/to/revision/scripts/v3'
HERE = os.path.dirname(os.path.abspath(__file__))
mode = sys.argv[1]; OUT = os.path.join(HERE, f'out_{mode}_hash0')
sys.path.insert(0, V3)
import netlib
nt = netlib.load_nodes(); E = netlib.load_edges()
if mode == 'nolegacy':
    E = [e for e in E if e['edge_type'] != 'miRNA_miRNA']
arcs = sorted({(e['source'], e['target']) for e in E if e['source'] != e['target']})
NODES = sorted(nt); N = len(NODES)
B = nx.Graph()
L = [('o', n) for n in NODES]; R = [('i', n) for n in NODES]
B.add_nodes_from(L, bipartite=0); B.add_nodes_from(R, bipartite=1)
B.add_edges_from((('o', u), ('i', v)) for u, v in arcs)
M = hopcroft_karp_matching(B, top_nodes=L)
size = sum(1 for k in M if k[0] == 'o')
ND = max(N - size, 1)
# even alternating reachability from exposed in-copies
exposed = [r for r in R if r not in M]
seen = set(exposed); stack = list(exposed)
while stack:
    r = stack.pop()
    for l in B[r]:                                   # non-matching edge right -> left
        if M.get(l) == r: continue
        r2 = M.get(l)                                # matching edge left -> right
        if r2 is not None and r2 not in seen:
            seen.add(r2); stack.append(r2)
indeg = collections.Counter(v for _, v in arcs); outdeg = collections.Counter(u for u, _ in arcs)
role = {}
for n in NODES:
    if indeg[n] == 0: role[n] = 'critical'
    elif ('i', n) in seen: role[n] = 'intermittent'
    else: role[n] = 'redundant'
ctrl = {r['node']: r for r in csv.DictReader(open(f'{OUT}/systems_controllability_nodes.csv'))}
# consistency with 03's sampled frequencies: a node ever sampled as a driver must be critical/intermittent
bad = [n for n in NODES if float(ctrl[n]['driver_frequency']) > 0 and role[n] == 'redundant']
assert not bad, bad
assert int(ctrl[NODES[0]]['ND_after_deletion']) >= 0
can = np.array([role[n] != 'redundant' for n in NODES])
deg = np.array([indeg[n] + outdeg[n] for n in NODES], float)
hub = np.array([ctrl[n]['is_FFL_hub'] == 'True' for n in NODES])
typ = [nt[n] for n in NODES]
X = sm.add_constant(np.column_stack([np.log1p(deg), [1.0 if t == 'TF' else 0.0 for t in typ],
                                     [1.0 if t == 'Gene' else 0.0 for t in typ], hub.astype(float)]))
fit = sm.Logit(can.astype(float), X).fit(disp=0)
rng = np.random.default_rng(11)
hub_i = np.where(hub)[0]; non_i = np.where(~hub)[0]; obs = can[hub_i].mean(); NB = 10000
nullr = np.empty(NB)
for b in range(NB):
    pick = []
    for i in hub_i:
        d = np.abs(deg[non_i] - deg[i]); cand = non_i[d <= max(2.0, d.min())]
        pick.append(rng.choice(cand))
    nullr[b] = can[pick].mean()
p_low = (np.sum(nullr <= obs) + 1) / (NB + 1); p_high = (np.sum(nullr >= obs) + 1) / (NB + 1)
module = ['hsa-miR-29a', 'hsa-miR-29b', 'hsa-miR-29c', 'COL1A1', 'COL3A1']
summ = dict(mode=mode, N=N, arcs=len(arcs), matching_size=size, N_D=ND, n_D_over_N=ND / N,
            roles_exact=dict(collections.Counter(role.values())),
            deletion_classes=dict(collections.Counter(ctrl[n]['deletion_class'] for n in NODES)),
            FFL_hubs=int(hub.sum()), can_be_driver_hubs=float(obs), can_be_driver_nonhubs=float(can[non_i].mean()),
            logistic_is_FFL_hub=dict(coef=float(fit.params[4]), se=float(fit.bse[4]), p=float(fit.pvalues[4])),
            degree_matched=dict(null_mean=float(nullr.mean()), null_sd=float(nullr.std(ddof=1)),
                                p_hubs_lower=float(p_low), p_hubs_higher=float(p_high), n_resamples=NB),
            mir29_collagen_module={n: dict(role=role[n], deletion_class=ctrl[n]['deletion_class'],
                                           N_D_after_deletion=int(ctrl[n]['ND_after_deletion'])) for n in module})
with open(os.path.join(HERE, f'exact_control_roles_{mode}.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['node', 'type', 'degree_total', 'is_FFL_hub', 'role_exact', 'deletion_class', 'driver_frequency_sampled_03'])
    for n in NODES: w.writerow([n, nt[n], int(deg[NODES.index(n)]), ctrl[n]['is_FFL_hub'], role[n], ctrl[n]['deletion_class'], ctrl[n]['driver_frequency']])
json.dump(summ, open(os.path.join(HERE, f'exact_control_roles_{mode}.json'), 'w'), indent=1)
print(json.dumps(summ, indent=1))
