#!/usr/bin/env python3
"""(B) STRUCTURAL CONTROLLABILITY of the canonical network.

Liu, Slotine & Barabasi (2011) Nature 473:167-173. For a linear time-invariant system
dx/dt = Ax + Bu on a directed network, the minimum number of driver nodes needed for full
structural controllability is
        N_D = max(N - |M*|, 1)
where M* is a maximum matching of the bipartite representation (out-copies -> in-copies).
The driver nodes of a given control configuration are exactly the nodes left UNMATCHED on
the in-side.

Two further classifications are computed:
  * Liu 2011 / Jia & Barabasi (2013) Nat Commun 4:2002 node roles - critical / intermittent /
    redundant - from the frequency with which a node is a driver across many DISTINCT
    maximum matchings (sampled by randomising the Hopcroft-Karp vertex order).
  * Vinayagam et al. (2016) PNAS 113:4976 - indispensable / neutral / dispensable, by whether
    deleting the node increases / leaves unchanged / decreases N_D.

The biological question: are the FFL hubs and the miR-29/collagen module drivers or
indispensable nodes?
"""
import os, sys, csv, json, collections, random
import numpy as np
import networkx as nx
from networkx.algorithms.bipartite import hopcroft_karp_matching
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES, FIG, ffl_cores

rng = random.Random(20260908)
nt = load_nodes()
E = load_edges()
P = [(e['source'], e['target']) for e in E]
NODES = sorted(nt)


def max_matching_size(pairs, nodes, shuffle=False):
    """Bipartite: '+'+u (out-copy) -> '-'+v (in-copy).  Returns (size, matched_in_set)."""
    B = nx.Graph()
    ns = list(nodes)
    if shuffle:
        ns = ns[:]; rng.shuffle(ns)
    B.add_nodes_from(['+' + n for n in ns], bipartite=0)
    B.add_nodes_from(['-' + n for n in ns], bipartite=1)
    pl = list(pairs)
    if shuffle:
        rng.shuffle(pl)
    B.add_edges_from([('+' + u, '-' + v) for u, v in pl])
    top = ['+' + n for n in ns]
    m = hopcroft_karp_matching(B, top_nodes=top)
    matched_in = {k[1:] for k in m if k.startswith('-')}
    return len(matched_in), matched_in


N = len(NODES)
size, matched = max_matching_size(P, NODES)
ND = max(N - size, 1)
drivers = set(NODES) - matched
print(f'N={N}  |M*|={size}  N_D={ND}  n_D/N={ND/N:.4f}  drivers listed={len(drivers)}')
assert len(drivers) == ND or ND == 1

# --------------------------------------------------- sample distinct maximum matchings
NSAMP = int(os.environ.get('NSAMP', '2000'))
freq = collections.Counter()
sizes = set()
for i in range(NSAMP):
    s, mi = max_matching_size(P, NODES, shuffle=True)
    sizes.add(s)
    for n in set(NODES) - mi:
        freq[n] += 1
assert len(sizes) == 1, f'matching size not invariant: {sizes}'
role = {}
for n in NODES:
    f = freq[n] / NSAMP
    role[n] = ('critical' if f > 0.999 else 'redundant' if f < 0.001 else 'intermittent')
print('roles:', collections.Counter(role.values()))

# --------------------------------------------------- node deletion -> indispensable etc.
adj_by_node = collections.defaultdict(list)
for u, v in P:
    adj_by_node[u].append((u, v)); adj_by_node[v].append((u, v))
cls = {}; ND_del = {}
for n in NODES:
    rest = [n2 for n2 in NODES if n2 != n]
    pr = [(u, v) for u, v in P if u != n and v != n]
    s2, _ = max_matching_size(pr, rest)
    nd2 = max(len(rest) - s2, 1)
    ND_del[n] = nd2
    cls[n] = ('indispensable' if nd2 > ND else 'dispensable' if nd2 < ND else 'neutral')
print('deletion classes:', collections.Counter(cls.values()))

# --------------------------------------------------- context tables
deg = collections.Counter(); ind = collections.Counter(); outd = collections.Counter()
for u, v in P:
    deg[u] += 1; deg[v] += 1; outd[u] += 1; ind[v] += 1
cores = ffl_cores(set(P), nt)
part = collections.Counter()
for R, M, T, c in cores:
    part[R] += 1; part[M] += 1; part[T] += 1

hubs = {}
with open('/path/to/revision/results/ffl_network_hubs.csv') as fh:
    for r in csv.DictReader(fh):
        hubs[r['node']] = r
MIR29 = ['hsa-miR-29a', 'hsa-miR-29b', 'hsa-miR-29c']
MODULE = MIR29 + ['COL1A1', 'COL3A1']

with open(f'{RES}/systems_controllability_nodes.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['node', 'type', 'degree_total', 'degree_in', 'degree_out', 'ffl_participation',
                'is_driver_one_config', 'driver_frequency', 'control_role',
                'ND_after_deletion', 'deletion_class', 'is_FFL_hub', 'is_mir29_module'])
    for n in NODES:
        w.writerow([n, nt[n], deg[n], ind[n], outd[n], part[n],
                    n in drivers, round(freq[n] / NSAMP, 4), role[n],
                    ND_del[n], cls[n],
                    hubs.get(n, {}).get('hub_ffl', 'FALSE') == 'TRUE',
                    n in MODULE])
print('wrote systems_controllability_nodes.csv')

# --------------------------------------------------- enrichment tests
def fisher(setA, prop_fn, labelA, labelP):
    a = sum(1 for n in NODES if n in setA and prop_fn(n))
    b = sum(1 for n in NODES if n in setA and not prop_fn(n))
    c = sum(1 for n in NODES if n not in setA and prop_fn(n))
    d = sum(1 for n in NODES if n not in setA and not prop_fn(n))
    orr, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
    return dict(group=labelA, property=labelP, n_group=a + b, k_group=a,
                frac_group=a / (a + b) if a + b else float('nan'),
                n_rest=c + d, k_rest=c, frac_rest=c / (c + d) if c + d else float('nan'),
                odds_ratio=float(orr), p=float(p))

ffl_hubs = {n for n, r in hubs.items() if r['hub_ffl'] == 'TRUE'}
deg_hubs = {n for n, r in hubs.items() if r['hub_degree'] == 'TRUE'}
top_ffl = set(sorted(NODES, key=lambda n: -part[n])[:20])
tests = []
for gname, gset in [('FFL_hub (hub_ffl==TRUE)', ffl_hubs), ('degree_hub', deg_hubs),
                    ('top20_FFL_participation', top_ffl),
                    ('miR-29 family', set(MIR29)), ('miR-29/collagen module', set(MODULE)),
                    ('TF', {n for n in NODES if nt[n] == 'TF'}),
                    ('miRNA', {n for n in NODES if nt[n] == 'miRNA'}),
                    ('Gene', {n for n in NODES if nt[n] == 'Gene'})]:
    tests.append(fisher(gset, lambda n: n in drivers, gname, 'driver_node'))
    tests.append(fisher(gset, lambda n: cls[n] == 'indispensable', gname, 'indispensable'))
    tests.append(fisher(gset, lambda n: cls[n] == 'dispensable', gname, 'dispensable'))
    tests.append(fisher(gset, lambda n: role[n] == 'critical', gname, 'critical_driver'))
with open(f'{RES}/systems_controllability_enrichment.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(tests[0].keys())); w.writeheader(); w.writerows(tests)
for t in tests:
    if t['group'] in ('FFL_hub (hub_ffl==TRUE)', 'miR-29/collagen module', 'miR-29 family',
                      'top20_FFL_participation'):
        print(f"{t['group']:<26} {t['property']:<16} {t['k_group']}/{t['n_group']} "
              f"({t['frac_group']:.3f}) vs {t['frac_rest']:.3f}  OR={t['odds_ratio']:.3g} p={t['p']:.3g}")

# per-node report for module + top hubs
print('\n--- module / hub nodes ---')
rows = []
for n in MODULE + sorted(ffl_hubs, key=lambda x: -part[x])[:15]:
    if n not in nt:
        continue
    rows.append([n, nt[n], deg[n], part[n], n in drivers, round(freq[n] / NSAMP, 3),
                 role[n], cls[n], ND_del[n] - ND])
    print(f'{n:<14} {nt[n]:<6} deg={deg[n]:<4} ffl={part[n]:<4} driver={n in drivers!s:<5} '
          f'freq={freq[n]/NSAMP:.3f} role={role[n]:<12} {cls[n]:<14} dND={ND_del[n]-ND:+d}')
with open(f'{RES}/systems_controllability_focus_nodes.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['node', 'type', 'degree', 'ffl_participation', 'is_driver', 'driver_frequency',
                'control_role', 'deletion_class', 'delta_ND_on_deletion'])
    w.writerows(rows)

# driver frequency vs FFL participation / degree
dfreq = np.array([freq[n] / NSAMP for n in NODES])
pv = np.array([part[n] for n in NODES], float)
dv = np.array([deg[n] for n in NODES], float)
from scipy.stats import spearmanr
r1 = spearmanr(dfreq, pv); r2 = spearmanr(dfreq, dv)
summary = dict(N=N, matching_size=size, N_D=ND, n_D_over_N=ND / N,
               n_drivers_listed=len(drivers),
               roles=dict(collections.Counter(role.values())),
               deletion_classes=dict(collections.Counter(cls.values())),
               spearman_driverfreq_vs_fflparticipation=[float(r1.statistic), float(r1.pvalue)],
               spearman_driverfreq_vs_degree=[float(r2.statistic), float(r2.pvalue)],
               n_matchings_sampled=NSAMP)
with open(f'{RES}/systems_controllability_summary.json', 'w') as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))
