#!/usr/bin/env python3
"""(A cont.) Are the architectural statistics anything other than a consequence of the
degree sequence?

Two nulls, matching the established motif-null construction in results/motif_significance_FINAL.md:
  NULL-A  curveball randomisation within each edge class (destroys TF<->miRNA reciprocity)
  NULL-B  as A but the 1,223 reciprocal TF<->miRNA dyads are held fixed
For each randomisation recompute: LSCC size, bow-tie sector sizes, flow hierarchy,
global reaching centrality, trophic incoherence F0, undirected transitivity, FFL cores.
"""
import os, sys, csv, json, collections
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES, ffl_cores, ffl_unique_count

rng = np.random.default_rng(20260908)
NRAND = int(sys.argv[1]) if len(sys.argv) > 1 else 200

nt = load_nodes()
E = load_edges()
cls_edges = collections.defaultdict(list)
for e in E:
    cls_edges[e['edge_type']].append((e['source'], e['target']))
is_tf = {n for n, t in nt.items() if t == 'TF'}


def curveball(edges, n_iter=None):
    adj = collections.defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
    keys = list(adj)
    if len(keys) < 2:
        return list(edges)
    n_iter = n_iter or 5 * len(keys)
    for _ in range(n_iter):
        a, b = rng.choice(len(keys), 2, replace=False)
        A_, B_ = keys[a], keys[b]
        sa, sb = adj[A_], adj[B_]
        common = sa & sb
        ex_a = list(sa - common); ex_b = list(sb - common)
        pool = ex_a + ex_b
        if not pool:
            continue
        rng.shuffle(pool)
        na = len(ex_a)
        adj[A_] = common | set(pool[:na])
        adj[B_] = common | set(pool[na:])
    return [(u, v) for u in adj for v in adj[u]]


tfm = set(cls_edges['TF_miRNA'])
mirtf = {(u, v) for u, v in cls_edges['miRNA_target'] if v in is_tf}
recip = {(t, m) for (t, m) in tfm if (m, t) in mirtf}
recip_rev = {(m, t) for (t, m) in recip}
print('reciprocal TF<->miRNA dyads:', len(recip))


def randomise(keep_recip):
    new = []
    for k, es in cls_edges.items():
        es = set(es)
        if keep_recip and k == 'TF_miRNA':
            new += list(recip) + curveball(list(es - recip))
        elif keep_recip and k == 'miRNA_target':
            new += list(recip_rev) + curveball(list(es - recip_rev))
        else:
            new += curveball(list(es))
    return set(new)


def trophic_F0(G):
    wcc = max(nx.weakly_connected_components(G), key=len)
    Gw = G.subgraph(wcc)
    ns = sorted(Gw.nodes()); idx = {n: i for i, n in enumerate(ns)}
    n = len(ns)
    A = np.zeros((n, n))
    for u, v in Gw.edges():
        A[idx[u], idx[v]] = 1.0
    kin = A.sum(0); kout = A.sum(1)
    Lam = np.diag(kin + kout) - (A + A.T)
    h, *_ = np.linalg.lstsq(Lam, kin - kout, rcond=None)
    num = sum((h[idx[v]] - h[idx[u]] - 1) ** 2 for u, v in Gw.edges())
    return float(np.sqrt(num / Gw.number_of_edges()))


def stats(pairs):
    G = nx.DiGraph(); G.add_nodes_from(nt); G.add_edges_from(pairs)
    sccs = sorted(nx.strongly_connected_components(G), key=len, reverse=True)
    L = sccs[0]
    rep = next(iter(L))
    IN = nx.ancestors(G, rep) - L
    OU = nx.descendants(G, rep) - L
    cores = ffl_cores(set(pairs), nt)
    return dict(LSCC=len(L), n_SCC=len(sccs), IN=len(IN), OUT=len(OU),
                flow_hierarchy=float(nx.flow_hierarchy(G)),
                GRC=float(nx.global_reaching_centrality(G)),
                F0=trophic_F0(G),
                transitivity=float(nx.transitivity(G.to_undirected())),
                reciprocity=float(nx.reciprocity(G)),
                ffl_unique=ffl_unique_count(cores), ffl_raw=len(cores))


obs = stats({(e['source'], e['target']) for e in E})
print('OBSERVED:', obs)
res = {'observed': obs, 'NRAND': NRAND}
rows = []
for null_name, keep in [('NULL_A', False), ('NULL_B', True)]:
    acc = collections.defaultdict(list)
    for i in range(NRAND):
        s = stats(randomise(keep))
        for k, v in s.items():
            acc[k].append(v)
        if (i + 1) % 50 == 0:
            print(null_name, i + 1, flush=True)
    res[null_name] = {}
    for k, v in acc.items():
        v = np.array(v, float)
        o = obs[k]
        z = (o - v.mean()) / v.std(ddof=1) if v.std(ddof=1) > 0 else float('nan')
        p_hi = (np.sum(v >= o) + 1) / (NRAND + 1)
        p_lo = (np.sum(v <= o) + 1) / (NRAND + 1)
        res[null_name][k] = dict(obs=float(o), null_mean=float(v.mean()),
                                 null_sd=float(v.std(ddof=1)), z=float(z),
                                 p_two=float(min(1.0, 2 * min(p_hi, p_lo))))
        rows.append([null_name, k, o, round(float(v.mean()), 4), round(float(v.std(ddof=1)), 4),
                     round(float(z), 3), round(float(min(1.0, 2 * min(p_hi, p_lo))), 4)])

with open(f'{RES}/systems_architecture_nulls.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['null_model', 'metric', 'observed', 'null_mean', 'null_sd', 'z', 'p_two_sided'])
    w.writerows(rows)
with open(f'{RES}/systems_architecture_nulls.json', 'w') as fh:
    json.dump(res, fh, indent=2)
for r in rows:
    print('\t'.join(str(x) for x in r))
