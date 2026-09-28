#!/usr/bin/env python3
"""(B cont.) Independent cross-check of the maximum-matching / driver-node computation.

The main analysis (03_controllability.py) uses networkx's Hopcroft-Karp implementation.
Here the same maximum matching is recomputed with scipy.sparse.csgraph.maximum_bipartite_matching
(a different algorithm and a different code base), and the Liu-Slotine-Barabasi driver count
N_D = max(N - |M*|, 1) is re-derived.  Any disagreement would invalidate section B.

Liu YY, Slotine JJ, Barabasi AL (2011) Nature 473:167-173.
"""
import os, sys, json, csv
import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from networkx.algorithms.bipartite import hopcroft_karp_matching

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES

nt = load_nodes()
E = load_edges()
P = [(e['source'], e['target']) for e in E]
NODES = sorted(nt)
idx = {n: i for i, n in enumerate(NODES)}
N = len(NODES)

# ---- scipy: rows = out-copies (sources), cols = in-copies (targets)
r = np.array([idx[u] for u, v in P])
c = np.array([idx[v] for u, v in P])
A = csr_matrix((np.ones(len(r), np.int8), (r, c)), shape=(N, N))
match_cols = maximum_bipartite_matching(A, perm_type='column')   # for each ROW, matched col
size_scipy = int((match_cols >= 0).sum())

# ---- networkx Hopcroft-Karp, as in 03_controllability.py
B = nx.Graph()
B.add_nodes_from(['+' + n for n in NODES], bipartite=0)
B.add_nodes_from(['-' + n for n in NODES], bipartite=1)
B.add_edges_from([('+' + u, '-' + v) for u, v in P])
m = hopcroft_karp_matching(B, top_nodes=['+' + n for n in NODES])
matched_in = {k[1:] for k in m if k.startswith('-')}
size_nx = len(matched_in)

ND_scipy = max(N - size_scipy, 1)
ND_nx = max(N - size_nx, 1)

# ---- sanity: the matching returned by networkx really is a matching
pairs = [(k[1:], v[1:]) for k, v in m.items() if k.startswith('+')]
assert len(set(a for a, b in pairs)) == len(pairs), 'a source is matched twice'
assert len(set(b for a, b in pairs)) == len(pairs), 'a target is matched twice'
assert all((a, b) in set(P) for a, b in pairs), 'matched pair is not an edge'

# ---- the 5 nodes that must be drivers in every configuration are the in-degree-0 nodes
indeg = {n: 0 for n in NODES}
for u, v in P:
    indeg[v] += 1
zero_in = sorted(n for n in NODES if indeg[n] == 0)

out = dict(N=N, matching_size_scipy=size_scipy, matching_size_networkx=size_nx,
           agree=bool(size_scipy == size_nx),
           N_D_scipy=ND_scipy, N_D_networkx=ND_nx,
           n_D_fraction=ND_nx / N,
           n_in_degree_zero=len(zero_in), in_degree_zero_nodes=zero_in,
           in_degree_zero_types=[nt[n] for n in zero_in])
json.dump(out, open(f'{RES}/systems_matching_crosscheck.json', 'w'), indent=2)
print(json.dumps(out, indent=2))
