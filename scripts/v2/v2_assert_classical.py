#!/usr/bin/env python3
"""Assert that at n=3 the formal definition D1-D4 reduces EXACTLY to the classical FFL
(R->M, R->T, M->T and no other arc among the three).  Full brute force in Python, set-equality
against the induced-transitive-triangle set, independent of the C enumerator."""
import sys, os, itertools
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v2_graph_final import make, reduce_graph
from v2_build_graph import NODE_TYPE

def d1d4(V, adj):
    arcs = [(x, y) for x in V for y in V if x != y and y in adj.get(x, ())]
    A = set(arcs)
    outd = defaultdict(int); ind = defaultdict(int)
    for x, y in arcs: outd[x] += 1; ind[y] += 1
    # D1 acyclic (graph already contracted)
    order = []; rem = dict((v, ind[v]) for v in V); left = set(V)
    while True:
        z = [v for v in left if rem[v] == 0]
        if not z: break
        for v in z:
            left.discard(v); order.append(v)
            for y in V:
                if (v, y) in A: rem[y] -= 1
    if left: return None
    # D2
    src = [v for v in V if ind[v] == 0]; snk = [v for v in V if outd[v] == 0]
    if len(src) != 1 or len(snk) != 1 or src[0] == snk[0]: return None
    s, t = src[0], snk[0]
    # D3 every vertex on an s->t path
    def reach(start, fwd):
        seen = {start}; st = [start]
        while st:
            u = st.pop()
            for w in V:
                e = (u, w) if fwd else (w, u)
                if e in A and w not in seen: seen.add(w); st.append(w)
        return seen
    if not (set(V) <= reach(s, True) & reach(t, False)): return None
    # D4 >=2 internally vertex-disjoint s->t paths  (enumerate all simple paths, n=3 only)
    paths = []
    def dfs(u, path):
        if u == t: paths.append(tuple(path)); return
        for w in V:
            if (u, w) in A and w not in path: dfs(w, path + [w])
    dfs(s, [s])
    ok = False
    for p1, p2 in itertools.combinations(paths, 2):
        if not (set(p1[1:-1]) & set(p2[1:-1])): ok = True; break
    if not ok: return None
    return (s, t, frozenset(A))

E = make("pub"); E2, comp = reduce_graph(E)
adj = defaultdict(set); und = defaultdict(set)
for (u, v) in E2:
    adj[u].add(v); und[u].add(v); und[v].add(u)

seen = set(); modules = {}
for a in und:
    na = und[a]
    for b in na:
        for c in (na | und[b]):
            if c == a or c == b: continue
            key = tuple(sorted((a, b, c)))
            if key in seen: continue
            seen.add(key)
            r = d1d4(list(key), adj)
            if r: modules[key] = r
print("connected triples examined:", len(seen))
print("D1-D4 modules at n=3:", len(modules))

# classical set: induced transitive triangles
classical = set()
for (R, M) in E2:
    if (M, R) in E2: continue
    for T in adj[R] & adj[M]:
        if R in adj.get(T, ()) or M in adj.get(T, ()): continue
        classical.add(tuple(sorted((R, M, T))))
print("classical FFL (R->M,R->T,M->T, induced) node-sets:", len(classical))
print("SET EQUALITY D1-D4 == classical :", set(modules) == classical)
bad = 0
for k, (s, t, A) in modules.items():
    if len(A) != 3: bad += 1
print("modules whose induced arc set is not exactly 3 arcs:", bad)
