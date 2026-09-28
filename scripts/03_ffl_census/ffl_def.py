"""
ffl_def.py -- INDEPENDENT reference implementation of the formal n-node FFL
definition.  Written separately from the C enumerator (scripts/03_ffl_census/ffl_enum.c) so
that the two can be cross-validated against each other.

An n-node FFL is a vertex-induced subgraph S on n nodes of the directed graph G
for which there exists an orientation of the reciprocal (mutual) arc pairs such
that
  (i)   S contains a 3-node FFL core R->M, R->T, M->T
  (ii)  S is a DAG
  (iii) S has exactly one source and exactly one sink
  (iv)  >=2 internally vertex-disjoint directed source->sink paths exist
Every directed path (the core's 2-step arm included) may traverse at most one
UNDIRECTED edge (miRNA-miRNA co-transcription; STRING association if enabled).

Orientations are enumerated as topological orders: fixing a linear order of S
with the source first and the sink last and orienting every edge forward
recovers exactly the set of admissible orientations.
"""
from itertools import permutations, combinations


def simple_paths(order_pos, fw, fwu, s, t):
    """all simple s->t paths in the oriented DAG; yields frozenset of internal nodes"""
    res = []
    stack = [(s, [], 0)]
    while stack:
        cur, internal, nund = stack.pop()
        for nxt in fw[cur]:
            k = nund + (1 if nxt in fwu[cur] else 0)
            if k > 1:
                continue
            if nxt == t:
                res.append(frozenset(internal))
            else:
                stack.append((nxt, internal + [nxt], k))
    return res


def test_order(S, arcs_in, und_in, order):
    """orient every edge forward along `order`; test (i)-(iv). -> (ok, ndisj)"""
    pos = {v: i for i, v in enumerate(order)}
    fw = {v: set() for v in S}
    bw = {v: set() for v in S}
    fwu = {v: set() for v in S}
    for u in S:
        for v in S:
            if u == v or (u, v) not in arcs_in:
                continue
            if pos[u] < pos[v]:
                fw[u].add(v); bw[v].add(u)
                if (u, v) in und_in:
                    fwu[u].add(v)
            elif (v, u) not in arcs_in:
                return False, 0          # strict arc runs backwards -> cyclic
    s, t = order[0], order[-1]
    if bw[s] or fw[t]:
        return False, 0
    for v in S:
        if v == s:
            if not fw[v]: return False, 0
        elif v == t:
            if not bw[v]: return False, 0
        else:
            if not fw[v] or not bw[v]: return False, 0
    # (i) 3-node core
    core = None
    for R in S:
        for T in fw[R]:
            for M in fw[R] & bw[T]:
                nu = (1 if M in fwu[R] else 0) + (1 if T in fwu[M] else 0)
                if nu <= 1:
                    core = (R, M, T); break
            if core: break
        if core: break
    if core is None:
        return False, 0
    # (iv) internally vertex-disjoint source->sink paths
    paths = simple_paths(pos, fw, fwu, s, t)
    best = 0
    for r in range(min(len(paths), len(S)), 1, -1):
        for comb in combinations(paths, r):
            ok = True
            for a in range(r):
                for b in range(a + 1, r):
                    if comb[a] & comb[b]:
                        ok = False; break
                if not ok: break
            if ok:
                best = r; break
        if best: break
    if best < 2:
        return False, 0
    return True, best


def test_set(S, arcs, und):
    """S: iterable of node names.  Returns dict or None."""
    S = list(S)
    n = len(S)
    arcs_in = set((u, v) for u in S for v in S if u != v and (u, v) in arcs)
    und_in = set(e for e in arcs_in if e in und)
    best = None
    for s in S:
        if any((v, s) in arcs_in and (s, v) not in arcs_in for v in S):
            continue                       # strict arc into s -> s not a source
        for t in S:
            if t == s: continue
            if any((t, v) in arcs_in and (v, t) not in arcs_in for v in S):
                continue
            mid = [v for v in S if v != s and v != t]
            for perm in permutations(mid):
                ok, nd = test_order(S, arcs_in, und_in, (s,) + perm + (t,))
                if ok:
                    cand = dict(source=s, sink=t, ndisj=nd, order=(s,) + perm + (t,))
                    if best is None or (cand["source"], cand["sink"]) < (best["source"], best["sink"]):
                        best = cand
                    break
    return best


def failing_condition(S, arcs, und):
    """diagnose WHICH condition a node set fails (first that fails)."""
    S = list(S)
    arcs_in = set((u, v) for u in S for v in S if u != v and (u, v) in arcs)
    und_in = set(e for e in arcs_in if e in und)
    # (i) any orientable transitive triangle?
    tri = False
    for a, b, c in permutations(S, 3):
        if (a, b) in arcs_in and (b, c) in arcs_in and (a, c) in arcs_in:
            nu = (1 if (a, b) in und_in else 0) + (1 if (b, c) in und_in else 0)
            if nu <= 1:
                tri = True; break
    # connectivity (weak)
    adj = {v: set() for v in S}
    for (u, v) in arcs_in:
        adj[u].add(v); adj[v].add(u)
    seen = {S[0]}; st = [S[0]]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    connected = len(seen) == len(S)
    # (ii)+(iii): any acyclic orientation with a unique source and unique sink?
    ok_iii = False
    for s in S:
        for t in S:
            if t == s: continue
            mid = [v for v in S if v != s and v != t]
            for perm in permutations(mid):
                order = (s,) + perm + (t,)
                pos = {v: i for i, v in enumerate(order)}
                bad = False
                fw = {v: set() for v in S}; bw = {v: set() for v in S}
                for (u, v) in arcs_in:
                    if pos[u] < pos[v]: fw[u].add(v); bw[v].add(u)
                    elif (v, u) not in arcs_in: bad = True; break
                if bad: continue
                if bw[s] or fw[t]: continue
                if all((fw[v] and bw[v]) for v in S if v not in (s, t)) and fw[s] and bw[t]:
                    ok_iii = True; break
            if ok_iii: break
        if ok_iii: break
    res = test_set(S, arcs, und)
    fails = []
    if not connected: fails.append("connected")
    if not tri: fails.append("(i) no 3-node FFL core")
    if not ok_iii: fails.append("(ii)/(iii) no acyclic orientation with a unique source and unique sink")
    if res is None and not fails:
        fails.append("(iv) fewer than two vertex-disjoint source->sink paths")
    return fails, res
