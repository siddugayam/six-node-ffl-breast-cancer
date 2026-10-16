#!/usr/bin/env python3
"""Independent check of census_reciprocal_variants.c at small n.

Enumerates every weakly connected induced n-vertex subgraph of the census graph with undirected ESU (Wernicke 2006) and
tests it with a separate implementation of the conditions (paper numbering): D2 acyclic; D3 one source, one sink, every
vertex on a source-sink path; D4 two internally vertex-disjoint source-sink paths (Menger: after removing a direct
source->sink arc, which is one path, the sink must still be reachable; otherwise no single internal vertex may separate
source from sink). Counts per variant:
  B  contracted census graph, D2-D4          U  uncontracted census graph, D2-D4
  M  uncontracted, SCCs of the induced subgraph contracted, then D3-D4 on the condensed module
  S  modules of the condensation DAG of the uncontracted graph (vertices = SCCs)
and the overlaps B&U, B&M. Prints one line per n; the summary script compares them with the C counts.
usage: FFL_ROOT=<repository checkout> python3 check_small_n.py 3 4
(the test functions are also imported by six_node_pattern_variants.py)"""
import os, sys, itertools
import networkx as nx

ROOT = os.environ["FFL_ROOT"]
G_DIR = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "graphs")


def load(fn):
    L = open(os.path.join(G_DIR, fn)).read().split("\n")
    n, m, _ = map(int, L[0].split())
    out = [set() for _ in range(n)]
    for line in L[2:2 + m]:
        a, b, _ = map(int, line.split()); out[a].add(b)
    return n, out


def reach_from(s, out, S, banned=()):
    seen = {s}; st = [s]
    while st:
        u = st.pop()
        for v in out[u]:
            if v in S and v not in seen and v not in banned: seen.add(v); st.append(v)
    return seen


def d234(S, out):
    """returns 0 if D2-D4 hold on the induced subgraph, else the first failing condition (2, 3 or 4)"""
    S = set(S)
    indeg = {v: sum(1 for u in S if v in out[u]) for v in S}
    order = [v for v in S if indeg[v] == 0]; k = 0
    while k < len(order):
        u = order[k]; k += 1
        for v in out[u]:
            if v in S:
                indeg[v] -= 1
                if indeg[v] == 0: order.append(v)
    if len(order) != len(S): return 2
    return d34(S, out)


def d34(S, out):
    src = [v for v in S if not any(v in out[u] for u in S)]
    snk = [v for v in S if not any(w in S for w in out[v])]
    if len(src) != 1 or len(snk) != 1 or src[0] == snk[0]: return 3
    s, t = src[0], snk[0]
    fw = reach_from(s, out, S)
    back = {t}; st = [t]
    while st:
        u = st.pop()
        for v in S:
            if u in out[v] and v not in back: back.add(v); st.append(v)
    if fw != S or back != S: return 3
    if t in out[s]:
        out2 = {v: (out[v] - {t} if v == s else out[v]) for v in S}
        return 0 if t in reach_from(s, out2, S) else 4
    for v in S - {s, t}:
        if t not in reach_from(s, out, S, banned=(v,)): return 4
    return 0


def condensed(S, out):
    H = nx.DiGraph(); H.add_nodes_from(S)
    H.add_edges_from((u, v) for u in S for v in out[u] if v in S)
    C = nx.condensation(H)
    cout = {c: set(C.successors(c)) for c in C.nodes()}
    return set(C.nodes()), cout


def mtest(S, out):
    Sc, cout = condensed(S, out)
    if len(Sc) == 1: return 3
    return d34(Sc, cout)


def esu(n_nodes, und, k, cb):
    for v in range(n_nodes):
        ext = {u for u in und[v] if u > v}
        _ext([v], ext, v, und, k, cb)


def _ext(sub, ext, v, und, k, cb):
    if len(sub) == k: cb(sub); return
    ext = set(ext)
    while ext:
        w = ext.pop()
        nb = set(sub) | {u for s in sub for u in und[s]}
        _ext(sub + [w], ext | {u for u in und[w] if u > v and u not in nb}, v, und, k, cb)


if __name__ == "__main__":
    n, cout_ = load("census_noexemplar_fine.txt")
    n2, uout = load("census_noexemplar_uncontracted.txt")
    und = [set() for _ in range(n)]
    for a in range(n):
        for b in uout[a]: und[a].add(b); und[b].add(a)
    for a in range(n):
        for b in cout_[a]: assert b in und[a]
    # condensation DAG of the uncontracted graph
    G = nx.DiGraph(); G.add_nodes_from(range(n)); G.add_edges_from((a, b) for a in range(n) for b in uout[a])
    CD = nx.condensation(G)
    cd_out = [set(CD.successors(c)) for c in range(CD.number_of_nodes())]
    cd_und = [set(CD.successors(c)) | set(CD.predecessors(c)) for c in range(CD.number_of_nodes())]
    for k in map(int, sys.argv[1:]):
        cnt = dict(B=0, U=0, M=0, BU=0, BM=0, S=0)

        def cb(sub):
            b = d234(sub, cout_) == 0; u = d234(sub, uout) == 0; m = mtest(sub, uout) == 0
            cnt["B"] += b; cnt["U"] += u; cnt["M"] += m; cnt["BU"] += b and u; cnt["BM"] += b and m
        esu(n, und, k, cb)

        def cbs(sub):
            cnt["S"] += d234(sub, cd_out) == 0
        esu(CD.number_of_nodes(), cd_und, k, cbs)
        print(f"n={k} " + " ".join(f"{a}={b}" for a, b in cnt.items()), flush=True)
