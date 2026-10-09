#!/usr/bin/env python3
"""Test the whole six-node exemplar circuit (its 10 molecules) against conditions D1-D4 (Methods 2.2) in two
graphs, from the deposited files of the public repository only:
  (1) the analysed network (Table S2 rows with in_analysed_network TRUE), each reciprocal TF<->miRNA pair replaced by its
      TF->miRNA arc;
  (2) the census graph: (1) plus the TRRUST arcs (data/network/layer_TF_target.tsv), the gene-gene layer
      (data/network/layer_gene_gene.tsv; STRING associations in the orientation listed) and the miRNA pairs co-clustered at
      10 kb (data/network/layer_miRNA_miRNA.tsv, threshold_10kb), each pair entered as two opposed arcs.
D1 induced and weakly connected; D2 acyclic; D3 one source, one sink, every vertex on a source-sink path; D4 at least two
internally vertex-disjoint source-sink paths. Prints the result and any directed cycle. No external packages.
Usage: python3 scripts/03_ffl_census/02i_exemplar_whole_circuit_check.py [repository root]"""
import csv, json, os, sys
from collections import defaultdict

R = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mm = json.load(open(f"{R}/data/network/name_map.json"))["merge_map"]
sif = [l.rstrip("\n").split("\t") for l in open(f"{R}/analyses/original_submission_code/SIF_files/6-TF.sif") if l.strip()]
V = sorted({mm.get(r[0], r[0]) for r in sif} | {mm.get(r[-1], r[-1]) for r in sif})
ntype = {r["name"]: r["type"] for r in csv.DictReader(open(f"{R}/data/network/canonical_nodes.tsv"), delimiter="\t")}
print("circuit molecules (%d):" % len(V), ", ".join(V))

net = set(); deposited = set()
for r in csv.DictReader(open(f"{R}/supplementary_tables/TableS2_all_interactions.csv")):
    if r["source"] in V and r["target"] in V:
        deposited.add((r["source"], r["target"]))
        if r["in_analysed_network"].upper() == "TRUE": net.add((r["source"], r["target"]))


def collapse(arcs):   # reciprocal TF<->miRNA pair -> TF->miRNA arc
    return {(a, b) for a, b in arcs if not (ntype.get(a) == "miRNA" and ntype.get(b) == "TF" and (b, a) in arcs)}


census = set(net)
for r in csv.DictReader(open(f"{R}/data/network/layer_TF_target.tsv"), delimiter="\t"):
    if r["source"] in V and r["target"] in V: census.add((r["source"], r["target"]))
for r in csv.DictReader(open(f"{R}/data/network/layer_gene_gene.tsv"), delimiter="\t"):
    if r["source"] in V and r["target"] in V: census.add((r["source"], r["target"]))
for r in csv.DictReader(open(f"{R}/data/network/layer_miRNA_miRNA.tsv"), delimiter="\t"):
    if r["threshold_10kb"].upper() == "TRUE" and r["miRNA_1"] in V and r["miRNA_2"] in V:
        census |= {(r["miRNA_1"], r["miRNA_2"]), (r["miRNA_2"], r["miRNA_1"])}


def cycle(arcs):
    adj = defaultdict(list)
    for a, b in sorted(arcs): adj[a].append(b)      # sorted: the same cycle is reported on every run
    color, stack = {}, []
    def dfs(u):
        color[u] = 1; stack.append(u)
        for w in adj[u]:
            if color.get(w) == 1: return stack[stack.index(w):] + [w]
            if not color.get(w):
                c = dfs(w)
                if c: return c
        color[u] = 2; stack.pop(); return None
    for v in V:
        if not color.get(v):
            c = dfs(v)
            if c: return c
    return None


def disjoint_paths(arcs, s, t):   # unit vertex capacities by node splitting, augmenting paths (BFS)
    cap = defaultdict(int); adj = defaultdict(set)
    def add(u, v, c): cap[(u, v)] += c; adj[u].add(v); adj[v].add(u)
    for v in V: add((v, "in"), (v, "out"), 1 if v not in (s, t) else 10**6)
    for a, b in arcs: add((a, "out"), (b, "in"), 1)
    src, snk, flow = (s, "out"), (t, "in"), 0
    while True:
        prev, q = {src: None}, [src]
        while q and snk not in prev:
            u = q.pop(0)
            for w in adj[u]:
                if w not in prev and cap[(u, w)] > 0: prev[w] = u; q.append(w)
        if snk not in prev: return flow
        w = snk
        while prev[w] is not None: u = prev[w]; cap[(u, w)] -= 1; cap[(w, u)] += 1; w = u
        flow += 1


def test(name, arcs):
    arcs = {(a, b) for a, b in arcs if a in V and b in V and a != b}
    und = defaultdict(set)
    for a, b in arcs: und[a].add(b); und[b].add(a)
    seen, q = {V[0]}, [V[0]]
    while q:
        u = q.pop()
        for w in und[u] - seen: seen.add(w); q.append(w)
    d1 = len(seen) == len(V)
    cyc = cycle(arcs); d2 = cyc is None
    indeg = {v: sum(1 for a, b in arcs if b == v) for v in V}; outdeg = {v: sum(1 for a, b in arcs if a == v) for v in V}
    srcs = [v for v in V if indeg[v] == 0]; snks = [v for v in V if outdeg[v] == 0]
    d3 = d4 = False; npaths = None
    if len(srcs) == 1 and len(snks) == 1:
        s, t = srcs[0], snks[0]
        fwd = defaultdict(set); bwd = defaultdict(set)
        for a, b in arcs: fwd[a].add(b); bwd[b].add(a)
        def reach(x, g):
            r, q = {x}, [x]
            while q:
                u = q.pop()
                for w in g[u] - r: r.add(w); q.append(w)
            return r
        d3 = reach(s, fwd) & reach(t, bwd) == set(V)
        if d2: npaths = disjoint_paths(arcs, s, t); d4 = npaths >= 2
    print(f"\n{name}: {len(arcs)} arcs; D1 connected {d1}; D2 acyclic {d2}"
          + (f" (cycle: {' -> '.join(cyc)})" if cyc else "")
          + f"; sources {srcs}; sinks {snks}; D3 {d3}; D4 {d4}" + (f" ({npaths} disjoint paths)" if npaths is not None else ""))
    print(f"  satisfies D1-D4: {d1 and d2 and d3 and d4}")
    if d1 and d2 and d3 and d4:
        no = arcs - {("COL1A1", "COL3A1")}
        print("  without COL1A1->COL3A1:", "still satisfies" if (cycle(no) is None and
              len([v for v in V if sum(1 for a, b in no if a == v) == 0]) == 1) else "fails (no single sink)")


test("analysed network as given", net)
test("analysed network, reciprocal pairs replaced by the TF->miRNA arc", collapse(net))
test("census graph (TRRUST, STRING, 10-kb co-transcription added; reciprocal pairs replaced)", collapse(census))
test("deposited network (with the 30 miRNA-miRNA exemplar edges), reciprocal pairs replaced", collapse(deposited))
