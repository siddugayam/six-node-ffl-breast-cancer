#!/usr/bin/env python3
"""Strongly connected components of the census graph and its condensation DAG.

Reads the uncontracted census graph (both arcs of each reciprocal TF<->miRNA pair kept; 9,226 arcs) and the contracted
one (8,003 arcs) from analyses/census_and_motif_nulls/graphs/ of the public repository (root: environment variable
FFL_ROOT). Writes:
  work/scc_uncontracted.txt   "NV NC NA" / component id of every vertex / NA arcs of the condensation DAG
                              (components numbered by their smallest vertex index), read by census_reciprocal_variants.c
  scc_structure.csv           size and node-type composition of every non-singleton SCC, and a summary row per graph
usage: FFL_ROOT=<repository checkout> python3 make_scc_file.py"""
import collections, csv, os
import networkx as nx

ROOT = os.environ["FFL_ROOT"]
G_DIR = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "graphs")
HERE = os.path.dirname(os.path.abspath(__file__))
TYPES = {0: "miRNA", 1: "TF", 2: "Gene"}


def load(fn):
    L = open(os.path.join(G_DIR, fn)).read().split("\n")
    n, m, _ = map(int, L[0].split()); ty = list(map(int, L[1].split()))
    G = nx.DiGraph(); G.add_nodes_from(range(n))
    for line in L[2:2 + m]:
        a, b, _ = map(int, line.split()); G.add_edge(a, b)
    assert G.number_of_edges() == m
    return G, ty


os.makedirs(os.path.join(HERE, "work"), exist_ok=True)
rows = []
for tag, fn in (("uncontracted", "census_noexemplar_uncontracted.txt"), ("contracted", "census_noexemplar_fine.txt")):
    G, ty = load(fn)
    comps = sorted((sorted(c) for c in nx.strongly_connected_components(G)), key=lambda c: c[0])
    cid = {v: i for i, c in enumerate(comps) for v in c}
    arcs = sorted({(cid[a], cid[b]) for a, b in G.edges() if cid[a] != cid[b]})
    big = [c for c in comps if len(c) > 1]
    for c in big:
        tc = collections.Counter(TYPES[ty[v]] for v in c)
        rows.append(dict(graph=tag, component_smallest_vertex=c[0], size=len(c), TF=tc["TF"], miRNA=tc["miRNA"],
                         Gene=tc["Gene"], note=""))
    tall = collections.Counter(TYPES[t] for t in ty)
    rows.append(dict(graph=tag, component_smallest_vertex="all", size=G.number_of_nodes(), TF=tall["TF"],
                     miRNA=tall["miRNA"], Gene=tall["Gene"],
                     note=f"{G.number_of_edges()} arcs; {len(comps)} SCCs ({len(big)} non-singleton); "
                          f"condensation DAG {len(comps)} vertices, {len(arcs)} arcs"))
    if tag == "uncontracted":
        with open(os.path.join(HERE, "work", "scc_uncontracted.txt"), "w") as fh:
            fh.write(f"{G.number_of_nodes()} {len(comps)} {len(arcs)}\n")
            fh.write(" ".join(str(cid[v]) for v in range(G.number_of_nodes())) + "\n")
            for a, b in arcs: fh.write(f"{a} {b}\n")
with open(os.path.join(HERE, "scc_structure.csv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
