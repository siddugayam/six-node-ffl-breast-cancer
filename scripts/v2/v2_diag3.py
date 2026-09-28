#!/usr/bin/env python3
"""Diagnostics: n=3 core counts under every plausible convention, per graph variant."""
import sys, os, csv, json
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v2_build_graph import build, NODE_TYPE, deposited_edges
from v2_core3 import reduce_graph, bitsets, VARIANTS

nodes = sorted(NODE_TYPE)
NIDX = {n: i for i, n in enumerate(nodes)}

def klass(s, m, comp):
    if (s, m) in comp: return "Composite-FFL"
    a, b = NODE_TYPE[s], NODE_TYPE[m]
    if a == "TF" and b == "miRNA": return "TF-FFL"
    if a == "miRNA" and b == "TF": return "miRNA-FFL"
    return "other(%s->%s)" % (a, b)

def counts(E, induced, reduced):
    E2, comp = reduce_graph(E)
    G = E2 if reduced else E
    idx, out, inn = bitsets(G, nodes)
    full = (1 << len(nodes)) - 1
    cls = Counter()
    for (s, m) in G:
        i, j = idx[s], idx[m]
        if induced:
            if (m, s) in G: continue
            cand = out[i] & out[j] & ~inn[i] & ~inn[j] & full
        else:
            cand = out[i] & out[j]
        k = bin(cand).count("1")
        if k: cls[klass(s, m, comp if reduced else set())] += k
    return cls

def typed(c):
    return sum(v for k, v in c.items() if not k.startswith("other"))

for name, v in VARIANTS.items():
    E, _, _ = build(v)
    for induced in (True, False):
        for reduced in (True, False):
            c = counts(E, induced, reduced)
            print("%-24s induced=%-5s reduced=%-5s total=%7d typed=%6d  comp=%6d tf=%5d mir=%5d" %
                  (name, induced, reduced, sum(c.values()), typed(c),
                   c.get("Composite-FFL", 0), c.get("TF-FFL", 0), c.get("miRNA-FFL", 0)))
    print()
