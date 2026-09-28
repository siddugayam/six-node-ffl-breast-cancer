#!/usr/bin/env python3
"""Sweep additional augmentation variants (incl. layer_TF_miRNA) looking for the
published 6,037 / 1,434 / 1,649 figures."""
import sys, os, csv, itertools
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2_build_graph as B
from v2_build_graph import NODE_TYPE, deposited_edges, layer_tf_target, layer_gene_gene, layer_mirna_mirna, canon, etype
from v2_core3 import reduce_graph, bitsets
from v2_diag3 import counts, typed

DATA = B.DATA

def layer_tf_mirna():
    out = []
    for row in csv.DictReader(open(os.path.join(DATA, "layer_TF_miRNA.tsv")), delimiter="\t"):
        s, t = canon(row["source"]), canon(row["target"])
        if s in NODE_TYPE and t in NODE_TYPE and s != t:
            out.append((s, t))
    return out

def mk(tt=0, gg=0, mm="none", st="none", tm=0):
    E, _ = deposited_edges()
    if tt:
        for e in layer_tf_target(): E.setdefault(e, etype(*e, "TF_target"))
    if gg:
        d, u = layer_gene_gene()
        for e in d: E.setdefault(e, etype(*e, "gene_gene"))
    if st != "none":
        d, u = layer_gene_gene()
        for (s, t) in u:
            for e in ([(s, t)] if st == "asis" else [(s, t), (t, s)]): E.setdefault(e, etype(*e, "gene_gene"))
    if mm != "none":
        for (s, t) in layer_mirna_mirna():
            for e in ([(s, t)] if mm == "asis" else [(s, t), (t, s)]): E.setdefault(e, etype(*e, "miRNA_miRNA"))
    if tm:
        for e in layer_tf_mirna(): E.setdefault(e, etype(*e, "TF_miRNA"))
    return E

print("tf_mirna layer edges in universe:", len(set(layer_tf_mirna())))
res = []
for tt in (0, 1):
 for gg in (0, 1):
  for tm in (0, 1):
   for mm in ("none", "asis", "both"):
    for st in ("none", "asis", "both"):
     E = mk(tt, gg, mm, st, tm)
     E2, comp = reduce_graph(E)
     for induced in (True, False):
      for reduced in (True,):
       c = counts(E, induced, reduced)
       tot = sum(c.values()); ty = typed(c)
       res.append((tot, ty, c.get("Composite-FFL",0), c.get("TF-FFL",0), c.get("miRNA-FFL",0),
                   "tt=%d gg=%d tm=%d mm=%s st=%s ind=%s" % (tt,gg,tm,mm,st,induced), len(E), len(comp)))
for r in sorted(res):
    print("total=%7d typed=%6d comp=%6d tf=%5d mir=%5d  |E|=%5d recip=%4d  %s" % (r[0],r[1],r[2],r[3],r[4],r[6],r[7],r[5]))
print()
print("--- any hitting 6037 / 1434 / 1649 ---")
for r in res:
    if r[0] in (6037,1649,1434) or r[1] in (6037,1649,1434) or r[2] in (1434,) :
        print(r)
