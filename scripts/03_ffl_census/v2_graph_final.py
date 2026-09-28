#!/usr/bin/env python3
"""Final graph definitions used by all v2 motif work.

PRIMARY ("pub"): the augmentation convention that exactly reproduces the published
exhaustive n=3 count of 6,037 -
   deposited canonical_edges.tsv
 + layer_TF_target.tsv (TRRUST; supersedes the 696 directed rows of layer_gene_gene.tsv)
 + layer_gene_gene.tsv STRING rows, kept in the single orientation listed in the file
 + layer_miRNA_miRNA.tsv @10 kb, added in BOTH orientations (co-transcription is symmetric)

Sensitivity graphs:
  "dep"       deposited edges only
  "sym"       as PRIMARY but STRING also symmetrised (both orientations)
  "nomm"      as PRIMARY but without the miRNA-miRNA co-transcription layer
"""
import sys, os, csv
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v2_build_graph import (deposited_edges, layer_tf_target, layer_gene_gene,
                            layer_mirna_mirna, etype, NODE_TYPE)

def make(name="pub"):
    E, _ = deposited_edges()
    if name == "dep":
        return E
    for e in layer_tf_target():
        E.setdefault(e, etype(*e, "TF_target"))
    gd, gu = layer_gene_gene()
    for e in gd:
        E.setdefault(e, etype(*e, "gene_gene"))
    if name in ("pub", "nomm"):
        for e in gu:
            E.setdefault(e, etype(*e, "gene_gene"))
    elif name == "sym":
        for (s, t) in gu:
            E.setdefault((s, t), etype(s, t, "gene_gene"))
            E.setdefault((t, s), etype(t, s, "gene_gene"))
    if name != "nomm":
        for (s, t) in layer_mirna_mirna():
            E.setdefault((s, t), "miRNA_miRNA")
            E.setdefault((t, s), "miRNA_miRNA")
    return E

def reduce_graph(E):
    comp = set(); drop = set()
    for (s, t) in E:
        if (t, s) in E and {NODE_TYPE[s], NODE_TYPE[t]} == {"TF", "miRNA"}:
            if NODE_TYPE[s] == "TF": comp.add((s, t))
            else: drop.add((s, t))
    return {k: v for k, v in E.items() if k not in drop}, comp

def fine_class(s, t):
    return NODE_TYPE[s] + "->" + NODE_TYPE[t]

if __name__ == "__main__":
    for nm in ("dep", "pub", "sym", "nomm"):
        E = make(nm); E2, comp = reduce_graph(E)
        print(nm, "|E|=", len(E), "|E_reduced|=", len(E2), "composite arcs=", len(comp))
        print("   coarse:", dict(Counter(E.values())))
        print("   fine (reduced graph):", dict(Counter(fine_class(*e) for e in E2)))
