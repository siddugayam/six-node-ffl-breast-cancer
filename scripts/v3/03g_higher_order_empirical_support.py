#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03g_higher_order_empirical_support.py
=====================================
HOW MUCH OF THE NETWORK ACTUALLY CARRIES THE EDGES THAT THE 4-, 5- AND 6-NODE
ARCHITECTURES ARE BUILT FROM?

The dynamical escalation in Part C adds, in order, a gene-gene edge, a
miRNA-miRNA co-transcription edge, and a TF-TF edge.  Before asking what those
additions DO, it is worth stating how often they EXIST.  This script counts them
directly in data/canonical_edges.tsv and in the enumerated higher-order module
listing (results/ffl_higher_order.csv), and reports the answer whichever way it
falls.

Out: results/v3/dynamics_higher_order_empirical_support.csv
"""
import numpy as np, pandas as pd

REV = "/path/to/revision"
RES = f"{REV}/results/v3"

E = pd.read_csv(f"{REV}/data/canonical_edges.tsv", sep="\t")
G = pd.read_csv(f"{REV}/results/graph_used_edges.tsv", sep="\t")
L = pd.read_csv(f"{REV}/data/layer_gene_gene.tsv", sep="\t")
H = pd.read_csv(f"{REV}/results/ffl_higher_order.csv")

rows = []


def add(**kw):
    rows.append(kw)


# ---- the project uses TWO edge sets on the SAME 587 nodes; state both ---------
add(quantity="edges in the SIGNED canonical network (coherence analysis)", value=len(E),
    detail="data/canonical_edges.tsv; every edge carries a sign")
add(quantity="edges in the graph used for the higher-order MODULE ENUMERATION",
    value=len(G), detail="results/graph_used_edges.tsv; a superset of the signed network "
    f"({len(G)-len(E)} additional edges, mostly unsigned gene-gene and miRNA-miRNA pairs)")
ne = set(E.source) | set(E.target); ng = set(G.source) | set(G.target)
add(quantity="nodes shared by the two edge sets", value=len(ne & ng),
    detail=f"signed network {len(ne)} nodes; enumeration graph {len(ng)} nodes; "
           f"{len(ng - ne)} nodes are unique to the enumeration graph")
for t, n in E.edge_type.value_counts().items():
    add(quantity=f"SIGNED network edges of type {t}", value=int(n), detail="")
for t, n in G.edge_type.value_counts().items():
    add(quantity=f"ENUMERATION graph edges of type {t}", value=int(n), detail="")

gg = E[E.edge_type == "gene_gene"]
add(quantity="gene-gene edges in the SIGNED canonical network", value=len(gg),
    detail="the 4-node step of the dynamical escalation has only this much signed support: "
    + "; ".join(f"{r.source}->{r.target} (sign {r.sign:+d}, present in {r.in_networks})"
                for r in gg.itertuples()))
add(quantity="gene-gene edges in the ENUMERATION graph",
    value=int((G.edge_type == "gene_gene").sum()),
    detail="unsigned; these are what the 4/5/6-node module listing actually uses")
add(quantity="candidate gene-gene edges in the source layer before network restriction",
    value=len(L), detail="data/layer_gene_gene.tsv: "
    + "; ".join(f"{k}={v}" for k, v in L.evidence.value_counts().items()))
add(quantity="TF->TF edges in the SIGNED canonical network",
    value=int(((E.edge_type == "TF_target") & (E.target_type == "TF")).sum()),
    detail="TF_target edges whose target is itself a TF")
add(quantity="miRNA-miRNA edges in the SIGNED canonical network",
    value=int((E.edge_type == "miRNA_miRNA").sum()),
    detail="genomic-cluster co-transcription pairs that survived into the signed network")
add(quantity="miRNA-miRNA edges in the ENUMERATION graph",
    value=int((G.edge_type == "miRNA_miRNA").sum()), detail="unsigned")

for n in (4, 5, 6):
    s = H[H.n_nodes == n]
    complete = (s.listing == "complete").all()
    add(quantity=f"{n}-node modules listed in results/ffl_higher_order.csv", value=int(len(s)),
        detail=("complete enumeration" if complete else
                "capped at 5,000 per architecture, so counts below are lower bounds"))
    for lab, pat in (("gene-gene", "gene-gene"), ("miRNA-miRNA", "miRNA-miRNA"), ("TF-TF", "TF-TF")):
        k = int(s.architecture.str.contains(pat).sum())
        add(quantity=f"{n}-node modules whose architecture contains a {lab} edge",
            value=k, detail=f"{100*k/max(len(s),1):.2f}% of the {n}-node modules listed")

# which nodes do the gene-gene-containing modules use?
sub = H[H.architecture.str.contains("gene-gene")]
if len(sub):
    genes = pd.Series([m for row in sub.members for m in row.split(";")]).value_counts()
    add(quantity="distinct nodes appearing in gene-gene-containing higher-order modules",
        value=int(len(genes)),
        detail="commonest: " + ", ".join(f"{g} ({c})" for g, c in genes.head(6).items()))
    both = sub.members.str.contains("COL1A1") & sub.members.str.contains("COL3A1")
    add(quantity="gene-gene-containing modules built on the ONE signed gene-gene edge "
                 "(COL1A1-COL3A1)", value=int(both.sum()),
        detail=f"{100*both.mean():.2f}% of the {len(sub)} gene-gene-containing modules; "
               "the rest use gene-gene edges that carry no sign in the canonical network")

R = pd.DataFrame(rows)
R.to_csv(f"{RES}/dynamics_higher_order_empirical_support.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 110)
print(R.to_string(index=False))
