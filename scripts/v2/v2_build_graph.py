#!/usr/bin/env python3
"""v2 independent graph builder for motif verification.

Builds the canonical directed network from data/canonical_edges.tsv and augments it
with the rebuilt regulatory layers (layer_TF_target.tsv, layer_gene_gene.tsv,
layer_miRNA_miRNA.tsv @ 10 kb).  Several augmentation variants are emitted so that the
sensitivity of the FFL census to the (under-specified) orientation conventions can be
measured.  Nothing here reuses scripts/03_ffl_census.py or scripts/30_motif_significance.py.
"""
import csv, json, os, sys, itertools
from collections import Counter, defaultdict

BASE = "/path/to/revision"
DATA = os.path.join(BASE, "data")

def load_nodes():
    nt = {}
    with open(os.path.join(DATA, "canonical_nodes.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            nt[row["name"]] = row["type"]
    return nt

def load_merge():
    m = json.load(open(os.path.join(DATA, "name_map.json")))
    return m["merge_map"], m["node_type"]

NODE_TYPE = load_nodes()
MERGE, NT2 = load_merge()
assert NODE_TYPE == NT2, "node_type in name_map.json disagrees with canonical_nodes.tsv"

def canon(x):
    return MERGE.get(x, x)

def deposited_edges():
    """(src,dst) -> edge_type from canonical_edges.tsv"""
    E = {}
    n = 0
    with open(os.path.join(DATA, "canonical_edges.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n += 1
            s, t = row["source"], row["target"]
            assert s in NODE_TYPE and t in NODE_TYPE, (s, t)
            assert NODE_TYPE[s] == row["source_type"] and NODE_TYPE[t] == row["target_type"]
            E[(s, t)] = row["edge_type"]
    return E, n

def layer_tf_target():
    out = []
    with open(os.path.join(DATA, "layer_TF_target.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s, t = canon(row["source"]), canon(row["target"])
            if s in NODE_TYPE and t in NODE_TYPE and s != t:
                out.append((s, t))
    return out

def layer_gene_gene():
    directed, undirected = [], []
    with open(os.path.join(DATA, "layer_gene_gene.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s, t = canon(row["source"]), canon(row["target"])
            if not (s in NODE_TYPE and t in NODE_TYPE) or s == t:
                continue
            if row["directed"].strip().upper() == "TRUE":
                directed.append((s, t))
            else:
                undirected.append((s, t))
    return directed, undirected

def layer_mirna_mirna(thr="threshold_10kb"):
    out = []
    with open(os.path.join(DATA, "layer_miRNA_miRNA.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row[thr].strip().upper() != "TRUE":
                continue
            s, t = canon(row["miRNA_1"]), canon(row["miRNA_2"])
            if s in NODE_TYPE and t in NODE_TYPE and s != t:
                out.append((s, t))
    return out

def etype(s, t, fallback):
    a, b = NODE_TYPE[s], NODE_TYPE[t]
    if a == "miRNA" and b in ("Gene", "TF"):
        return "miRNA_target"
    if a == "miRNA" and b == "miRNA":
        return "miRNA_miRNA"
    if a == "TF" and b == "miRNA":
        return "TF_miRNA"
    if a == "TF" and b in ("Gene", "TF"):
        return "TF_target"
    if a == "Gene":
        return "gene_gene" if b == "Gene" else "gene_" + b
    return fallback

def build(variant):
    """variant: dict of switches."""
    E, nraw = deposited_edges()
    prov = {k: "deposited" for k in E}
    added = Counter()
    if variant.get("tf_target"):
        for (s, t) in layer_tf_target():
            if (s, t) not in E:
                E[(s, t)] = etype(s, t, "TF_target"); prov[(s, t)] = "TRRUST_TF_target"; added["TF_target"] += 1
    gg_dir, gg_und = layer_gene_gene()
    if variant.get("gg_directed"):
        for (s, t) in gg_dir:
            if (s, t) not in E:
                E[(s, t)] = etype(s, t, "gene_gene"); prov[(s, t)] = "TRRUST_gene_gene"; added["gg_directed"] += 1
    mode = variant.get("gg_string", "none")
    if mode != "none":
        for (s, t) in gg_und:
            pairs = [(s, t)] if mode == "asis" else [(s, t), (t, s)]
            for (a, b) in pairs:
                if (a, b) not in E:
                    E[(a, b)] = etype(a, b, "gene_gene"); prov[(a, b)] = "STRING"; added["gg_string"] += 1
    mode = variant.get("mirna_mirna", "none")
    if mode != "none":
        for (s, t) in layer_mirna_mirna():
            pairs = [(s, t)] if mode == "asis" else [(s, t), (t, s)]
            for (a, b) in pairs:
                if (a, b) not in E:
                    E[(a, b)] = etype(a, b, "miRNA_miRNA"); prov[(a, b)] = "miRBase_10kb"; added["mirna_mirna"] += 1
    return E, prov, added

def summarise(E):
    c = Counter(E.values())
    recip = set()
    for (s, t) in E:
        if (t, s) in E:
            a, b = NODE_TYPE[s], NODE_TYPE[t]
            if {a, b} == {"TF", "miRNA"}:
                recip.add(frozenset((s, t)))
    return c, recip

if __name__ == "__main__":
    print("nodes:", len(NODE_TYPE), Counter(NODE_TYPE.values()))
    print("raw labels in merge_map:", len(MERGE))
    E0, nraw = deposited_edges()
    print("deposited rows:", nraw, "unique directed:", len(E0))
    c, r = summarise(E0)
    print("  classes:", dict(c))
    print("  reciprocal TF<->miRNA pairs:", len(r))
    tt = layer_tf_target(); gd, gu = layer_gene_gene(); mm = layer_mirna_mirna()
    print("layer_TF_target rows mapped into node universe:", len(tt), "unique:", len(set(tt)))
    print("layer_gene_gene directed:", len(gd), "undirected:", len(gu))
    print("layer_miRNA_miRNA 10kb pairs in universe:", len(mm), "unique:", len(set(mm)))
