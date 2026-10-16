#!/usr/bin/env python3
"""Inputs of the tier- and sign-aware null models: graphs in both typings and one attribute line per arc.

Reads the public repository (root: environment variable FFL_ROOT):
  analyses/census_and_motif_nulls/graphs/null_dep_noexemplar.txt (analysed network) and null_pub_noexemplar.txt (census
  graph), with their .meta.json; analyses/six_node_pattern/S1/graph_noexemplar_labels.txt (layer of each census arc:
  0 TF_miRNA, 1 miRNA_target, 2 TF_target analysed, 3 TRRUST layer, 4 gene_gene analysed, 5 STRING, 6 10 kb);
  data/network/edge_evidence_tier.tsv; supplementary_tables/TableS2_all_interactions.csv (trrust_mode);
  data/network/layer_TF_target.tsv (mode); analyses/six_node_pattern/E/E5/null_dep_noexemplar_retyped22.txt (check).
Writes work/<graph>_<typing>.txt (null format) and work/<graph>_<typing>.attr ("tier sign" per arc; tier 0 strong,
1 weak, 2 predicted_only for miRNA -> {Gene, TF} arcs, else -1; sign 0 Activation, 1 Repression, 2 unannotated for arcs
from a TF-typed node to a TF or Gene, else -1), and arc_attribute_counts.csv.
usage: FFL_ROOT=<repository checkout> python3 make_tier_sign_inputs.py"""
import collections, csv, json, os

ROOT = os.environ["FFL_ROOT"]
GR = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "graphs")
HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(HERE, "work"); os.makedirs(W, exist_ok=True)
T22 = ("APEX1 BMI1 BRCA1 CREBBP CTNNB1 EP300 EZH2 HDAC1 HDAC2 HDAC3 HDAC4 HDAC9 ILF3 MEN1 MKL1 MTA1 NCOR1 NF1 RB1 SIRT1 "
       "SUZ12 VHL").split()
TN = {0: "miRNA", 1: "TF", 2: "Gene"}
TIER = {"strong": 0, "weak": 1, "predicted_only": 2}
SIGN = {"Activation": 0, "Repression": 1}

tier = {(r["source"], r["target"]): r["tier"] for r in
        csv.DictReader(open(os.path.join(ROOT, "data", "network", "edge_evidence_tier.tsv")), delimiter="\t")}
s2 = {(r["source"], r["target"]): r["trrust_mode"] for r in
      csv.DictReader(open(os.path.join(ROOT, "supplementary_tables", "TableS2_all_interactions.csv")))
      if r["edge_type"] == "TF_target" and r["in_analysed_network"] == "TRUE"}
assert len(s2) == 770
trrust = {(r["source"], r["target"]): r["mode"] for r in
          csv.DictReader(open(os.path.join(ROOT, "data", "network", "layer_TF_target.tsv")), delimiter="\t")}
labels = [int(x) for x in open(os.path.join(ROOT, "analyses", "six_node_pattern", "S1",
                                            "graph_noexemplar_labels.txt")).read().split()]


def read(fn):
    L = open(os.path.join(GR, fn)).read().rstrip("\n").split("\n")
    nv, ne, nc = map(int, L[0].split()); ty = list(map(int, L[1].split()))
    meta = json.load(open(os.path.join(GR, fn + ".meta.json")))
    arcs = [tuple(map(int, l.split())) for l in L[2:2 + ne]]
    return nv, ne, nc, ty, meta, arcs


counts = []
for gname, fn in (("analysed", "null_dep_noexemplar.txt"), ("census", "null_pub_noexemplar.txt")):
    nv, ne, nc, ty0, meta, arcs0 = read(fn)
    names = meta["nodes"]; ix = {n: i for i, n in enumerate(names)}
    lab = labels if gname == "census" else None
    if lab is not None: assert len(lab) == ne
    for typing in ("published", "retyped"):
        ty = list(ty0)
        if typing == "retyped":
            for x in T22: assert ty[ix[x]] == 1; ty[ix[x]] = 2
            cls = sorted({f"{TN[ty[u]]}->{TN[ty[v]]}" for u, v, _ in arcs0})
            arcs = [(u, v, cls.index(f"{TN[ty[u]]}->{TN[ty[v]]}")) for u, v, _ in arcs0]; ncl = len(cls)
        else:
            arcs, ncl, cls = arcs0, nc, meta["classes"]
        text = f"{nv} {ne} {ncl}\n" + " ".join(map(str, ty)) + "\n" + "\n".join(f"{u} {v} {c}" for u, v, c in arcs) + "\n"
        out = os.path.join(W, f"{gname}_{typing}.txt")
        open(out, "w").write(text)
        if gname == "analysed" and typing == "retyped":
            ref = open(os.path.join(ROOT, "analyses", "six_node_pattern", "E", "E5", "null_dep_noexemplar_retyped22.txt")).read()
            assert text == ref, "retyped analysed network differs from E5"
        if gname == "analysed" and typing == "published":
            assert text == open(os.path.join(GR, fn)).read(), "copy differs"
        attr, c = [], collections.Counter()
        for k, (u, v, _) in enumerate(arcs):
            a, b = names[u], names[v]; su, sv = TN[ty[u]], TN[ty[v]]
            t = s = -1
            if su == "miRNA" and sv in ("Gene", "TF"):
                t = TIER[tier[(a, b)]]
            if su == "TF" and sv in ("Gene", "TF"):
                L = 2 if lab is None else lab[k]
                if L == 2: s = SIGN.get(s2[(a, b)], 2)
                elif L == 3: s = SIGN.get(trrust[(a, b)], 2)
                elif L == 5: s = 2
                else: raise SystemExit(f"TF arc {a}->{b} with layer {L}")
            attr.append(f"{t} {s}")
            c[("tier", t, f"{su}->{sv}")] += t >= 0; c[("sign", s, f"{su}->{sv}")] += s >= 0
        open(os.path.join(W, f"{gname}_{typing}.attr"), "w").write("\n".join(attr) + "\n")
        for (kind, code, cl), n in sorted(c.items()):
            if code >= 0 and n:
                name = {"tier": ["strong", "weak", "predicted_only"], "sign": ["+", "-", "unannotated"]}[kind][code]
                counts.append(dict(graph=gname, typing=typing, attribute=kind, value=name, arc_class=cl, arcs=n))
        print(gname, typing, "classes", cls)
with open(os.path.join(HERE, "arc_attribute_counts.csv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(counts[0])); w.writeheader(); w.writerows(counts)
for r in counts: print(r)
