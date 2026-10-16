#!/usr/bin/env python3
"""Census graphs with the STRING associations re-oriented, and the gene-gene layer variants.

Re-implements build_graph() of scripts/03_ffl_census/03_ffl_census_noexemplarmirna.py from the public network files
(root: environment variable FFL_ROOT): analysed network (data/network/canonical_edges.tsv without the 30 exemplar
miRNA-miRNA edges), then the TRRUST layer (layer_TF_target.tsv), the gene-gene file (layer_gene_gene.tsv: STRING and
TRRUST rows) and the 10 kb miRNA pairs (layer_miRNA_miRNA.tsv, both directions), each arc added only if absent; finally
each reciprocal TF<->miRNA pair is replaced by its TF->miRNA arc.

Checks (the script stops if one fails): listed orientation == analyses/census_and_motif_nulls/graphs/
census_noexemplar_fine.txt arc for arc; no STRING == the arcs of census_noexemplar_nostring_edgeclass.txt.
Writes work/graphs/: listed.txt, nostring.txt, nogenegene.txt, boot_000.txt ... boot_499.txt (the census "fine"
format), work/string_pairs.txt (the 941 associations as vertex pairs), work/string_arcs_listed.txt (the STRING arcs of
the listed-orientation graph) and string_bootstrap_graphs.csv (arcs per replicate).
usage: FFL_ROOT=<repository checkout> python3 string_orientation_graphs.py [replicates, default 500]"""
import csv, os, sys
import numpy as np

ROOT = os.environ["FFL_ROOT"]
NET = os.path.join(ROOT, "data", "network")
GR = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "graphs")
HERE = os.path.dirname(os.path.abspath(__file__))
R = int(sys.argv[1]) if len(sys.argv) > 1 else 500
SEED0 = 20261015
TMAP = {"miRNA": 0, "TF": 1, "Gene": 2}

rd = lambda f: list(csv.DictReader(open(os.path.join(NET, f)), delimiter="\t"))
nodes = {r["name"]: r["type"] for r in rd("canonical_nodes.tsv")}
V = sorted(nodes); idx = {v: i for i, v in enumerate(V)}
CANON = [r for r in rd("canonical_edges.tsv") if r["edge_type"] != "miRNA_miRNA"]
TFT = rd("layer_TF_target.tsv")
GG = rd("layer_gene_gene.tsv")
MM = [r for r in rd("layer_miRNA_miRNA.tsv") if str(r.get("threshold_10kb", "TRUE")).upper() in ("TRUE", "1")]
STRING_ROWS = [i for i, r in enumerate(GG) if r["evidence"] == "STRING"]
assert len(STRING_ROWS) == 941


def build(flip=None, string=True, trrust_layer=True, gene_gene_file=True):
    """flip: boolean array over the 941 STRING rows (True = enter as target->source)."""
    edges, layer = {}, {}

    def add(a, b, lab):
        if a not in nodes or b not in nodes or a == b: return
        if (a, b) not in edges: edges[(a, b)] = 1; layer[(a, b)] = lab
    for r in CANON: edges[(r["source"], r["target"])] = 1; layer[(r["source"], r["target"])] = "analysed"
    if trrust_layer:
        for r in TFT: add(r["source"], r["target"], "TRRUST")
    if gene_gene_file:
        k = 0
        for i, r in enumerate(GG):
            if r["evidence"] == "STRING":
                f = flip is not None and flip[k]; k += 1
                if string: add(*((r["target"], r["source"]) if f else (r["source"], r["target"])), "STRING")
            else:
                add(r["source"], r["target"], "gene_gene_TRRUST")
    for r in MM:
        add(r["miRNA_1"], r["miRNA_2"], "10kb"); add(r["miRNA_2"], r["miRNA_1"], "10kb")
    nrec = 0
    for (a, b) in list(edges):
        if nodes.get(a) == "TF" and nodes.get(b) == "miRNA" and (b, a) in edges:
            edges.pop((b, a)); layer.pop((b, a)); nrec += 1
    assert nrec == 1223, nrec
    return set(edges), layer


def write(fn, arcs):
    with open(fn, "w") as fh:
        fh.write(f"{len(V)} {len(arcs)} 1\n" + " ".join(str(TMAP[nodes[v]]) for v in V) + "\n")
        for a, b in sorted((idx[a], idx[b]) for a, b in arcs): fh.write(f"{a} {b} 0\n")


def read_arcs(fn, fine=True):
    L = open(os.path.join(GR, fn)).read().split("\n")
    n, m = map(int, L[0].split()[:2]); start = 2 if fine else 1 + n
    return {tuple(map(int, l.split()[:2])) for l in L[start:start + m]}


os.makedirs(os.path.join(HERE, "work", "graphs"), exist_ok=True)
listed, lay = build()
assert {(idx[a], idx[b]) for a, b in listed} == read_arcs("census_noexemplar_fine.txt"), "listed graph differs"
nostr, _ = build(string=False)
assert {(idx[a], idx[b]) for a, b in nostr} == read_arcs("census_noexemplar_nostring_edgeclass.txt", fine=False), \
    "no-STRING graph differs"
nogg, _ = build(string=False, trrust_layer=False, gene_gene_file=False)
print(f"listed {len(listed)} arcs (STRING {sum(1 for e in listed if lay[e] == 'STRING')}); "
      f"no STRING {len(nostr)}; no gene-gene layer {len(nogg)}; checks passed")
G = os.path.join(HERE, "work", "graphs")
write(os.path.join(G, "listed.txt"), listed); write(os.path.join(G, "nostring.txt"), nostr)
write(os.path.join(G, "nogenegene.txt"), nogg)
with open(os.path.join(HERE, "work", "string_pairs.txt"), "w") as fh:
    fh.write(f"{len(STRING_ROWS)}\n")
    for i in STRING_ROWS: fh.write(f"{idx[GG[i]['source']]} {idx[GG[i]['target']]}\n")
sa = sorted((idx[a], idx[b]) for a, b in listed if lay[(a, b)] == "STRING")
with open(os.path.join(HERE, "work", "string_arcs_listed.txt"), "w") as fh:
    fh.write(f"{len(sa)}\n" + "".join(f"{a} {b}\n" for a, b in sa))
rows = []
for r in range(R):
    flip = np.random.default_rng(SEED0 + r).random(941) < 0.5
    arcs, ly = build(flip=flip)
    write(os.path.join(G, f"boot_{r:03d}.txt"), arcs)
    rows.append(dict(replicate=r, seed=SEED0 + r, n_flipped=int(flip.sum()), arcs=len(arcs),
                     string_arcs_added=sum(1 for e in arcs if ly[e] == "STRING")))
with open(os.path.join(HERE, "string_bootstrap_graphs.csv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"{R} bootstrap graphs written; STRING arcs added: min {min(r['string_arcs_added'] for r in rows)}, "
      f"max {max(r['string_arcs_added'] for r in rows)}")
