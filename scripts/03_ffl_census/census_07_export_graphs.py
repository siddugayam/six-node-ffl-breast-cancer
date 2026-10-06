#!/usr/bin/env python3
# Writes a census graph in the plain-text format that the C census programs read (counts of nodes, arcs and
# edge classes; node types; one line per arc with its class) and a .meta.json with node names and classes.
import sys, os, json
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_06_graph_definitions import make, reduce_graph, fine_class
from census_03_build_graph import NODE_TYPE

def export(name, path, classing="fine"):
    E = make(name)
    E2, comp = reduce_graph(E)
    nodes = sorted(NODE_TYPE)
    idx = {n: i for i, n in enumerate(nodes)}
    tmap = {"miRNA": 0, "TF": 1, "Gene": 2}
    if classing == "fine":
        classes = sorted({fine_class(*e) for e in E2})
        cof = lambda e: classes.index(fine_class(*e))
    else:
        classes = sorted(set(E2.values()))
        cof = lambda e: classes.index(E2[e])
    with open(path, "w") as fh:
        fh.write("%d %d %d\n" % (len(nodes), len(E2), len(classes)))
        fh.write(" ".join(str(tmap[NODE_TYPE[n]]) for n in nodes) + "\n")
        for e in E2:
            fh.write("%d %d %d\n" % (idx[e[0]], idx[e[1]], cof(e)))
    json.dump({"nodes": nodes, "classes": classes}, open(path + ".meta.json", "w"))
    print(name, classing, "-> ", path, "|V|", len(nodes), "|E'|", len(E2), "classes", classes)

if __name__ == "__main__":
    out = "/path/to/revision/results/v2"
    for nm in ("dep", "pub", "sym", "nomm"):
        export(nm, os.path.join(out, "graph_%s_fine.txt" % nm), "fine")
    export("pub", os.path.join(out, "graph_pub_coarse.txt"), "coarse")

def export_original(name, path):
    """Original (un-contracted) graph, fine classes -- input to the null models."""
    E = make(name)
    nodes = sorted(NODE_TYPE)
    idx = {n: i for i, n in enumerate(nodes)}
    tmap = {"miRNA": 0, "TF": 1, "Gene": 2}
    classes = sorted({fine_class(*e) for e in E})
    with open(path, "w") as fh:
        fh.write("%d %d %d\n" % (len(nodes), len(E), len(classes)))
        fh.write(" ".join(str(tmap[NODE_TYPE[n]]) for n in nodes) + "\n")
        for e in E:
            fh.write("%d %d %d\n" % (idx[e[0]], idx[e[1]], classes.index(fine_class(*e))))
    json.dump({"nodes": nodes, "classes": classes}, open(path + ".meta.json", "w"))
    print("ORIG", name, "->", path, "|E|", len(E), "classes", classes)
