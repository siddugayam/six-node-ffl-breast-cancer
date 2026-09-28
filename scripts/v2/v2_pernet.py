#!/usr/bin/env python3
import sys, os, csv
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v2_build_graph import deposited_edges, NODE_TYPE
from v2_core3 import census3, reduce_graph

DATA = "/path/to/revision/data"
rows = list(csv.DictReader(open(os.path.join(DATA, "canonical_edges.tsv")), delimiter="\t"))
nets = ["3-miR", "3-TF", "3-Comp", "4-TF", "5-TF", "6-TF"]
for net in nets:
    E = {}
    for r in rows:
        if net in r["in_networks"].split(";"):
            E[(r["source"], r["target"])] = r["edge_type"]
    tot, cls, ordG, ordR = census3(E)
    nodes = set()
    for (s, t) in E:
        nodes.add(s); nodes.add(t)
    print("%-7s |V|=%4d |E|=%5d   D1-D4_n3=%6d  ordered_G=%6d ordered_reduced=%6d" %
          (net, len(nodes), len(E), tot, ordG, ordR))
    print("        ", dict(cls))
