#!/usr/bin/env python3
# Three-node FFL count (conditions D1-D4) for each of the six networks of the original submission, selected by
# the in_networks column of data/canonical_edges.tsv; prints the counts.
import sys, os, csv
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_03_build_graph import deposited_edges, NODE_TYPE
from census_04_n3_graph_variants import census3, reduce_graph

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
