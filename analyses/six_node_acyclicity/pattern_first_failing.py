#!/usr/bin/env python3
"""Part (a): the six-node composite pattern against the topological definition (paper numbering: D1 induced and weakly connected; D2
acyclic; D3 one source and one sink, every vertex on a source-sink path; D4 two internally vertex-disjoint source-sink
paths).

The pattern instances are listed by the compiled six_node_d2prime_nulls (mode OBS, which runs bhat6() of s1_null6.c
unchanged) and checked against the deposited list analyses/six_node_pattern/S1/obs/noexemplar_bhat6_composite_instances.tsv.
Each of the distinct node sets is tested on the census graph (contracted, as the census uses:
analyses/census_and_motif_nulls/graphs/census_noexemplar_fine.txt). For the sets that fail D2, the directed cycles are
classified (2-cycles by the node types of the pair; cycles of length >= 3 that remain once the 2-cycle arcs are
removed), and D3-D4 are tested with D2 waived, on the cyclic module as it is and with its strongly connected components
condensed. Writes pattern_first_failing_condition.csv.
usage: FFL_ROOT=<repository checkout> python3 pattern_first_failing.py (after run_d2prime_nulls.sh compiled the program)"""
import collections, csv, os, subprocess
import networkx as nx

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(ROOT, "analyses", "six_node_pattern", "S1")
GR = os.path.join(ROOT, "analyses", "census_and_motif_nulls", "graphs")
os.makedirs(os.path.join(HERE, "work"), exist_ok=True)
pref = os.path.join(HERE, "work", "obs")
subprocess.run([os.path.join(HERE, "work", "bin", "d2p"), os.path.join(S1, "graph_noexemplar_null.txt"),
                os.path.join(S1, "graph_noexemplar_labels.txt"), "OBS", "0", "20250908", "100", "1", "1", "0", pref],
               check=True, capture_output=True)
inst = open(pref + "_bhat6_composite_instances.tsv").read()
assert inst == open(os.path.join(S1, "obs", "noexemplar_bhat6_composite_instances.tsv")).read(), "instance list differs"
rows = [tuple(map(int, l.split("\t"))) for l in inst.strip().split("\n")]
sets = sorted({tuple(sorted(r)) for r in rows})
names = open(os.path.join(S1, "node_names.txt")).read().split("\n")
L = open(os.path.join(GR, "census_noexemplar_edgeclass.txt")).read().split("\n")
assert [l.split("\t")[0] for l in L[1:588]] == names[:587], "node orders differ"
F = open(os.path.join(GR, "census_noexemplar_fine.txt")).read().split("\n")
nv, ne, _ = map(int, F[0].split()); ty = list(map(int, F[1].split()))
out = collections.defaultdict(set)
for l in F[2:2 + ne]:
    a, b, _ = map(int, l.split()); out[a].add(b)
TN = {0: "miRNA", 1: "TF", 2: "Gene"}
print(f"pattern: {len(rows)} instances, {len(sets)} node sets; instance list identical to the deposit")


def d34(H):
    src = [v for v in H if H.in_degree(v) == 0]; snk = [v for v in H if H.out_degree(v) == 0]
    if len(src) != 1 or len(snk) != 1 or src[0] == snk[0]: return "D3"
    s, t = src[0], snk[0]
    if (nx.descendants(H, s) | {s}) != set(H) or (nx.ancestors(H, t) | {t}) != set(H): return "D3"
    if H.has_edge(s, t):
        H2 = H.copy(); H2.remove_edge(s, t)
        return "pass" if nx.has_path(H2, s, t) else "D4"
    return "pass" if nx.node_connectivity(H, s, t) >= 2 else "D4"


res = collections.Counter()
for S in sets:
    H = nx.DiGraph(); H.add_nodes_from(S); H.add_edges_from((u, v) for u in S for v in out[u] if v in S)
    if not nx.is_weakly_connected(H): res[("D1", "", "", "")] += 1; continue
    if nx.is_directed_acyclic_graph(H):
        r = d34(H); res[(r if r != "pass" else "none (passes the topological definition)", "", "", "")] += 1; continue
    two = sorted({"-".join(sorted((TN[ty[u]], TN[ty[v]]))) for u, v in H.edges() if H.has_edge(v, u)})
    H3 = H.copy(); H3.remove_edges_from([(u, v) for u, v in H.edges() if H.has_edge(v, u)])
    longer = not nx.is_directed_acyclic_graph(H3)
    cls = ("2-cycles: " + ", ".join(two) if two else "no 2-cycle") + ("; a cycle of length >= 3" if longer else "")
    lit = d34(H)
    C = nx.condensation(H)
    con = d34(C) if C.number_of_nodes() > 1 else "D3"
    res[("D2", cls, lit, con)] += 1
with open(os.path.join(HERE, "pattern_first_failing_condition.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["first_failing_condition", "cycle_class", "with_D2_waived_literal", "with_D2_waived_condensed", "node_sets"])
    for k, v in sorted(res.items(), key=lambda kv: -kv[1]): w.writerow(list(k) + [v])
for k, v in sorted(res.items(), key=lambda kv: -kv[1]): print(k, v)
