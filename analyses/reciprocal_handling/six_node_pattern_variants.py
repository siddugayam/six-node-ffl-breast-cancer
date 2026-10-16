#!/usr/bin/env python3
"""The six-node composite pattern under the four treatments of reciprocal pairs.

Recounts the pattern with analyses/six_node_pattern/S1/s1_null6.c, unchanged (mode OBS, census graph without the
exemplar edges; its graph is uncontracted), and checks the instance list against the deposited one
(analyses/six_node_pattern/S1/obs/noexemplar_bhat6_composite_instances.tsv). Then, for each of its distinct node sets:
  B  topological-definition module of the contracted census graph       U  topological-definition module of the uncontracted census graph
  S  its six nodes lie in six distinct SCCs of the census graph (needed for an image on the condensation DAG), and
     whether that image is a module of the DAG
  M  valid under module-level condensation
The tests are those of check_small_n.py. Writes six_node_pattern_variants.csv.
usage: FFL_ROOT=<repository checkout> python3 six_node_pattern_variants.py"""
import csv, os, subprocess
import networkx as nx
from check_small_n import load, d234, mtest

ROOT = os.environ["FFL_ROOT"]
HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(ROOT, "analyses", "six_node_pattern", "S1")
W = os.path.join(HERE, "work"); os.makedirs(os.path.join(W, "bin"), exist_ok=True)

exe = os.path.join(W, "bin", "s1_null6")
subprocess.run(["gcc", "-O2", "-o", exe, os.path.join(S1, "s1_null6.c")], check=True)
pref = os.path.join(W, "s1obs")
obs = subprocess.run([exe, os.path.join(S1, "graph_noexemplar_null.txt"), os.path.join(S1, "graph_noexemplar_labels.txt"),
                      "OBS", "0", "20250908", "100", "1", "1", "0", pref], check=True, capture_output=True, text=True).stdout
hdr, val = [l.split("\t") for l in obs.strip().split("\n")]
o = dict(zip(hdr, val))
inst = open(pref + "_bhat6_composite_instances.tsv").read()
same = inst == open(os.path.join(S1, "obs", "noexemplar_bhat6_composite_instances.tsv")).read()
rows = [tuple(map(int, l.split("\t"))) for l in inst.strip().split("\n")]
sets = sorted({tuple(sorted(r)) for r in rows})
print("recount:", o["bhat6_comp_inst"], "instances,", o["bhat6_comp_sets"], "sets; instance list identical to deposit:", same)

names = open(os.path.join(S1, "node_names.txt")).read().split("\n")
EC = open(os.path.join(ROOT, "analyses", "census_and_motif_nulls", "graphs", "census_noexemplar_edgeclass.txt")).read().split("\n")
assert [l.split("\t")[0] for l in EC[1:588]] == names[:587], "node orders of the S1 and census graphs differ"
n, cout_ = load("census_noexemplar_fine.txt"); _, uout = load("census_noexemplar_uncontracted.txt")
G = nx.DiGraph(); G.add_nodes_from(range(n)); G.add_edges_from((a, b) for a in range(n) for b in uout[a])
comps = list(nx.strongly_connected_components(G)); cid = {v: i for i, c in enumerate(comps) for v in c}
CD = nx.condensation(G, comps); cd_out = [set(CD.successors(c)) for c in range(CD.number_of_nodes())]
cnt = dict(B=0, U=0, M=0, S_distinct=0, S_module=0)
firstB = {}
for S in sets:
    rb = d234(S, cout_); firstB[rb] = firstB.get(rb, 0) + 1
    cnt["B"] += rb == 0; cnt["U"] += d234(S, uout) == 0; cnt["M"] += mtest(S, uout) == 0
    img = {cid[v] for v in S}
    if len(img) == 6:
        cnt["S_distinct"] += 1; cnt["S_module"] += d234(img, cd_out) == 0
inst_sets = {}
for r in rows: inst_sets.setdefault(tuple(sorted(r)), 0); inst_sets[tuple(sorted(r))] += 1
inst_M = sum(c for s, c in inst_sets.items() if mtest(s, uout) == 0)
out = [dict(treatment="pattern recount (s1_null6 OBS, uncontracted census graph)", instances=int(o["bhat6_comp_inst"]),
            node_sets=int(o["bhat6_comp_sets"]), note=f"instance list identical to the deposited one: {same}"),
       dict(treatment="U no contraction", instances=int(o["bhat6_comp_inst"]), node_sets=int(o["bhat6_comp_sets"]),
            note="the pattern already reads both arcs of every reciprocal pair, so its graph is the uncontracted one"),
       dict(treatment="S graph-level SCC condensation: node sets with six distinct SCCs", instances="",
            node_sets=cnt["S_distinct"], note=f"of which topological-definition modules of the condensation DAG: {cnt['S_module']}"),
       dict(treatment="M module-level condensation: node sets that are valid modules", instances=inst_M,
            node_sets=cnt["M"], note=""),
       dict(treatment="B baseline census (contracted): node sets that are topological-definition modules", instances="",
            node_sets=cnt["B"], note="first failing condition: " + "; ".join(
                f"{ {0: 'none', 2: 'D2', 3: 'D3', 4: 'D4'}[k]} {v}" for k, v in sorted(firstB.items()))),
       dict(treatment="U: node sets that are topological-definition modules of the uncontracted graph", instances="", node_sets=cnt["U"],
            note="")]
with open(os.path.join(HERE, "six_node_pattern_variants.csv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for r in out: print(r)
