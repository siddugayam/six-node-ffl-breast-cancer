#!/usr/bin/env python3
"""
02h_exemplar_permissive.py -- re-tests the authors' deposited exemplars on the
PERMISSIVE graph, in which the 941 undirected STRING associations are admitted as
orientable edges (so COL1A1--COL3A1, which has no directed evidence anywhere, does
become available).  Also tests, directly, the higher-order wiring the manuscript
states in Methods (MS.md lines 217-226).

Result recorded in results/ffl_definition_check.md sections 4 and 5.
"""
import sys, itertools, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ffl_graph as fg, ffl_def

SIF = "/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/SIF_files"
n2, ntype, arcs2, _ = fg.build(include_string=True, verbose=False)
A2 = set(arcs2); und2 = set(e for e, a in arcs2.items() if a["undirected"])
print("PERMISSIVE GRAPH: canonical arcs + undirected STRING associations. arcs=%d" % len(A2))

for n, fn, tgt in ((4, "4-TF.sif", {'TF': 1, 'miRNA': 1, 'Gene': 2}),
                   (5, "5-TF.sif", {'TF': 1, 'miRNA': 2, 'Gene': 2}),
                   (6, "6-TF.sif", {'TF': 2, 'miRNA': 2, 'Gene': 2})):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(SIF, fn)) if l.strip()]
    V = sorted(set([r[0] for r in rows] + [r[2] for r in rows]))
    allok = specn = specok = 0
    fails = Counter(); shown = []
    for S in itertools.combinations(V, n):
        r = ffl_def.test_set(list(S), A2, und2)
        if r: allok += 1
        c = Counter(ntype[x] for x in S)
        if all(c.get(k, 0) == v for k, v in tgt.items()):
            specn += 1
            if r:
                specok += 1; shown.append((S, r))
            else:
                fc, _ = ffl_def.failing_condition(list(S), A2, und2)
                for x in fc: fails[x] += 1
    print("%s n=%d : %d/%d induced subsets valid; manuscript composition %s : %d/%d valid"
          % (fn, n, allok, len(list(itertools.combinations(V, n))), tgt, specok, specn))
    print("    failure reasons for manuscript-composition subsets:", dict(fails))
    for S, r in shown[:8]:
        print("      ", ",".join(S), "-> order:", " -> ".join(r['order']), "ndisj=%d" % r['ndisj'])

print("\nDIRECT TEST of the wiring stated in Methods (MS.md lines 217-226):")
def mk(edges, und=()):
    A = set(edges); U = set()
    for a, b in und:
        A.add((a, b)); A.add((b, a)); U.add((a, b)); U.add((b, a))
    return A, U
tests = [
 ("4-node  TF->miRNA, TF->Gene1, miRNA->Gene1, Gene1->Gene2",
  [("TF","miR"),("TF","G1"),("miR","G1"),("G1","G2")], [], ["TF","miR","G1","G2"]),
 ("5-node  + miRNA1--miRNA2 co-transcription, both -> Gene1",
  [("TF","m1"),("TF","m2"),("TF","G1"),("m1","G1"),("m2","G1"),("G1","G2")], [("m1","m2")],
  ["TF","m1","m2","G1","G2"]),
 ("6-node  + TF1->TF2 (miRNA-miRNA-TF-TF-gene-gene, Bhat et al. 2024 architecture)",
  [("TF1","TF2"),("TF1","m1"),("TF1","m2"),("TF2","G1"),("m1","G1"),("m2","G1"),("G1","G2")],
  [("m1","m2")], ["TF1","TF2","m1","m2","G1","G2"]),
]
for lab, e, u, S in tests:
    A, U = mk(e, u); r = ffl_def.test_set(S, A, U); f, _ = ffl_def.failing_condition(S, A, U)
    print("  %-72s -> %s" % (lab, ("VALID src=%s sink=%s ndisj=%d" % (r['source'], r['sink'], r['ndisj']))
                             if r else "INVALID: " + "; ".join(f)))
