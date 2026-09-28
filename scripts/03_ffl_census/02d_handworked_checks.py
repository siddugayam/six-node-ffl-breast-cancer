#!/usr/bin/env python3
"""02d_handworked_checks.py -- the hand-worked sanity checks of the formal FFL
definition (section 3 of results/ffl_definition_check.md), including H5, which is
'cascade with a side branch' objection reproduced formally."""
import sys, itertools
sys.path.insert(0,"/path/to/revision/scripts")
import ffl_def

def mk(edges, und=()):
    arcs=set(); u2=set()
    for (a,b) in edges: arcs.add((a,b))
    for (a,b) in und:
        arcs.add((a,b)); arcs.add((b,a)); u2.add((a,b)); u2.add((b,a))
    return arcs,u2

print("=== HAND-WORKED SANITY CHECKS ===")
cases=[
 ("H1 classical TF-FFL  TF->miR, TF->G, miR->G",
   [("TF","miR"),("TF","G"),("miR","G")],[], ["TF","miR","G"], True),
 ("H2 classical miRNA-FFL  miR->TF, miR->G, TF->G",
   [("miR","TF"),("miR","G"),("TF","G")],[], ["miR","TF","G"], True),
 ("H3 composite FFL  TF<->miR, TF->G, miR->G",
   [("TF","miR"),("miR","TF"),("TF","G"),("miR","G")],[], ["TF","miR","G"], True),
 ("H4 simple 3-cascade (NOT an FFL)  TF->miR->G, no TF->G",
   [("TF","miR"),("miR","G")],[], ["TF","miR","G"], False),
 ("H5 AUTHORS' 4-NODE  TF->miR, TF->G1, miR->G1, G1->G2",
   [("TF","miR"),("TF","G1"),("miR","G1"),("G1","G2")],[], ["TF","miR","G1","G2"], False),
 ("H6 4-node parallel fan, no core  TF->m1,TF->m2,m1->G,m2->G",
   [("TF","m1"),("TF","m2"),("m1","G"),("m2","G")],[], ["TF","m1","m2","G"], False),
 ("H7 4-node fan + miRNA-miRNA co-transcription (m1--m2)",
   [("TF","m1"),("TF","m2"),("m1","G"),("m2","G")],[("m1","m2")], ["TF","m1","m2","G"], True),
 ("H8 4-node fan + TF->G direct arc",
   [("TF","m1"),("TF","m2"),("m1","G"),("m2","G"),("TF","G")],[], ["TF","m1","m2","G"], True),
 ("H9 pure chain of 4 (cascade)  TF->m->G1->G2",
   [("TF","m"),("m","G1"),("G1","G2")],[], ["TF","m","G1","G2"], False),
 ("H10 3 co-transcribed miRNAs (2 undirected arcs on one path) + shared target",
   [("m1","G"),("m2","G"),("m3","G")],[("m1","m2"),("m2","m3")], ["m1","m2","m3","G"], None),
 ("H11 6-node authors' architecture TF1->TF2, TF1->m1, m1--m2, m1->G1, TF2->G1, G1->G2",
   [("TF1","TF2"),("TF1","m1"),("m1","G1"),("TF2","G1"),("G1","G2")],[("m1","m2")],
   ["TF1","TF2","m1","m2","G1","G2"], None),
]
for label,e,u,S,expect in cases:
    arcs,und=mk(e,u)
    r=ffl_def.test_set(S,arcs,und)
    fails,_=ffl_def.failing_condition(S,arcs,und)
    got = r is not None
    verdict = "" if expect is None else ("  [OK]" if got==expect else "  [*** MISMATCH ***]")
    print("%-70s -> %s %s%s"%(label, "VALID" if got else "INVALID",
          ("src=%s sink=%s ndisj=%d"%(r['source'],r['sink'],r['ndisj'])) if r else "fails: "+"; ".join(fails), verdict))
