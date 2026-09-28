#!/usr/bin/env python3
"""Resolve the composite-FFL counting ambiguity definitively.

A composite FFL is (TF, miRNA, gene) where TF->miRNA AND miRNA->TF both exist, and both act on
the gene. The question is whether that is ONE motif instance or TWO. It is one: the reciprocal
pair is a single mutual regulatory relationship comprising two directed arcs, exactly the point
This script reports both numbers and the exact factor, so the manuscript can state one."""
import csv, collections, itertools
REV='/path/to/revision'

nodes={r['name']:r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t')}
E=set(); ET={}
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
    E.add((r['source'],r['target'])); ET[(r['source'],r['target'])]=r['edge_type']
out=collections.defaultdict(set)
for a,b in E: out[a].add(b)

recip = {(t,m) for (t,m) in E if nodes.get(t)=='TF' and nodes.get(m)=='miRNA' and (m,t) in E}
print(f"reciprocal TF<->miRNA pairs: {len(recip)}")

# --- count 1: one instance per unordered {TF,miRNA} pair per shared gene target
once=0; twice=0; targets_seen=set()
for (t,m) in recip:
    shared = {g for g in (out[t] & out[m]) if nodes.get(g) != 'miRNA' and g not in (t,m)}
    once += len(shared)
    twice += 2*len(shared)          # counting the pair from both arc orientations
    for g in shared: targets_seen.add((t,m,g))
print(f"composite cores counted ONCE per reciprocal pair : {once}")
print(f"composite cores counted TWICE (per arc)          : {twice}")
print(f"distinct (TF, miRNA, gene) triples               : {len(targets_seen)}")
print(f"ratio                                            : {twice/once:.2f}")

# --- non-composite classes for completeness
tf_ffl=0
for (t,m) in E:
    if nodes.get(t)=='TF' and nodes.get(m)=='miRNA' and (t,m) not in recip:
        tf_ffl += len({g for g in (out[t] & out[m]) if nodes.get(g)!='miRNA' and g not in (t,m)})
mir_ffl=0
for (m,t) in E:
    if nodes.get(m)=='miRNA' and nodes.get(t)=='TF' and (t,m) not in E:
        mir_ffl += len({g for g in (out[m] & out[t]) if nodes.get(g)!='miRNA' and g not in (m,t)})
print(f"\nTF-FFL (non-composite)   : {tf_ffl}")
print(f"miRNA-FFL (non-composite): {mir_ffl}")
print(f"TOTAL, composite counted once : {once+tf_ffl+mir_ffl}")
print(f"TOTAL, composite counted twice: {twice+tf_ffl+mir_ffl}")
print("\nRESOLUTION: report the ONCE convention. Each reciprocal TF<->miRNA relationship is a")
print("single composite motif comprising two directed arcs. The 2,868")
print("figure in one intermediate output is the twice convention and must not be used.")
