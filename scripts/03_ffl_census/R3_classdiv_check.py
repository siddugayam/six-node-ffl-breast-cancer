#!/usr/bin/env python3
"""Test the specific 5-node module reported to carry SEVEN edge classes, against the
paper's own _cls taxonomy and its own is_ffl predicate."""
import sys, importlib.util, collections, itertools
REV='/path/to/revision'
spec=importlib.util.spec_from_file_location('cen', f'{REV}/scripts/03_ffl_census/03_ffl_census.py')
cen=importlib.util.module_from_spec(spec); sys.argv=['x','3']
try: spec.loader.exec_module(cen)
except SystemExit: pass
nodes,edges,nrecip=cen.build_graph(True)

def cls(a,b):
    t=edges[(a,b)]
    if t=='TF_target': return 'TF->TF' if nodes[b]=='TF' else 'TF->Gene'
    if t=='miRNA_target': return 'miRNA->TF' if nodes[b]=='TF' else 'miRNA->Gene'
    return {'TF_miRNA':'TF->miRNA','gene_gene':'Gene->Gene','miRNA_miRNA':'miRNA-miRNA'}.get(t,t)

print("distinct classes possible under the paper's own _cls map:",
      sorted(set(cls(a,b) for (a,b) in edges)), "->", len(set(cls(a,b) for (a,b) in edges)))

cand=['EGR1','FOSB','JUND','hsa-miR-455','hsa-miR-29b']
print("\ncandidate 5-node module:",cand)
missing=[x for x in cand if x not in nodes]
print(" missing nodes:",missing)
if not missing:
    S=set(cand)
    out=collections.defaultdict(set); inn=collections.defaultdict(set)
    for (a,b) in edges: out[a].add(b); inn[b].add(a)
    os_={v:(out[v]&S) for v in S}; is_={v:(inn[v]&S) for v in S}
    arcs=[(a,b) for a in S for b in os_[a]]
    print(" induced arcs:")
    for a,b in arcs: print(f"   {a:14s} -> {b:14s}  [{edges[(a,b)]:12s}] class={cls(a,b)}")
    ec=set(cls(a,b) for a,b in arcs)
    print(" distinct edge classes =",len(ec),sorted(ec))
    r=cen.is_ffl(S,os_,is_)
    print(" passes D1-D4 (paper's own is_ffl)?",r)
