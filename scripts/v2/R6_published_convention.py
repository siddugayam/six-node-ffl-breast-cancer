#!/usr/bin/env python3
"""Reproduce the counting convention behind the published 1,434 / 206 / 9 / 1,649,
directly from data/canonical_edges.tsv."""
import csv, collections
REV='/path/to/revision'
nodes={r['name']:r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t')}
E=set(); out=collections.defaultdict(set)
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
    E.add((r['source'],r['target'])); out[r['source']].add(r['target'])
TFs=[v for v,t in nodes.items() if t=='TF']; MIRs=[v for v,t in nodes.items() if t=='miRNA']
nonmir=set(v for v,t in nodes.items() if t!='miRNA')
comp_pair=comp_arc=tffl=mirfl=0
for tf in TFs:
    for m in out[tf]:
        if nodes.get(m)!='miRNA': continue
        shared=(out[tf]&out[m])&nonmir
        if (m,tf) in E:                      # reciprocal -> composite
            comp_pair+=len(shared); comp_arc+=len(shared)
        else:
            tffl+=len(shared)
for m in MIRs:
    for tf in out[m]:
        if nodes.get(tf)!='TF': continue
        shared=(out[m]&out[tf])&nonmir
        if (tf,m) in E:
            comp_arc+=len(shared)            # the reverse arc of the same composite
        else:
            mirfl+=len(shared)
print("PUBLISHED CONVENTION on the deposited network (target restricted to non-miRNA nodes,")
print("no acyclicity/induced test; composite indexed once per reciprocal (TF,miRNA) pair):")
print(f"  Composite-FFL (once per pair) = {comp_pair}")
print(f"  Composite-FFL (once per arc)  = {comp_arc}")
print(f"  TF-FFL  (non-reciprocal)      = {tffl}")
print(f"  miRNA-FFL (non-reciprocal)    = {mirfl}")
print(f"  Any-3-node total (pair conv.) = {comp_pair+tffl+mirfl}")
print(f"  Any-3-node total (arc conv.)  = {comp_arc+tffl+mirfl}")
