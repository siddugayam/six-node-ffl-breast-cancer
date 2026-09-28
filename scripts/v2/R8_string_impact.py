#!/usr/bin/env python3
"""How much of the n=3 census depends on STRING-orientation arcs?"""
import sys, importlib.util, collections, csv
REV='/path/to/revision'
spec=importlib.util.spec_from_file_location('cen', f'{REV}/scripts/03_ffl_census.py')
cen=importlib.util.module_from_spec(spec); sys.argv=['x','3']
try: spec.loader.exec_module(cen)
except SystemExit: pass
nodes,edges,_=cen.build_graph(True)
# which arcs came from the STRING rows of layer_gene_gene?
string_pairs=set(); trrust_pairs=set()
for r in csv.DictReader(open(f'{REV}/data/layer_gene_gene.tsv'),delimiter='\t'):
    a,b=r['source'],r['target']
    if a not in nodes or b not in nodes or a==b: continue
    (string_pairs if r['evidence']=='STRING' else trrust_pairs).add((a,b))
dep=set()
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
    dep.add((r['source'],r['target']))
string_arcs={e for e in edges if edges[e]=='gene_gene' and e in string_pairs and e not in dep}
gg_arcs={e for e in edges if edges[e]=='gene_gene'}
print(f"gene_gene arcs in the augmented graph: {len(gg_arcs)}")
print(f"  of which STRING-orientation arcs absent from the deposited network: {len(string_arcs)}")
print(f"  Gene->TF arcs (direction is an artefact of row order): "
      f"{sum(1 for (a,b) in gg_arcs if nodes[a]!='TF' and nodes[b]=='TF')}")
out=collections.defaultdict(set); inn=collections.defaultdict(set); und=collections.defaultdict(set)
for (a,b) in edges: out[a].add(b); inn[b].add(a); und[a].add(b); und[b].add(a)
V=sorted(nodes); idx={v:i for i,v in enumerate(V)}
tot=0; with_gg=0; with_string=0
for v in V:
    for u in [x for x in und[v] if idx[x]>idx[v]]:
        for w in set(x for x in und[v]|und[u] if idx[x]>idx[v] and x!=u):
            if idx[w]<idx[u] and w in und[v]: continue
            S={v,u,w}
            os_={x:(out[x]&S) for x in S}; is_={x:(inn[x]&S) for x in S}
            if not cen.is_ffl(S,os_,is_): continue
            tot+=1
            arcs=[(a,b) for a in S for b in os_[a]]
            if any(edges[e]=='gene_gene' for e in arcs): with_gg+=1
            if any(e in string_arcs for e in arcs): with_string+=1
print(f"\nn=3 FFL modules total                        : {tot}")
print(f"  containing >=1 gene_gene-layer arc         : {with_gg} ({100*with_gg/tot:.1f}%)")
print(f"  containing >=1 STRING-orientation arc      : {with_string} ({100*with_string/tot:.1f}%)")
