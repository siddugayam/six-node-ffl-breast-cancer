#!/usr/bin/env python3
"""n=3 census with the undirected STRING tier excluded — the number the Methods sentence promises."""
import sys, importlib.util, collections, csv
REV='/path/to/revision'
spec=importlib.util.spec_from_file_location('cen', f'{REV}/scripts/03_ffl_census.py')
cen=importlib.util.module_from_spec(spec); sys.argv=['x','3']
try: spec.loader.exec_module(cen)
except SystemExit: pass
nodes,edges,nrecip=cen.build_graph(True)
string_pairs=set()
for r in csv.DictReader(open(f'{REV}/data/layer_gene_gene.tsv'),delimiter='\t'):
    if r['evidence']=='STRING' and r['source'] in nodes and r['target'] in nodes and r['source']!=r['target']:
        string_pairs.add((r['source'],r['target']))
dep={(r['source'],r['target']) for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t')}
drop={e for e in edges if edges[e]=='gene_gene' and e in string_pairs and e not in dep}
for e in drop: del edges[e]
print(f"dropped {len(drop)} STRING-orientation arcs -> {len(edges)} arcs remain")
out=collections.defaultdict(set); inn=collections.defaultdict(set); und=collections.defaultdict(set)
for (a,b) in edges: out[a].add(b); inn[b].add(a); und[a].add(b); und[b].add(a)
V=sorted(nodes); idx={v:i for i,v in enumerate(V)}
tot=0; sub=0; comp=collections.Counter()
for v in V:
    for u in [x for x in und[v] if idx[x]>idx[v]]:
        for w in set(x for x in und[v]|und[u] if idx[x]>idx[v] and x!=u):
            if idx[w]<idx[u] and w in und[v]: continue
            S={v,u,w}; sub+=1
            os_={x:(out[x]&S) for x in S}; is_={x:(inn[x]&S) for x in S}
            if not cen.is_ffl(S,os_,is_): continue
            tot+=1
            src=[x for x in S if not is_[x]][0]; snk=[x for x in S if not os_[x]][0]
            med=[x for x in S if x not in (src,snk)][0]
            comp[(nodes[src],nodes[med])]+=1
print(f"connected 3-subgraphs = {sub}")
print(f"n=3 FFL modules (STRING excluded) = {tot}")
for k,v in comp.most_common(): print(f"   source {k[0]:6s} -> mediator {k[1]:6s} : {v}")
