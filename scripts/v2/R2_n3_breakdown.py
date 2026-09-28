#!/usr/bin/env python3
"""Break the n=3 FFL census down by node-type composition, on both graphs."""
import sys, importlib.util, collections
REV='/path/to/revision'
spec=importlib.util.spec_from_file_location('cen', f'{REV}/scripts/03_ffl_census.py')
cen=importlib.util.module_from_spec(spec); sys.argv=['x','3']
try: spec.loader.exec_module(cen)
except SystemExit: pass

def run(tag, extra):
    nodes,edges,nrecip=cen.build_graph(extra_layers=extra)
    out=collections.defaultdict(set); inn=collections.defaultdict(set); und=collections.defaultdict(set)
    for (a,b) in edges: out[a].add(b); inn[b].add(a); und[a].add(b); und[b].add(a)
    V=sorted(nodes); idx={v:i for i,v in enumerate(V)}
    cnt=collections.Counter(); arch=collections.Counter()
    tot=0
    for v in V:
        for u in [x for x in und[v] if idx[x]>idx[v]]:
            for w in set(x for x in und[v]|und[u] if idx[x]>idx[v] and x!=u):
                if idx[w]<idx[u] and w in und[v]: continue
                S={v,u,w}
                os_={x:(out[x]&S) for x in S}; is_={x:(inn[x]&S) for x in S}
                if not cen.is_ffl(S,os_,is_): continue
                tot+=1
                src=[x for x in S if not is_[x]][0]; snk=[x for x in S if not os_[x]][0]
                med=[x for x in S if x not in (src,snk)][0]
                arch[(nodes[src],nodes[med],nodes[snk])]+=1
                R,M=nodes[src],nodes[med]
                if R=='TF' and M=='miRNA': cnt['TF->miRNA-mediated (TF-FFL/composite)']+=1
                elif R=='miRNA' and M=='TF': cnt['miRNA->TF-mediated (miRNA-FFL)']+=1
                elif R=='TF' and M=='TF': cnt['TF->TF-mediated']+=1
                else: cnt[f'other: {R}->{M}']+=1
    print(f"--- {tag}: total n=3 FFLs = {tot}")
    for k,v in cnt.most_common(): print(f"    {k:42s} {v}")
    print("    top source->mediator->sink type triples:")
    for k,v in arch.most_common(12): print(f"      {k}: {v}")
    return nodes,edges

na,ea=run('AUGMENTED',True)
print()
nd,ed=run('DEPOSITED',False)
# how many of the augmented TF->miRNA-mediated cores are composite (reciprocal pair contracted)?
