#!/usr/bin/env python3
"""Reconciliation pass: rebuild the census graph exactly as 03_ffl_census.py does,
report its size, and re-enumerate n=3 FFLs exhaustively under D1-D4 on both the
augmented graph and the deposited-only graph."""
import sys, os, importlib.util, collections, itertools
REV='/path/to/revision'
spec=importlib.util.spec_from_file_location('cen', f'{REV}/scripts/03_ffl_census.py')
cen=importlib.util.module_from_spec(spec)
sys.argv=['x','3']
try:
    spec.loader.exec_module(cen)
except SystemExit:
    pass
except Exception as e:
    print("module exec note:", type(e).__name__, e)

def report(tag, extra):
    nodes,edges,nrecip=cen.build_graph(extra_layers=extra)
    print(f"[{tag}] nodes={len(nodes)} edges_after_contraction={len(edges)} reciprocal_contracted={nrecip}")
    cls=collections.Counter()
    for (a,b),et in edges.items(): cls[et]+=1
    print(f"[{tag}] edge_type counts: {dict(cls)}")
    tcls=collections.Counter((nodes[a],nodes[b]) for (a,b) in edges)
    print(f"[{tag}] (source_type,target_type) counts: {dict(tcls)}")
    out=collections.defaultdict(set); inn=collections.defaultdict(set)
    for (a,b) in edges: out[a].add(b); inn[b].add(a)
    und=collections.defaultdict(set)
    for (a,b) in edges: und[a].add(b); und[b].add(a)
    V=sorted(nodes); idx={v:i for i,v in enumerate(V)}
    n_conn=0; n_ffl=0
    ffls=[]
    # exhaustive connected 3-subgraph enumeration (min-label ESU)
    for v in V:
        ext=[u for u in und[v] if idx[u]>idx[v]]
        for i,u in enumerate(ext):
            # extend with neighbours of {v,u} greater than v, not already used
            cand=set(x for x in und[v]|und[u] if idx[x]>idx[v] and x!=u)
            for w in cand:
                if idx[w]<idx[u] and w in und[v]:
                    continue  # avoid double count: standard ESU exclusive neighbourhood
                S={v,u,w}
                n_conn+=1
                os_={x:(out[x]&S) for x in S}; is_={x:(inn[x]&S) for x in S}
                r=cen.is_ffl(S,os_,is_)
                if r: n_ffl+=1; ffls.append(tuple(sorted(S)))
    print(f"[{tag}] connected 3-subgraphs enumerated={n_conn}  n=3 FFLs={n_ffl}  distinct={len(set(ffls))}")
    return nodes,edges,ffls

print("=== AUGMENTED (extra_layers=True) — the graph 03_ffl_census.py actually uses ===")
na,ea,fa=report('augmented',True)
print()
print("=== DEPOSITED ONLY (extra_layers=False) ===")
nd,ed,fd=report('deposited',False)
