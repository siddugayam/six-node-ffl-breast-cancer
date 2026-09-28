#!/usr/bin/env python3
"""02f_crossvalidate.py -- independent cross-validation of scripts/ffl_enum.c.
Enumerates candidate modules from underlying triangles (condition (i) is
necessary, so every valid module contains one) rather than from (source,sink)
pairs, and tests them with the reference implementation in ffl_def.py."""
import sys, itertools, time
sys.path.insert(0,"/path/to/revision/scripts")
import ffl_graph as fg, ffl_def
nodes,ntype,arcs,meta=fg.build(verbose=False)
und=set(e for e,a in arcs.items() if a["undirected"])
A=set(arcs)
# underlying undirected adjacency
nb={n:set() for n in nodes}
for (u,v) in A: nb[u].add(v); nb[v].add(u)
idx={n:i for i,n in enumerate(nodes)}
# INDEPENDENT generator: every valid module must contain an (underlying) triangle,
# because condition (i) requires arcs R->M, R->T, M->T.
t0=time.time()
tris=[]
for u in nodes:
    for v in nb[u]:
        if idx[v]<=idx[u]: continue
        for w in nb[u]&nb[v]:
            if idx[w]<=idx[v]: continue
            tris.append((u,v,w))
print("underlying triangles:",len(tris),"  %.1fs"%(time.time()-t0))
# n=3
t0=time.time(); n3=set()
for T in tris:
    if ffl_def.test_set(T,A,und) is not None: n3.add(frozenset(T))
print("PY n=3 valid modules:",len(n3),"  %.1fs"%(time.time()-t0))
# n=4 : triangle + 1 node adjacent to it
t0=time.time(); n4=set(); tested=set()
for T in tris:
    cand=(nb[T[0]]|nb[T[1]]|nb[T[2]])-set(T)
    for x in cand:
        S=frozenset(T)|{x}
        if S in tested: continue
        tested.add(S)
        if ffl_def.test_set(sorted(S),A,und) is not None: n4.add(S)
print("PY n=4 candidate sets:",len(tested),"valid:",len(n4),"  %.1fs"%(time.time()-t0))
import pickle
pickle.dump((n3,n4),open("/path/to/scratch/py_n34.pkl","wb"))
