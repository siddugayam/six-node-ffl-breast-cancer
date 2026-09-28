#!/usr/bin/env python3
"""02c_candidate_space.py -- computes, exactly, the number of candidate vertex
sets the exhaustive enumeration must examine for each module size n (the figures
in section 8 of results/ffl_definition_check.md).  This is what establishes that
n=7 (3.67e13) and n=8 (4.09e15) are intractable while n<=6 is not."""
import sys, math
sys.path.insert(0,"/path/to/revision/scripts")
import ffl_graph as fg
nodes, ntype, arcs, meta = fg.build(verbose=False)
out, inn = fg.adjacency(nodes, arcs)
A=set(arcs)
def bfs(adj,s,L):
    d={s:0}; fr=[s]
    for st in range(1,L+1):
        nx=[]
        for u in fr:
            for v in adj[u]:
                if v not in d: d[v]=st; nx.append(v)
        fr=nx
        if not fr: break
    return d
def cnt(a,b,c,d,k,rmin):
    tot=0
    for ia in range(0,min(a,k)+1):
     for ib in range(0,min(b,k-ia)+1):
      for ic in range(0,min(c,k-ia-ib)+1):
        idd=k-ia-ib-ic
        if idd<0 or idd>d: continue
        if ia+ib<rmin or ia+ic<rmin: continue
        tot+=math.comb(a,ia)*math.comb(b,ib)*math.comb(c,ic)*math.comb(d,idd)
    return tot
for n in (4,5,6,7,8):
    L=n-1;k=n-2
    fwd={s:bfs(out,s,L) for s in nodes}
    bwd={t:bfs(inn,t,L) for t in nodes}
    tot=0;pairs=0
    for s in nodes:
        fs=fwd[s]
        for t in nodes:
            if t==s or t not in fs: continue
            if (t,s) in A and (s,t) not in A: continue
            bt=bwd[t]
            U=[v for v in fs if v!=s and v!=t and v in bt and fs[v]+bt[v]<=L]
            U=[v for v in U if not((v,s) in A and (s,v) not in A) and not((t,v) in A and (v,t) not in A)]
            if len(U)<k: continue
            aa=bb=cc=dd=0
            for v in U:
                o=(s,v) in A;i=(v,t) in A
                if o and i: aa+=1
                elif o: bb+=1
                elif i: cc+=1
                else: dd+=1
            rmin = 1 if (s,t) in A else 2
            tot+=cnt(aa,bb,cc,dd,k,rmin); pairs+=1
    print("n=%d pairs=%d  candidate_subsets=%.5e"%(n,pairs,tot))
