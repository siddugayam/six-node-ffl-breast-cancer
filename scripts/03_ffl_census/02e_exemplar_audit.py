#!/usr/bin/env python3
"""02e_exemplar_audit.py -- tests the exemplar circuits of the original submission
(SIF_files/4-TF.sif, 5-TF.sif, 6-TF.sif) against the formal definition:
edge support, whole-structure size, exhaustive test of every induced n-subset,
and the subsets matching the composition given in the Methods."""
import sys, itertools, csv
sys.path.insert(0,"/path/to/revision/scripts")
import ffl_graph as fg, ffl_def
SIF="/path/to/revision/analyses/original_submission_code/SIF_files"
nodes,ntype,arcs,meta=fg.build(verbose=False)
und=set(e for e,a in arcs.items() if a["undirected"])
A=set(arcs)
# permissive graph (STRING association edges added) for a second opinion
n2,t2,arcs2,_=fg.build(include_string=True,verbose=False)
A2=set(arcs2); und2=set(e for e,a in arcs2.items() if a["undirected"])

for n,fn in ((4,"4-TF.sif"),(5,"5-TF.sif"),(6,"6-TF.sif")):
    rows=[l.rstrip("\n").split("\t") for l in open(f"{SIF}/{fn}") if l.strip()]
    E=[(r[0],r[2]) for r in rows]
    V=sorted(set([a for a,b in E]+[b for a,b in E]))
    print("="*78)
    print("%s : %d deposited edges, %d DISTINCT NODES (the manuscript calls it a %d-node FFL)"%(fn,len(E),len(V),n))
    dup=len(E)-len(set(frozenset(e) for e in E))
    unsupported=[e for e in E if (e[0],e[1]) not in A and (e[1],e[0]) not in A]
    unsup2=[e for e in E if (e[0],e[1]) not in A2 and (e[1],e[0]) not in A2]
    print("  duplicate undirected edges: %d ; edges with NO support in the canonical graph: %d/%d"
          %(dup,len(unsupported),len(E)))
    print("  ... still unsupported when STRING associations are allowed: %d"%len(unsup2))
    for e in unsupported: print("        UNSUPPORTED:",e[0],"--",e[1])
    missing=[v for v in V if v not in ntype]
    if missing: print("  nodes absent from the node universe:",missing)
    # whole deposited structure as a module of its own size
    r=ffl_def.test_set(V,A,und) if len(V)<=8 else None
    if len(V)<=8:
        f,_=ffl_def.failing_condition(V,A,und)
        print("  whole deposited structure as one module: %s %s"%("VALID" if r else "INVALID", "" if r else f))
    else:
        print("  whole deposited structure has %d nodes -> it is not an %d-node module at all"%(len(V),n))
    # all n-subsets of the exemplar node set
    ok=[]; fails={}
    for S in itertools.combinations(V,n):
        res=ffl_def.test_set(list(S),A,und)
        if res: ok.append((S,res))
        else:
            fc,_=ffl_def.failing_condition(list(S),A,und)
            for c in fc: fails[c]=fails.get(c,0)+1
    print("  induced %d-node subsets of the exemplar: %d ; satisfying the definition: %d"
          %(n,len(list(itertools.combinations(V,n))),len(ok)))
    print("  failure reasons (subsets may fail several):",dict(sorted(fails.items())))
    # composition of the valid ones
    from collections import Counter
    comp=Counter()
    for S,res in ok:
        comp[tuple(sorted(Counter(ntype[x] for x in S).items()))]+=1
    for k,v in comp.most_common(): print("     valid composition",dict(k),"->",v)
    for S,res in ok[:6]:
        print("     e.g.",",".join(S),"| source=%s sink=%s ndisj=%d"%(res['source'],res['sink'],res['ndisj']))
    # the architecture given in the Methods
    if n==4: want=[("TF",1),("miRNA",1),("Gene",2)]
    elif n==5: want=[("TF",1),("miRNA",2),("Gene",2)]
    else: want=[("TF",2),("miRNA",2),("Gene",2)]
    tgt=dict(want); nspec=0; nspec_ok=0
    for S in itertools.combinations(V,n):
        c=Counter(ntype[x] for x in S)
        if all(c.get(k,0)==v for k,v in tgt.items()):
            nspec+=1
            if ffl_def.test_set(list(S),A,und): nspec_ok+=1
    print("  subsets with the EXACT composition the manuscript specifies %s: %d, of which valid: %d"
          %(dict(tgt),nspec,nspec_ok))
