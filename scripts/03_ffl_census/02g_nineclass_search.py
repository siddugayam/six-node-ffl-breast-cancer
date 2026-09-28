#!/usr/bin/env python3
"""02g_nineclass_search.py -- exact targeted search for modules realising ALL NINE
interaction classes present in G.  Complete for the property tested: any such
module must contain a gene->gene, a gene->TF and a gene->miRNA arc (15/21/16 arcs
in the whole network) plus a TF-TF arc and a miRNA-miRNA edge, so the search can
be driven from those rare arcs instead of enumerating all n-node subsets."""
import sys, pickle, itertools
sys.path.insert(0,"/path/to/revision/scripts")
import ffl_graph as fg, ffl_def
nodes,ntype,arcs,_=fg.build(verbose=False)
A=set(arcs); und=set(e for e,a in arcs.items() if a["undirected"])
ICL={("TF","TF"):"TF-TF",("TF","Gene"):"TF-gene",("TF","miRNA"):"TF-miRNA",
     ("miRNA","Gene"):"miRNA-gene",("miRNA","TF"):"miRNA-TF",("miRNA","miRNA"):"miRNA-miRNA",
     ("Gene","TF"):"gene-TF",("Gene","miRNA"):"gene-miRNA",("Gene","Gene"):"gene-gene"}
def classes(S): return set(ICL[(ntype[u],ntype[v])] for u in S for v in S if u!=v and (u,v) in A)
gg=[(u,v) for (u,v) in A if ntype[u]=="Gene" and ntype[v]=="Gene"]
gT=[(u,v) for (u,v) in A if ntype[u]=="Gene" and ntype[v]=="TF"]
gM=[(u,v) for (u,v) in A if ntype[u]=="Gene" and ntype[v]=="miRNA"]
tfnb={};mnb={}
for (u,v) in A:
    if ntype[u]=="TF" and ntype[v]=="TF": tfnb.setdefault(u,set()).add(v); tfnb.setdefault(v,set()).add(u)
    if ntype[u]=="miRNA" and ntype[v]=="miRNA": mnb.setdefault(u,set()).add(v); mnb.setdefault(v,set()).add(u)
print("rare arcs: gene-gene=%d gene-TF=%d gene-miRNA=%d"%(len(gg),len(gT),len(gM)),flush=True)
# --- exhaustive over sets whose gene part is exactly the endpoints of a gg arc plus
#     (optionally) the sources of the gT / gM arcs
for nsize in (6,7,8):
    seen=set(); hits=[]
    for a1 in gg:
        for a2 in gT:
            for a3 in gM:
                base=frozenset(a1)|frozenset(a2)|frozenset(a3)
                if len(base)>nsize: continue
                for T2 in tfnb.get(a2[1],()):
                    for M2 in mnb.get(a3[1],()):
                        S=base|{T2,M2}
                        if len(S)!=nsize or S in seen: continue
                        seen.add(S)
                        if len(classes(S))<9: continue
                        r=ffl_def.test_set(sorted(S),A,und)
                        if r: hits.append((sorted(S),r))
    print("  n=%d exact-size candidate sets built from the five constraining arc classes: %d ; "
          "carrying all NINE classes AND valid: %d"%(nsize,len(seen),len(hits)),flush=True)
    for S,r in hits[:5]:
        print("      ",",".join(S)," src=%s sink=%s ndisj=%d"%(r['source'],r['sink'],r['ndisj']),flush=True)
# --- second route: extend the 8-class 6-node modules by one node supplying gene-miRNA
v6=pickle.load(open("/path/to/scratch/classcover6.pkl","rb"))
eight=[S for S,r in v6 if len(classes(S))==8]
print("8-class 6-node modules to extend: %d"%len(eight),flush=True)
cand=set(); 
for S in eight:
    genes=[x for x in S if ntype[x]=="Gene"]; mirs=[x for x in S if ntype[x]=="miRNA"]
    add=set()
    for g in genes:
        for (u,v) in A:
            if u==g and ntype[v]=="miRNA": add.add(v)
    for m in mirs:
        for (u,v) in A:
            if v==m and ntype[u]=="Gene": add.add(u)
    for x in add:
        if x in S: continue
        cand.add(frozenset(S)|{x})
print("7-node candidate supersets that could supply gene-miRNA: %d"%len(cand),flush=True)
hits7=[]
for S in cand:
    if len(classes(S))<9: continue
    r=ffl_def.test_set(sorted(S),A,und)
    if r: hits7.append((sorted(S),r))
print("  of which carry all NINE classes AND satisfy the definition: %d"%len(hits7),flush=True)
for S,r in hits7[:5]:
    print("      ",",".join(S)," src=%s sink=%s ndisj=%d"%(r['source'],r['sink'],r['ndisj']),flush=True)
