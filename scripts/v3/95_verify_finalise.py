#!/usr/bin/env python3
"""Part 5: re-run the (A) architecture checks so they are written to file, verify the
targeted-set knockouts, consolidate all verification parts, and record the corrections."""
import csv, json, collections, os
import numpy as np, networkx as nx
REV='/path/to/revision'; D=REV+'/data'; R3=REV+'/results/v3'
V=[]
def rec(sec,item,mine,stored,ok=None,note=''):
    if ok is None:
        try: ok=abs(float(mine)-float(stored))<1e-6
        except Exception: ok=str(mine)==str(stored)
    V.append(dict(section=sec,item=item,independent_value=mine,stored_value=stored,agree=bool(ok),note=note))
    print(f"[{'OK ' if ok else 'DIFF'}] {sec} | {item}: mine={mine} stored={stored} {note}",flush=True)

ntype={r['name']:r['type'] for r in csv.DictReader(open(D+'/canonical_nodes.tsv'),delimiter='\t')}
edges=[(r['source'],r['target'],int(r['sign'])) for r in csv.DictReader(open(D+'/canonical_edges.tsv'),delimiter='\t')]
NODES=sorted(ntype); idx={n:i for i,n in enumerate(NODES)}; N=len(NODES)
P=[(u,v) for u,v,_ in edges]; sgn={(u,v):s for u,v,s in edges}
G=nx.DiGraph(); G.add_nodes_from(NODES); G.add_edges_from(P)
succ={n:set(G.successors(n)) for n in NODES}
rec('A','nodes',G.number_of_nodes(),587); rec('A','directed edges',G.number_of_edges(),6859)
rec('A','duplicate directed edges',len(P)-len(set(P)),0)
rec('A','weakly connected components',nx.number_weakly_connected_components(G),1)
sccs=sorted(nx.strongly_connected_components(G),key=len,reverse=True)
rec('A','strongly connected components',len(sccs),224)
rec('A','largest SCC = bow-tie CORE',len(sccs[0]),362)
rec('A','singleton SCCs',sum(1 for c in sccs if len(c)==1),221)
core=sccs[0]; a=next(iter(core))
OUT=nx.descendants(G,a)-core; IN=nx.ancestors(G,a)-core
rec('A','bow-tie IN',len(IN),7); rec('A','bow-tie OUT',len(OUT),215)
rec('A','tubes + tendrils',N-len(core)-len(IN)-len(OUT),3)
rec('A','CORE composition (TF/miRNA/Gene)',
    f"{sum(1 for n in core if ntype[n]=='TF')}/{sum(1 for n in core if ntype[n]=='miRNA')}/{sum(1 for n in core if ntype[n]=='Gene')}",
    '152/210/0')
rec('A','OUT composition (Gene/TF/miRNA)',
    f"{sum(1 for n in OUT if ntype[n]=='Gene')}/{sum(1 for n in OUT if ntype[n]=='TF')}/{sum(1 for n in OUT if ntype[n]=='miRNA')}",
    '207/2/6')
rec('A','edges inside the core',sum(1 for u,v in P if u in core and v in core),2871)
mut=sum(1 for u,v in P if (v,u) in sgn)
rec('A','reciprocity',round(mut/len(P),4),0.3566)
rec('A','mutual TF<->miRNA dyads',mut//2,1223)
rec('A','flow hierarchy',round(nx.flow_hierarchy(G),4),0.5808)
rec('A','undirected transitivity',round(nx.transitivity(G.to_undirected()),4),0.035,
    ok=abs(nx.transitivity(G.to_undirected())-0.035)<5e-4)
lr={n:len(nx.descendants(G,n))/(N-1) for n in NODES}; mx=max(lr.values())
rec('A','global reaching centrality',round(sum(mx-v for v in lr.values())/(N-1),4),0.3688)
rec('A','out-degree-0 nodes (all genes)',
    f"{sum(1 for n in NODES if G.out_degree(n)==0)} ({sum(1 for n in NODES if G.out_degree(n)==0 and ntype[n]=='Gene')} genes)",
    '206 (206 genes)')
rec('A','in-degree-0 nodes',sum(1 for n in NODES if G.in_degree(n)==0),5)
comp=other=0
for Rn in NODES:
    tR=ntype[Rn]
    for M in succ[Rn]:
        tM=ntype[M]
        if not ((tR=='miRNA' and tM=='TF') or (tR=='TF' and tM=='miRNA')): continue
        ic=Rn in succ[M]
        for T in succ[Rn]&succ[M]:
            if T in (Rn,M) or ntype[T]=='miRNA': continue
            if ic: comp+=1
            else: other+=1
rec('A','3-node FFL cores (unique)',comp//2+other,1649)
rec('A','3-node FFL cores (raw directed)',comp+other,3083)
# regulator-only subnetwork (the actual hierarchy test)
reg=[n for n in NODES if ntype[n] in ('TF','miRNA')]
Gr=G.subgraph(reg)
rec('A','regulator subnetwork N / M',f'{Gr.number_of_nodes()} / {Gr.number_of_edges()}','380 / 2911')
sr=sorted(nx.strongly_connected_components(Gr),key=len,reverse=True)
rec('A','regulator subnetwork largest SCC',len(sr[0]),362)
rec('A','regulator subnetwork fraction in the giant SCC',round(len(sr[0])/Gr.number_of_nodes(),4),0.9526,
    ok=abs(len(sr[0])/Gr.number_of_nodes()-0.9526)<1e-3)
rec('A','regulator subnetwork flow hierarchy',round(nx.flow_hierarchy(Gr),4),0.0124,
    ok=abs(nx.flow_hierarchy(Gr)-0.012366884)<1e-4,
    note='KEY: among regulators the network is almost entirely recurrent, not hierarchical')
C=nx.condensation(Gr)
rec('A','regulator condensation DAG depth (layers)',nx.dag_longest_path_length(C)+1,3)
Cf=nx.condensation(G)
rec('A','full-network condensation DAG depth (layers)',nx.dag_longest_path_length(Cf)+1,4)

# ---- targeted set knockouts
LAM=0.15; COLS=['COL1A1','COL3A1']
A_s=np.zeros((N,N))
for (u,v),s in sgn.items(): A_s[idx[u],idx[v]]=s
def ffl_unique(alive):
    s2={n:succ[n]&alive for n in alive}; c=o=0
    for Rn in alive:
        tR=ntype[Rn]
        for M in s2[Rn]:
            tM=ntype[M]
            if not ((tR=='miRNA' and tM=='TF') or (tR=='TF' and tM=='miRNA')): continue
            ic=Rn in s2[M]
            for T in s2[Rn]&s2[M]:
                if T in (Rn,M) or ntype[T]=='miRNA': continue
                if ic: c+=1
                else: o+=1
    return c//2+o
ALL=set(NODES); FFL_WT=ffl_unique(ALL)
desc={n:nx.descendants(G,n) for n in NODES}
def cinf(alive):
    ai=sorted(idx[n] for n in alive); W=A_s[np.ix_(ai,ai)].copy()
    od=np.abs(W).sum(1); nz=od>0; W[nz]=W[nz]/od[nz,None]
    S=np.linalg.inv(np.eye(len(ai))-(1-LAM)*W); np.fill_diagonal(S,0.0)
    pos={NODES[j]:k for k,j in enumerate(ai)}
    return sum(float(np.abs(S[:,pos[c]]).sum()) for c in COLS if c in pos)
INF=cinf(ALL)
def dmg(ko):
    ko=set(ko); alive=ALL-ko
    dF=(FFL_WT-ffl_unique(alive))/FFL_WT
    g=G.subgraph(alive)
    rw=sum(len(desc[n]-ko) for n in alive); rk=sum(len(nx.descendants(g,n)) for n in alive)
    dR=(rw-rk)/rw if rw else 0.0
    dC=float('nan') if ko&set(COLS) else (INF-cinf(alive))/INF
    return dF,dR,dC,float(np.mean([abs(dF),abs(dR),0.0 if np.isnan(dC) else abs(dC)]))
sing={r['node']:float(r['composite_damage']) for r in csv.DictReader(open(R3+'/systems_single_knockout.csv'))}
bs=max(sing.values())
for r in csv.DictReader(open(R3+'/systems_targeted_set_knockouts.csv')):
    mem=r['members'].split(';')
    dF,dR,dC,c=dmg(mem)
    rec('D',f"targeted set: {r['set_name']}",round(c,6),float(r['composite_damage']),
        ok=abs(c-float(r['composite_damage']))<1e-5,
        note=f"ratio to best single {c/bs:.4f} vs stored {r['ratio_to_best_single']}")

with open(R3+'/systems_independent_verification_part1.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=['section','item','independent_value','stored_value','agree','note'])
    w.writeheader(); w.writerows(V)

# ---- consolidate
allrows=[]
for p in ['part1','part2','part3','part4']:
    fn=f'{R3}/systems_independent_verification_{p}.csv'
    if os.path.exists(fn):
        for r in csv.DictReader(open(fn)): r['part']=p; allrows.append(r)
with open(R3+'/systems_verification_summary.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=['part','section','item','independent_value','stored_value','agree','note'])
    w.writeheader(); w.writerows(allrows)
nok=sum(1 for r in allrows if r['agree']=='True')
print(f"\n===== CONSOLIDATED: {nok}/{len(allrows)} independently reproduced =====")
for r in allrows:
    if r['agree']!='True': print('  NOT REPRODUCED:',r['section'],'|',r['item'],'|',r['independent_value'],'vs',r['stored_value'])
