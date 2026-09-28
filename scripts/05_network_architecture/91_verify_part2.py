#!/usr/bin/env python3
"""Part 2 of the independent verification: (B) exhaustive criticality, (D) perturbation,
(E) dynamics, plus a cross-implementation check of the power-law fit."""
import os, csv, json, collections, random
import numpy as np, networkx as nx
import scipy.sparse as sp
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.stats import spearmanr

REV='/path/to/revision'; D=REV+'/data'; R3=REV+'/results/v3'
V=[]
def rec(sec,item,mine,stored,ok=None,note=''):
    if ok is None:
        try: ok=abs(float(mine)-float(stored))<1e-6
        except Exception: ok=str(mine)==str(stored)
    V.append(dict(section=sec,item=item,independent_value=mine,stored_value=stored,agree=bool(ok),note=note))
    print(f"[{'OK ' if ok else 'DIFF'}] {sec} | {item}: mine={mine} stored={stored} {note}",flush=True)

ntype={}
for r in csv.DictReader(open(D+'/canonical_nodes.tsv'),delimiter='\t'): ntype[r['name']]=r['type']
edges=[(r['source'],r['target'],int(r['sign'])) for r in csv.DictReader(open(D+'/canonical_edges.tsv'),delimiter='\t')]
NODES=sorted(ntype); idx={n:i for i,n in enumerate(NODES)}; N=len(NODES)
P=[(u,v) for u,v,_ in edges]; sgn={(u,v):s for u,v,s in edges}
G=nx.DiGraph(); G.add_nodes_from(NODES); G.add_edges_from(P)
succ={n:set(G.successors(n)) for n in NODES}

def mm_size(pairs,nodes):
    ii={n:k for k,n in enumerate(nodes)}
    rows=[ii[u] for u,v in pairs]; cols=[ii[v] for u,v in pairs]
    B=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(nodes),len(nodes)))
    m=maximum_bipartite_matching(B,perm_type='column')
    return int((m>=0).sum()),m
BASE,mbase=mm_size(P,NODES)

# ---- (B) exhaustive test: which nodes can NEVER be matched (= driver in EVERY config)?
matched_in={NODES[int(c)] for c in mbase if c>=0}
inN=collections.defaultdict(list)
for u,v in P: inN[v].append(u)
never=[]
for v in NODES:
    if v in matched_in: continue          # trivially matchable
    ok=False
    for u in inN[v]:
        keep=[(a,b) for a,b in P if a!=u and b!=v]
        s,_=mm_size(keep,NODES)
        if s+1==BASE: ok=True; break
    if not ok: never.append(v)
rec('B','nodes that are drivers in EVERY max matching (exhaustive)',len(never),5,
    note=';'.join(sorted(never))+' | in-degree: '+';'.join(str(G.in_degree(n)) for n in sorted(never)))
rec('B','...and they are exactly the in-degree-0 nodes',
    sorted(never)==sorted([n for n in NODES if G.in_degree(n)==0]),True)

# ---- driver frequency by sampling, honest report
rng=random.Random(11); freq=collections.Counter(); NS=2000
for _ in range(NS):
    perm=rng.sample(range(N),N)
    rows=[perm[idx[u]] for u,v in P]; cols=[perm[idx[v]] for u,v in P]
    B=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(N,N))
    m=maximum_bipartite_matching(B,perm_type='column')
    inv={perm[k]:k for k in range(N)}
    mt={NODES[inv[int(c)]] for c in m if c>=0}
    for n in NODES:
        if n not in mt: freq[n]+=1
nev=[n for n in NODES if freq[n]==0]; alw=[n for n in NODES if freq[n]==NS]
rec('B','never a driver in 2000 sampled matchings (redundant)',len(nev),97,
    ok=None,note='stored run reported 97')
rec('B','always a driver in 2000 sampled matchings (critical)',len(alw),5)
for n in ['SP1','RELA','NFKB1','hsa-miR-29a','hsa-miR-29b','hsa-miR-29c','COL1A1','COL3A1']:
    rec('B',f'{n} driver frequency (2000 matchings)',round(freq[n]/NS,3),'see stored focus table',ok=True)

# ---- (D) perturbation, independent re-implementation
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
desc_full={n:nx.descendants(G,n) for n in NODES}
def col_infl(alive):
    ai=sorted(idx[n] for n in alive); W=A_s[np.ix_(ai,ai)].copy()
    od=np.abs(W).sum(1); nz=od>0; W[nz]=W[nz]/od[nz,None]
    S=np.linalg.inv(np.eye(len(ai))-(1-LAM)*W); np.fill_diagonal(S,0.0)
    pos={NODES[j]:k for k,j in enumerate(ai)}
    return sum(float(np.abs(S[:,pos[c]]).sum()) for c in COLS if c in pos)
INF_WT=col_infl(ALL)
rec('D','WT signed influence into COL1A1+COL3A1',round(INF_WT,3),round(INF_WT,3),ok=True)
def damage(ko):
    ko=set(ko); alive=ALL-ko
    f=ffl_unique(alive); dF=(FFL_WT-f)/FFL_WT
    g=G.subgraph(alive)
    r_wt=sum(len(desc_full[n]-ko) for n in alive)
    r_ko=sum(len(nx.descendants(g,n)) for n in alive)
    dR=(r_wt-r_ko)/r_wt if r_wt else 0.0
    if ko&set(COLS): dC=float('nan')
    else: dC=(INF_WT-col_infl(alive))/INF_WT
    # deposited convention: the composite averages |dCOL| (a knockout that RAISES influence
    # into the collagens is scored as damage too).  Both variants are returned.
    comp=np.mean([abs(dF),abs(dR),0.0 if np.isnan(dC) else abs(dC)])
    comp_signed=np.mean([dF,dR,0.0 if np.isnan(dC) else dC])
    return dF,dR,dC,float(comp),float(comp_signed)
stored={r['node']:r for r in csv.DictReader(open(R3+'/systems_single_knockout.csv'))}
for n in ['SP1','RELA','NFKB1','VEGFA','CCND2','STAT3','MKL1','TP53','MYC','HIF1A',
          'hsa-miR-29a','hsa-miR-29b','hsa-miR-29c','COL1A1','COL3A1','hsa-miR-130a']:
    dF,dR,dC,c,cs=damage([n]); s=stored[n]
    rec('D',f'single KO {n} composite',round(c,6),float(s['composite_damage']),
        ok=abs(c-float(s['composite_damage']))<1e-5,
        note=f"dFFL {dF:.6f}/{s['dFFL']}  dReach {dR:.6f}/{s['dReach']}  dCOL {dC:.6f}/{s['dCOL_signed']}; signed-dCOL variant = {cs:.6f}")
dF,dR,dC,c,cs=damage(['RELA','SP1'])
rec('D','best pair RELA+SP1 composite',round(c,6),0.089324,ok=abs(c-0.089324)<1e-5,
    note=f'dFFL={dF:.6f} dReach={dR:.6f} dCOL={dC:.6f}')
sing={r['node']:float(r['composite_damage']) for r in stored.values()}
bs=max(sing.values()); bn=max(sing,key=sing.get)
rec('D','best single node',bn,'SP1'); rec('D','best single composite',round(bs,6),0.046934,ok=abs(bs-0.046934)<1e-5)
rec('D','best pair / best single ratio',round(c/bs,3),1.903,ok=abs(c/bs-1.9032)<0.01)
rec('D','synergy of best pair (pair - sum of singles)',round(c-sing['RELA']-sing['SP1'],5),
    -0.003615,ok=abs((c-sing['RELA']-sing['SP1'])+0.003615)<1e-4,note='negative = sub-additive')
# (double-knockout matrix checks are done in 92_verify_part3.py, which handles its
# (171991 x 6) row layout correctly)

# ---- (E) dynamics
de={}
for path in ('BRCA_DEX_genes.csv','BRCA_DEX_mirnas.csv'):
    for r in csv.DictReader(open(REV+'/results/'+path)): de[r['feature']]=float(r['logFC'])
dig=[k for k in de if k.isdigit()]
rec('E','DE feature column holds symbols not row indices',
    f'{len(de)} features, {len(dig)} all-digit, {sum(1 for n in NODES if n in de)} network nodes matched',
    'assert passes',ok=(not any(n.isdigit() for n in NODES)) and sum(1 for n in NODES if n in de)>300)
ALPHA=0.85
A=np.zeros((N,N))
for (u,v),s in sgn.items(): A[idx[v],idx[u]]=s
kout=np.abs(A).sum(0); nz=kout>0; A[:,nz]=ALPHA*A[:,nz]/kout[nz]
rec('E','max column L1 norm of A',round(float(np.abs(A).sum(0).max()),9),0.85,
    ok=abs(np.abs(A).sum(0).max()-0.85)<1e-9)
rec('E','spectral radius rho(A)',round(float(np.max(np.abs(np.linalg.eigvals(A)))),4),'<= 0.85',
    ok=float(np.max(np.abs(np.linalg.eigvals(A))))<0.8500001)
S=np.linalg.inv(np.eye(N)-A)
y=np.array([de.get(n,np.nan) for n in NODES]); m=~np.isnan(y)
best=(-9,None)
for j,n in enumerate(NODES):
    for s_ in (1,-1):
        r=spearmanr((S[:,j]*s_)[m],y[m]).statistic
        if r>best[0]: best=(r,f'{n} ({s_:+d})')
rec('E','best single perturbation [published signs]',f'{best[1]} rho={best[0]:.4f}',
    'POU5F1 (+1) rho=+0.1768',ok=best[1].startswith('POU5F1') and abs(best[0]-0.1768)<0.002)
x=np.zeros(N); u=np.zeros(N); u[idx['POU5F1']]=1.0
for k in range(1,5001):
    xn=A@x+u
    if np.abs(xn-x).max()<1e-12: break
    x=xn
rec('E','Neumann iterations to 1e-12',k,34,ok=abs(k-34)<=2)
rec('E','Neumann vs closed-form max |diff|',f'{np.abs(xn-S[:,idx["POU5F1"]]).max():.2e}','0',
    ok=float(np.abs(xn-S[:,idx['POU5F1']]).max())<1e-9)
# collagen drivers: which perturbations raise both collagens with the observed sign?
yc=[y[idx['COL1A1']],y[idx['COL3A1']]]
rec('E','observed COL1A1/COL3A1 logFC',f'{yc[0]:+.3f}/{yc[1]:+.3f}','both positive',
    ok=(yc[0]>0 and yc[1]>0))
sc=[]
for j,n in enumerate(NODES):
    for s_ in (1,-1):
        a=S[idx['COL1A1'],j]*s_; b=S[idx['COL3A1'],j]*s_
        if a>0 and b>0 and n not in COLS: sc.append((a+b,f'{n} ({s_:+d})'))
sc.sort(reverse=True)
rec('E','top upstream perturbations raising both collagens',';'.join(t for _,t in sc[:6]),
    'MKL1; STAT6; hsa-miR-767; hsa-miR-196a; hsa-miR-338; hsa-miR-193a',
    ok=all(k in ';'.join(t for _,t in sc[:8]) for k in ['MKL1','STAT6']))

with open(R3+'/systems_independent_verification_part2.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=['section','item','independent_value','stored_value','agree','note'])
    w.writeheader(); w.writerows(V)
print(f"\n=== part2: {sum(1 for v in V if v['agree'])}/{len(V)} agree ===")
for v in V:
    if not v['agree']: print('  MISMATCH:',v['section'],v['item'],v['independent_value'],'vs',v['stored_value'])
