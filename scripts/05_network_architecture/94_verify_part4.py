#!/usr/bin/env python3
"""Part 4: degree-adjusted controllability control, Boolean attractor, C-section spot checks."""
import csv, json, collections
import numpy as np, networkx as nx
from scipy.stats import fisher_exact, binomtest, spearmanr
import scipy.sparse as sp
from scipy.sparse.csgraph import maximum_bipartite_matching
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

# ---- FFL participation and hub definition, recomputed
part=collections.Counter()
for Rn in NODES:
    tR=ntype[Rn]
    for M in succ[Rn]:
        tM=ntype[M]
        if not ((tR=='miRNA' and tM=='TF') or (tR=='TF' and tM=='miRNA')): continue
        for T in succ[Rn]&succ[M]:
            if T in (Rn,M) or ntype[T]=='miRNA': continue
            part[Rn]+=1; part[M]+=1; part[T]+=1
stored={r['node']:r for r in csv.DictReader(open(R3+'/systems_controllability_nodes.csv'))}
mx=max(abs(part[n]-int(stored[n]['ffl_participation'])) for n in NODES)
rec('B','FFL participation per node matches stored',f'max diff {mx}',0,ok=mx==0)
hub=[n for n in NODES if stored[n].get('is_FFL_hub','').upper()=='TRUE']
rec('B','FFL hubs defined in the stored table',len(hub),60,ok=len(hub)==60,
    note='definition taken from the stored table; participation itself independently reproduced')

# ---- driver set from an independent maximum matching, then degree-adjusted logistic
ii={n:k for k,n in enumerate(NODES)}
B=sp.csr_matrix((np.ones(len(P)),([ii[u] for u,v in P],[ii[v] for u,v in P])),shape=(N,N))
m=maximum_bipartite_matching(B,perm_type='column')
matched={NODES[int(c)] for c in m if c>=0}
drv=[n for n in NODES if n not in matched]
rec('B','driver nodes in one maximum matching',len(drv),206)
a=sum(1 for n in hub if n in drv); b=len(hub)-a
c=len(drv)-a; d=N-len(hub)-c
orr,pf=fisher_exact([[a,b],[c,d]])
rec('B','FFL hub -> driver node (independent matching)',f'{a}/{len(hub)} OR={orr:.3f} p={pf:.2g}',
    '8/60 OR=0.256 p=9.8e-05',ok=(a/len(hub)<0.25 and pf<0.05),
    note='the identity of drivers is matching-dependent; the DEPLETION is what replicates')
# degree-adjusted logistic
try:
    import statsmodels.api as smapi
    deg=np.array([G.degree(n) for n in NODES],float)
    X=np.column_stack([np.log(deg+1),[1.0 if ntype[n]=='TF' else 0 for n in NODES],
                       [1.0 if ntype[n]=='miRNA' else 0 for n in NODES],
                       [1.0 if n in hub else 0 for n in NODES]])
    X=smapi.add_constant(X)
    yv=np.array([1.0 if n in drv else 0.0 for n in NODES])
    fit=smapi.Logit(yv,X).fit(disp=0)
    rec('B','logistic driver ~ log(degree)+type+is_FFL_hub: coefficient on is_FFL_hub',
        f'{fit.params[-1]:+.3f} p={fit.pvalues[-1]:.3f}','-0.559 p=0.248',
        ok=(fit.pvalues[-1]>0.05),
        note='KEY CONTROL: the FFL-hub driver depletion is NOT significant once degree is adjusted for')
except Exception as e:
    rec('B','degree-adjusted logistic',f'FAILED: {e}','-0.559 p=0.248',ok=False)

# ---- (C) information flow spot checks against the stored table
inf={r['node']:r for r in csv.DictReader(open(R3+'/systems_information_flow.csv'))}
print('   information-flow columns:',list(next(iter(inf.values())).keys()))
U=G.to_undirected()
Ul=U.subgraph(max(nx.connected_components(U),key=len)).copy()
cfb=nx.current_flow_betweenness_centrality(Ul,normalized=True)
nl=sorted(Ul.nodes())
col=[k for k in next(iter(inf.values())) if 'current' in k.lower()]
if col:
    st=[float(inf[n][col[0]]) for n in nl]
    r=spearmanr(st,[cfb[n] for n in nl]).statistic
    rec('C','stored current-flow betweenness vs independent recomputation',round(r,4),1.0,
        ok=r>0.999,note=f'column {col[0]}')

# ---- (E) Boolean threshold model, independent re-implementation
de={};fdr={}
for p in ('BRCA_DEX_genes.csv','BRCA_DEX_mirnas.csv'):
    for r in csv.DictReader(open(REV+'/results/'+p)):
        de[r['feature']]=float(r['logFC']); fdr[r['feature']]=float(r['adj.P.Val'])
y=np.array([de.get(n,np.nan) for n in NODES]); f=np.array([fdr.get(n,np.nan) for n in NODES])
meas=~np.isnan(y); sig=meas&(f<0.05)&(np.abs(y)>0.5)
rec('E','strongly DE network nodes (FDR<0.05, |logFC|>0.5)',int(sig.sum()),336)
rows=np.array([idx[u] for u,v in P]); cols=np.array([idx[v] for u,v in P])
vals=np.array([float(sgn[(u,v)]) for u,v in P])
kout=np.bincount(rows,weights=np.abs(vals),minlength=N)
w=vals/np.where(kout[rows]>0,kout[rows],1.0)
indeg0=np.bincount(cols,minlength=N)==0
def brun(seed,maxsteps=200):
    x=seed.copy(); hist={}
    for t in range(maxsteps):
        k=x.tobytes()
        if k in hist: return x,t,t-hist[k]
        hist[k]=t
        ax=np.bincount(cols,weights=w*x[rows],minlength=N)
        xn=np.sign(ax); xn[indeg0]=x[indeg0]; x=xn
    return x,maxsteps,-1
s0=np.zeros(N); s0[sig]=np.sign(y[sig])
att,steps,period=brun(s0)
rec('E','Boolean attractor: steps / period',f'{steps} / {period}','12 / 4',
    ok=(steps==12 and period==4))
ag=(np.sign(att[sig])==np.sign(y[sig])); nz=np.abs(att[sig])>0
rec('E','attractor / phenotype sign agreement',f'{ag[nz].mean():.4f} on n={int(nz.sum())}',
    '0.7113 on n=336',ok=abs(ag[nz].mean()-0.7113)<0.002)
bt=binomtest(int(ag[nz].sum()),int(nz.sum()),0.5)
rec('E','binomial p vs 0.5',f'{bt.pvalue:.2e}','5.68e-15',ok=bt.pvalue<1e-12)
# permuted-seed null, independent
rng=np.random.default_rng(5); nullv=[]
sidx=np.where(sig)[0]
for _ in range(2000):
    sp_=np.zeros(N); sp_[sidx]=rng.permutation(np.sign(y[sig]))
    a2,_,_=brun(sp_)
    g=(np.sign(a2[sig])==np.sign(y[sig])); z=np.abs(a2[sig])>0
    if z.sum(): nullv.append(g[z].mean())
nullv=np.array(nullv)
rec('E','permuted-seed null mean/sd',f'{nullv.mean():.3f}/{nullv.std():.3f}','0.503/0.110',
    ok=abs(nullv.mean()-0.503)<0.03)
rec('E','p vs permuted-seed null',f'{(1+(nullv>=ag[nz].mean()).sum())/(1+len(nullv)):.4f}','0.0005',
    ok=(1+(nullv>=ag[nz].mean()).sum())/(1+len(nullv))<0.01)

with open(R3+'/systems_independent_verification_part4.csv','w',newline='') as fh:
    wtr=csv.DictWriter(fh,fieldnames=['section','item','independent_value','stored_value','agree','note'])
    wtr.writeheader(); wtr.writerows(V)
print(f"\n=== part4: {sum(1 for v in V if v['agree'])}/{len(V)} agree ===")
for v in V:
    if not v['agree']: print('  MISMATCH:',v['section'],v['item'],v['independent_value'],'vs',v['stored_value'])
