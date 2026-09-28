#!/usr/bin/env python3
"""
Motif over-representation against randomised networks.

Two null models are reported, because they answer different questions:

  NULL-A  degree- and class-preserving randomisation. Within each edge class the bipartite
          source->target structure is randomised by the curveball algorithm (Strona et al.
          2014), which is an exact, uniform sampler over the space of binary matrices with
          the given row and column sums. In/out degree per class is preserved exactly.
          This null DESTROYS reciprocal TF<->miRNA pairs.

  NULL-B  as NULL-A, but the reciprocal TF<->miRNA pairs are held fixed and only the
          non-reciprocal TF->miRNA and miRNA->TF edges are randomised. This is the
          appropriate null for asking whether composite motifs are over-represented GIVEN
          the reciprocity that the study's own miRNA-TF filtering step imposes; NULL-A
          confounds the two.

Counting is exact and vectorised: for a boolean adjacency A, the number of common targets
of every ordered pair is (A @ A.T), so the number of FFL cores over an edge set S is
sum_{(u,v) in S} (A @ A.T)[u,v], corrected for the degenerate cases T==u and T==v.
"""
import csv, json, sys, collections
import numpy as np

REV='/path/to/revision'
NRAND=int(sys.argv[1]) if len(sys.argv)>1 else 1000
rng=np.random.default_rng(20260908)

nodes={}; order=[]
for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t'):
    nodes[r['name']]=r['type']; order.append(r['name'])
idx={n:i for i,n in enumerate(order)}; N=len(order)
typ=np.array([nodes[n] for n in order])
is_mir=(typ=='miRNA'); is_tf=(typ=='TF'); is_gene=(typ=='Gene')

cls_edges=collections.defaultdict(list)
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
    if r['edge_type']=='miRNA_miRNA': continue        # F5 sandbox: analysed network (6,829 edges)
    cls_edges[r['edge_type']].append((idx[r['source']],idx[r['target']]))
print("edge classes:", {k:len(v) for k,v in cls_edges.items()})

def build_A(cls):
    A=np.zeros((N,N),dtype=bool)
    for es in cls.values():
        for u,v in es: A[u,v]=True
    return A

def curveball(edges, n_iter=None):
    """Uniform degree-preserving randomisation of a bipartite edge list."""
    adj=collections.defaultdict(set)
    for u,v in edges: adj[u].add(v)
    keys=list(adj)
    if len(keys)<2: return [(u,v) for u in adj for v in adj[u]]
    n_iter = n_iter or 5*len(keys)
    for _ in range(n_iter):
        a,b=rng.choice(len(keys),2,replace=False)
        A_,B_=keys[a],keys[b]
        sa,sb=adj[A_],adj[B_]
        common=sa & sb
        ex_a=list(sa-common); ex_b=list(sb-common)
        pool=ex_a+ex_b
        if not pool: continue
        rng.shuffle(pool)
        na=len(ex_a)
        adj[A_]=common | set(pool[:na])
        adj[B_]=common | set(pool[na:])
    return [(u,v) for u in adj for v in adj[u]]

def count_motifs(A, cls):
    """Return dict of FFL-core counts by class, given adjacency A and class edge lists."""
    M=(A.astype(np.int16) @ A.astype(np.int16).T)          # common out-neighbour counts
    res={}
    # reciprocal / non-reciprocal TF<->miRNA
    tfm=np.array(cls.get('TF_miRNA',[]),dtype=int).reshape(-1,2)
    mt =np.array(cls.get('miRNA_target',[]),dtype=int).reshape(-1,2)
    mir_tf = mt[is_tf[mt[:,1]]] if len(mt) else mt          # miRNA -> TF edges
    tfm_set=set(map(tuple,tfm.tolist()))
    mirtf_set=set(map(tuple,mir_tf.tolist()))
    recip={(t,m) for (t,m) in tfm_set if (m,t) in mirtf_set}
    def core_count(pairs, require_gene_target=True):
        tot=0
        for u,v in pairs:
            c=int(M[u,v])
            if A[u,v] and A[v,v]: pass
            if A[u,u] and A[v,u]: c-=1                      # T == u
            if A[u,v] and A[v,v]: c-=1                      # T == v
            if require_gene_target:
                c=int(np.count_nonzero(A[u] & A[v] & ~is_mir))
                if A[u,v] and A[v,v]: c-=1
            tot+=c
        return tot
    res['Composite-FFL']=core_count(recip)
    res['TF-FFL']       =core_count(tfm_set-recip)
    res['miRNA-FFL']    =core_count({(m,t) for (m,t) in mirtf_set if (t,m) not in tfm_set})
    res['Any-3node-FFL']=res['Composite-FFL']+res['TF-FFL']+res['miRNA-FFL']
    return res

A0=build_A(cls_edges)
real=count_motifs(A0, cls_edges)
print("\nOBSERVED:", real)

def randomise(cls, keep_recip):
    new={}
    if not keep_recip:
        for k,es in cls.items(): new[k]=curveball(es)
        return new
    tfm=set(map(tuple,cls.get('TF_miRNA',[])))
    mt =list(map(tuple,cls.get('miRNA_target',[])))
    mirtf={(u,v) for u,v in mt if is_tf[v]}
    recip={(t,m) for (t,m) in tfm if (m,t) in mirtf}
    recip_rev={(m,t) for (t,m) in recip}
    for k,es in cls.items():
        es=set(map(tuple,es))
        if k=='TF_miRNA':
            free=es-recip; new[k]=list(recip)+curveball(list(free))
        elif k=='miRNA_target':
            free=es-recip_rev; new[k]=list(recip_rev)+curveball(list(free))
        else:
            new[k]=curveball(list(es))
    return new

rows=[]
for null_name, keep in [('NULL-A_full_randomisation',False), ('NULL-B_reciprocity_preserved',True)]:
    acc=collections.defaultdict(list)
    for i in range(NRAND):
        rc=randomise(cls_edges, keep)
        Ar=build_A(rc)
        c=count_motifs(Ar, rc)
        for k,v in c.items(): acc[k].append(v)
        if (i+1)%200==0: print(f"  {null_name} {i+1}/{NRAND}")
    for k in real:
        arr=np.array(acc[k],dtype=float)
        mu,sd=arr.mean(),arr.std(ddof=1)
        z=(real[k]-mu)/sd if sd>0 else np.nan
        p_over=(np.sum(arr>=real[k])+1)/(NRAND+1)
        p_under=(np.sum(arr<=real[k])+1)/(NRAND+1)
        rows.append(dict(null_model=null_name, motif_class=k, n_real=real[k],
                         rand_mean=round(mu,2), rand_sd=round(sd,2), Z=round(z,3),
                         fold_change=round(real[k]/mu,4) if mu>0 else np.nan,
                         p_over=p_over, p_under=p_under,
                         rand_min=int(arr.min()), rand_max=int(arr.max()), n_rand=NRAND))
        print(f"{null_name:<32} {k:<16} real={real[k]:<7} rand={mu:8.1f}+-{sd:6.1f}  Z={z:8.2f}  fold={real[k]/mu if mu>0 else float('nan'):6.3f}  p_over={p_over:.4f}")

import os
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'motif_significance_1000_analysed_network.csv'),'w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\nwrote motif_significance_1000_analysed_network.csv (F5 sandbox)")
