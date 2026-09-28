#!/usr/bin/env python3
"""
n-node feed-forward loop census by connected-induced-subgraph enumeration.

DEFINITION. In the directed regulatory graph G (reciprocal TF<->miRNA pairs contracted to
their transcriptional arc), a vertex-induced connected subgraph S on n vertices is an
n-node feed-forward loop iff
  (D1) S is acyclic;
  (D2) S has exactly one source and exactly one sink;
  (D3) every vertex of S lies on a directed source->sink path in S;
  (D4) S contains >=2 internally vertex-disjoint source->sink paths.
n=3 recovers the classical FFL.

Enumeration uses ESU (Wernicke 2006): exhaustive for small n, uniform-probability
sampling (RAND-ESU) for larger n with an unbiased estimator of the total count.
Usage: 03_ffl_census.py <n> [sampling_prob_per_depth, comma separated] [seed]
"""
import csv, collections, itertools, json, sys, os, random

REV='/path/to/revision'

def build_graph(extra_layers=True, mirna_thresh='threshold_10kb'):
    nodes={r['name']:r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t')}
    edges={}
    for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
        edges[(r['source'],r['target'])]=r['edge_type']
    if extra_layers:
        def add(path,sc,tc,et,undirected=False,filt=None):
            if not os.path.exists(path): return 0
            n=0
            for r in csv.DictReader(open(path),delimiter='\t'):
                if filt and not filt(r): continue
                a,b=r.get(sc),r.get(tc)
                if a not in nodes or b not in nodes or a==b: continue
                for e in ([(a,b),(b,a)] if undirected else [(a,b)]):
                    if e not in edges: edges[e]=et; n+=1
            return n
        add(f'{REV}/data/layer_TF_target.tsv','source','target','TF_target')
        add(f'{REV}/data/layer_gene_gene.tsv','source','target','gene_gene')
        add(f'{REV}/data/layer_miRNA_miRNA.tsv','miRNA_1','miRNA_2','miRNA_miRNA',
            undirected=True, filt=lambda r: str(r.get(mirna_thresh,'TRUE')).upper() in ('TRUE','1'))
    # contract reciprocal TF<->miRNA
    nrecip=0
    for (a,b) in list(edges):
        if nodes.get(a)=='TF' and nodes.get(b)=='miRNA' and (b,a) in edges:
            if edges.pop((b,a),None) is not None: nrecip+=1
    return nodes, edges, nrecip

def is_ffl(S, out_s, in_s):
    src=[v for v in S if not in_s[v]]; snk=[v for v in S if not out_s[v]]
    if len(src)!=1 or len(snk)!=1: return None
    s,t=src[0],snk[0]
    # acyclic
    colour={}
    def visit(u):
        colour[u]=1
        for v in out_s[u]:
            c=colour.get(v)
            if c==1: return False
            if c is None and not visit(v): return False
        colour[u]=2; return True
    for v in S:
        if colour.get(v) is None and not visit(v): return None
    # D3: every vertex reachable from s and reaching t
    reach={s}; st=[s]
    while st:
        u=st.pop()
        for v in out_s[u]:
            if v not in reach: reach.add(v); st.append(v)
    back={t}; st=[t]
    while st:
        u=st.pop()
        for v in in_s[u]:
            if v not in back: back.add(v); st.append(v)
    if not S <= (reach & back): return None
    # D4: >=2 internally vertex-disjoint s->t paths  == max-flow >= 2 with unit node caps
    # small n: test by finding one path, deleting its internals, seeking another
    def find_path(banned):
        prev={s:None}; st=[s]
        while st:
            u=st.pop()
            if u==t: break
            for v in out_s[u]:
                if v in prev or (v in banned and v!=t): continue
                prev[v]=u; st.append(v)
        if t not in prev: return None
        p=[]; u=t
        while u is not None: p.append(u); u=prev[u]
        return p[::-1]
    p1=find_path(set())
    if p1 is None: return None
    if find_path(set(p1[1:-1])) is None: return None
    return (s,t)

def esu(nodes_list, adj_u, k, probs, rng, callback):
    """RAND-ESU. probs[d] = probability of expanding at depth d (1.0 = exhaustive)."""
    idx={v:i for i,v in enumerate(nodes_list)}
    for v in nodes_list:
        if probs[0]<1.0 and rng.random()>=probs[0]: continue
        ext={u for u in adj_u[v] if idx[u]>idx[v]}
        _extend([v], ext, v, adj_u, idx, k, probs, rng, callback)

def _extend(sub, ext, v, adj_u, idx, k, probs, rng, callback):
    if len(sub)==k:
        callback(sub); return
    d=len(sub)
    ext=set(ext)
    while ext:
        w=ext.pop()
        if probs[d]<1.0 and rng.random()>=probs[d]: continue
        excl=set(sub)|{u for s in sub for u in adj_u[s]}
        newext=ext | {u for u in adj_u[w] if idx[u]>idx[v] and u not in excl}
        _extend(sub+[w], newext, v, adj_u, idx, k, probs, rng, callback)

def census(k, probs, seed=1, extra_layers=True, verbose=True):
    nodes, edges, nrecip = build_graph(extra_layers)
    out=collections.defaultdict(set); inn=collections.defaultdict(set); adj_u=collections.defaultdict(set)
    for a,b in edges:
        out[a].add(b); inn[b].add(a); adj_u[a].add(b); adj_u[b].add(a)
    nl=[v for v in nodes if adj_u[v]]
    if verbose:
        print(f"graph: {len(nl)} connected nodes, {len(edges)} directed edges, {nrecip} reciprocal pairs contracted")
    rng=random.Random(seed)
    stats=collections.Counter(); comp=collections.Counter(); found=[]
    classdiv=collections.Counter(); classsets=collections.Counter()
    def cb(sub):
        S=set(sub)
        os_={v:{u for u in out[v] if u in S} for v in S}
        is_={v:{u for u in inn[v] if u in S} for v in S}
        stats['subgraphs']+=1
        r=is_ffl(S, os_, is_)
        if r:
            stats['ffl']+=1
            tc=collections.Counter(nodes[v] for v in S)
            comp[(tc['TF'],tc['miRNA'],tc['Gene'])]+=1
            def _cls(a,b):
                t=edges[(a,b)]
                if t=='TF_target': return 'TF->TF' if nodes[b]=='TF' else 'TF->Gene'
                if t=='miRNA_target': return 'miRNA->TF' if nodes[b]=='TF' else 'miRNA->Gene'
                return {'TF_miRNA':'TF->miRNA','gene_gene':'Gene->Gene','miRNA_miRNA':'miRNA-miRNA'}.get(t,t)
            ec=frozenset(_cls(a,b) for a in S for b in os_[a])
            classdiv[len(ec)]+=1; classsets[tuple(sorted(ec))]+=1
            if len(found)<2000: found.append((r[0],r[1],sorted(S)))
    esu(nl, adj_u, k, probs, rng, cb)
    return stats, comp, found, nodes, classdiv, classsets

if __name__=='__main__':
    k=int(sys.argv[1])
    probs=[float(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [1.0]*k
    seed=int(sys.argv[3]) if len(sys.argv)>3 else 1
    st,comp,found,nodes,classdiv,classsets=census(k,probs,seed)
    scale=1.0
    for p in probs: scale/=p
    print(f"\nk={k}  probs={probs}  sampled connected subgraphs={st['subgraphs']}  FFLs among them={st['ffl']}")
    print(f"estimated total connected induced {k}-subgraphs: {st['subgraphs']*scale:,.0f}")
    print(f"estimated total {k}-node FFLs:                   {st['ffl']*scale:,.0f}")
    print("\ncomposition (TF, miRNA, Gene) -> estimated count:")
    for kk,v in comp.most_common(12):
        print(f"   TF={kk[0]} miRNA={kk[1]} Gene={kk[2]}  ->  {v*scale:,.0f}   (sampled {v})")
    print("\ndistinct edge classes per FFL module (n_classes -> estimated modules):")
    for nc in sorted(classdiv): print(f"   {nc} classes -> {classdiv[nc]*scale:,.0f}  (sampled {classdiv[nc]})")
    print("\ntop edge-class combinations:")
    for cs,v in classsets.most_common(8): print(f"   {'+'.join(cs):<62} {v*scale:,.0f}")
    json.dump({'k':k,'class_diversity':{str(a):b*scale for a,b in classdiv.items()},
               'class_sets':{'+'.join(a):b*scale for a,b in classsets.items()},'probs':probs,'sampled_subgraphs':st['subgraphs'],'sampled_ffl':st['ffl'],
               'scale':scale,'est_subgraphs':st['subgraphs']*scale,'est_ffl':st['ffl']*scale,
               'composition':{f"TF{a}_miR{b}_G{c}":v*scale for (a,b,c),v in comp.items()}},
              open(f'{REV}/results/ffl_census_k{k}.json','w'), indent=1)
    with open(f'{REV}/results/ffl_examples_k{k}.csv','w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['source','sink','members','types'])
        for s,t,S in found[:2000]:
            w.writerow([s,t,';'.join(S),';'.join(nodes[x] for x in S)])
    print(f"\nwrote results/ffl_census_k{k}.json and ffl_examples_k{k}.csv")
