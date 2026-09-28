#!/usr/bin/env python3
"""
Formal enumeration of n-node feed-forward loops.

DEFINITION (final). Let G be the directed regulatory graph in which each reciprocal
TF<->miRNA pair is represented by its transcriptional arc TF->miRNA (the reciprocity is
retained as an edge attribute). A vertex-induced subgraph S on n vertices is an
*n-node feed-forward loop* iff

  (D1) S is a directed acyclic graph;
  (D2) S has exactly one source s (in-degree 0 in S) and exactly one sink t (out-degree 0 in S);
  (D3) every vertex of S lies on at least one directed s->t path within S;
  (D4) S contains at least two internally vertex-disjoint directed s->t paths.

n = 3 recovers the classical FFL exactly. (D3) excludes the "cascade with a side branch"
routes from a single input to a single output.

Enumeration strategy: a subgraph satisfying D2-D4 is exactly the vertex set of a union of
>=2 internally vertex-disjoint s->t paths. We therefore enumerate simple s->t paths
(bounded length), combine disjoint path sets, and verify D1-D3 on the induced subgraph.
"""
import csv, collections, itertools, json, sys, os
REV='/path/to/revision'
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 6

# ---------------------------------------------------------------- build graph
nodes={r['name']:r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t')}
edges={}
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
    edges[(r['source'],r['target'])]={'etype':r['etype'] if 'etype' in r else r['edge_type']}

def add_layer(path, s_col, t_col, etype, undirected=False, filt=None):
    if not os.path.exists(path): return 0
    n=0
    for r in csv.DictReader(open(path),delimiter='\t'):
        if filt and not filt(r): continue
        a,b=r[s_col],r[t_col]
        if a not in nodes or b not in nodes or a==b: continue
        for e in ([(a,b),(b,a)] if undirected else [(a,b)]):
            if e not in edges: edges[e]={'etype':etype}; n+=1
    return n
n_tt = add_layer(f'{REV}/data/layer_TF_target.tsv','source','target','TF_target')
n_gg = add_layer(f'{REV}/data/layer_gene_gene.tsv','source','target','gene_gene')
n_mm = add_layer(f'{REV}/data/layer_miRNA_miRNA.tsv','miRNA_1','miRNA_2','miRNA_miRNA',
                 undirected=True, filt=lambda r: r.get('threshold_10kb','TRUE') in ('TRUE','True','1'))
print(f"layers added: TF_target +{n_tt}  gene_gene +{n_gg}  miRNA_miRNA(10kb) +{n_mm}")

# contract reciprocal TF<->miRNA to the transcriptional arc
recip=set()
for (a,b) in list(edges):
    if nodes.get(a)=='TF' and nodes.get(b)=='miRNA' and (b,a) in edges:
        recip.add((a,b)); edges.pop((b,a),None)
print(f"reciprocal TF<->miRNA pairs contracted: {len(recip)}")

out=collections.defaultdict(set); inn=collections.defaultdict(set)
for a,b in edges: out[a].add(b); inn[b].add(a)
print(f"graph: {len(nodes)} nodes, {len(edges)} directed edges")
outdeg=collections.Counter({n:len(out[n]) for n in nodes})
by=collections.defaultdict(list)
for n in nodes: by[nodes[n]].append(len(out[n]))
for k,v in by.items(): print(f"   out-degree {k:<6} mean={sum(v)/len(v):.1f} max={max(v)} zero={sum(1 for x in v if x==0)}")

# ---------------------------------------------------------------- enumerate
MAXINT = NMAX-2                      # max internal nodes across all paths of a module
paths_by_st=collections.defaultdict(list)
def dfs(s, cur, seen):
    u=cur[-1]
    for v in out[u]:
        if v in seen: continue
        if len(cur)-1 >= MAXINT+1:    # cur has s + internals; cap path edge count
            continue
        p=cur+[v]
        paths_by_st[(s,v)].append(tuple(p))
        if len(p)-2 < MAXINT:
            dfs(s, p, seen|{v})
sys.setrecursionlimit(10000)
for s in nodes:
    if out[s]: dfs(s,[s],{s})
tot_paths=sum(len(v) for v in paths_by_st.values())
print(f"simple paths enumerated (<= {MAXINT+1} edges): {tot_paths} over {len(paths_by_st)} (s,t) pairs")

def induced_ok(S, s, t):
    sub={(a,b) for (a,b) in edges if a in S and b in S}
    ind=collections.defaultdict(set); ino=collections.Counter(); outo=collections.Counter()
    for a,b in sub: ind[a].add(b); outo[a]+=1; ino[b]+=1
    if [v for v in S if ino[v]==0]!=[s] and set(v for v in S if ino[v]==0)!={s}: return None
    if set(v for v in S if outo[v]==0)!={t}: return None
    # acyclicity
    colour={}
    def visit(u):
        colour[u]=1
        for v in ind[u]:
            if colour.get(v)==1: return False
            if colour.get(v) is None and not visit(v): return False
        colour[u]=2; return True
    for v in S:
        if colour.get(v) is None and not visit(v): return None
    return sub

census=collections.Counter(); arch=collections.Counter(); members=collections.defaultdict(set)
examples=collections.defaultdict(list)
for (s,t),plist in paths_by_st.items():
    if len(plist)<2: continue
    plist=sorted(set(plist), key=len)
    for k in (2,3):
        for combo in itertools.combinations(plist,k):
            ints=[set(p[1:-1]) for p in combo]
            if any(ints[i]&ints[j] for i in range(k) for j in range(i+1,k)): continue
            S=set().union(*[set(p) for p in combo])
            n=len(S)
            if n<3 or n>NMAX: continue
            sub=induced_ok(S,s,t)
            if sub is None: continue
            key=frozenset(S)
            if key in members[n]: continue
            members[n].add(key)
            census[n]+=1
            et=collections.Counter(edges[e]['etype'] for e in sub)
            arch[(n,tuple(sorted(et)))]+=1
            if len(examples[n])<5: examples[n].append((s,t,sorted(S)))

print("\n=== n-node FFL census ===")
for n in range(3,NMAX+1):
    ns=set().union(*members[n]) if members[n] else set()
    print(f"   n={n}: {census[n]:>8} modules   distinct nodes involved: {len(ns)}")
print(f"   TOTAL: {sum(census.values())}")
print("\nedge-class composition by size (top architectures):")
for (n,cls),v in sorted(arch.items(), key=lambda x:(-x[1]))[:14]:
    print(f"   n={n} {'+'.join(cls):<52} {v}")
for n in range(3,NMAX+1):
    if examples[n]: print(f"\n example n={n}: source={examples[n][0][0]} sink={examples[n][0][1]} nodes={examples[n][0][2]}")

json.dump({str(n):census[n] for n in range(3,NMAX+1)}, open(f'{REV}/results/ffl_census_crosscheck.json','w'), indent=1)
with open(f'{REV}/results/ffl_modules_crosscheck.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['n_nodes','members'])
    for n in range(3,NMAX+1):
        for S in members[n]: w.writerow([n,';'.join(sorted(S))])
print("\nwrote results/ffl_census_crosscheck.json and ffl_modules_crosscheck.csv")
