#!/usr/bin/env python3
"""Hub stability by non-parametric bootstrap
. 1000 replicates; edges resampled with
replacement; degree and betweenness hub rankings recomputed in each replicate; selection
frequency in the top 20 within each node class reported."""
import csv, collections, random, sys
import numpy as np, networkx as nx

REV='/path/to/revision'
NB=int(sys.argv[1]) if len(sys.argv)>1 else 1000
rng=random.Random(20260908)

nodes={r['name']:r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t')}
edges=[(r['source'],r['target']) for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t') if True]
print(f"nodes={len(nodes)} edges={len(edges)} replicates={NB}")

def ranks(el):
    G=nx.DiGraph(); G.add_nodes_from(nodes); G.add_edges_from(el)
    deg=dict(G.degree())
    btw=nx.betweenness_centrality(G, k=min(120,len(nodes)), seed=7, normalized=True)
    out={}
    for metric,vals in (('degree',deg),('betweenness',btw)):
        for t in ('TF','miRNA','Gene'):
            sub=[(n,vals.get(n,0)) for n in nodes if nodes[n]==t]
            sub.sort(key=lambda x:-x[1])
            out[(metric,t)]=[n for n,_ in sub[:20]]
    return out

real=ranks(edges)
hits=collections.Counter()
for b in range(NB):
    samp=[edges[rng.randrange(len(edges))] for _ in range(len(edges))]
    r=ranks(list(set(samp)))
    for k,v in r.items():
        for n in v: hits[(k,n)]+=1
    if (b+1)%100==0: print(f"  {b+1}/{NB}")

rows=[]
for (metric,t),top in real.items():
    for rank,n in enumerate(top,1):
        rows.append(dict(node=n, node_type=t, metric=metric, real_rank=rank,
                         stability_pct=round(100*hits[((metric,t),n)]/NB,1), n_boot=NB))
rows.sort(key=lambda r:(r['metric'],r['node_type'],r['real_rank']))
with open('/path/to/revision/INBOX_2026-09-27/legacy_audit/hubboot/hub_bootstrap_stability_1000_full.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

print("\nmedian stability by metric x class:")
agg=collections.defaultdict(list)
for r in rows: agg[(r['metric'],r['node_type'])].append(r['stability_pct'])
for k,v in sorted(agg.items()):
    print(f"   {k[0]:<12}{k[1]:<7} median={np.median(v):5.1f}%  min={min(v):5.1f}%  n>=80%: {sum(1 for x in v if x>=80)}/20")
print("\ntop TF hubs (degree):")
for r in rows:
    if r['metric']=='degree' and r['node_type']=='TF' and r['real_rank']<=10:
        print(f"   {r['real_rank']:>2}. {r['node']:<14} {r['stability_pct']:5.1f}%")
print("\nwrote results/hub_bootstrap_stability_1000.csv")
