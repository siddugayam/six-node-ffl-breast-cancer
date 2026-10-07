#!/usr/bin/env python3
"""Enumerate 3-node FFL cores WITHIN each of the six networks of the original submission
separately, after applying only the identifier harmonisation. This mirrors the per-network
analysis of the original submission and gives the per-network census."""
import csv, collections, json
REV='/path/to/revision'
REPO='/path/to/revision/analyses/original_submission_code'
nm=json.load(open(f'{REV}/data/name_map.json'))
merge=nm['merge_map']; ntype=nm['node_type']
def T(n): return ntype.get(merge.get(n,n),'Gene')
def C(n): return merge.get(n,n)

print(f"{'network':<10}{'edges':>8}{'nodes':>7}  {'miRNA-FFL':>10}{'TF-FFL':>8}{'Composite':>11}{'TOTAL':>8}{'nodes_in_FFL':>14}")
summary={}
for pfx in ['3-miR','3-TF','3-Comp','4-TF','5-TF','6-TF']:
    E=set()
    for l in open(f'{REPO}/SIF_files/{pfx}.sif'):
        p=l.rstrip('\n').split('\t')
        if len(p)<3: continue
        a,b=C(p[0]),C(p[2])
        if a!=b: E.add((a,b))
    out=collections.defaultdict(set)
    for a,b in E: out[a].add(b)
    cls=collections.Counter(); members=set()
    for R in list(out):
        for M in out[R]:
            for Tg in (out[R] & out.get(M,set())):
                if Tg in (R,M) or T(Tg)=='miRNA': continue
                if T(R)=='miRNA' and T(M)=='TF':
                    c='Composite-FFL' if (M,R) in E else 'miRNA-FFL'
                elif T(R)=='TF' and T(M)=='miRNA':
                    c='Composite-FFL' if (M,R) in E else 'TF-FFL'
                else: continue
                cls[c]+=1; members.update((R,M,Tg))
    tot=sum(cls.values())
    nodes=len({x for e in E for x in e})
    print(f"{pfx:<10}{len(E):>8}{nodes:>7}  {cls['miRNA-FFL']:>10}{cls['TF-FFL']:>8}{cls['Composite-FFL']:>11}{tot:>8}{len(members):>14}")
    summary[pfx]=dict(edges=len(E),nodes=nodes,**cls,total=tot,ffl_nodes=len(members))
json.dump(summary,open(f'{REV}/results/ffl_census_per_network_crosscheck.json','w'),indent=1)
print("\nNOTE: counts are after identifier harmonisation (precursor/mature collapse) and")
print("de-duplication; they are therefore NOT comparable to the raw row counts in Table 1.")
