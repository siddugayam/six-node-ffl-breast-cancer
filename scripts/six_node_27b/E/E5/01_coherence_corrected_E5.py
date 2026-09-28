#!/usr/bin/env python3
"""Coherent / incoherent typing (Mangan & Alon, PNAS 2003) on the CORRECTLY counted core set.

Numbering, keyed on (sign X->Y, sign Y->Z, sign X->Z):
  C1 (+,+,+)  C2 (-,-,+)  C3 (-,+,-)  C4 (+,-,-)
  I1 (+,-,+)  I2 (-,+,+)  I3 (-,-,-)  I4 (+,+,-)
A core is coherent iff sign(X->Z) == sign(X->Y)*sign(Y->Z).
Each composite core is counted ONCE (see 00_counting_convention.py)."""
import csv, collections, os
REV='/path/to/revision'

nodes={r['name']:r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t')}
# ---- E5 sandbox: optional re-typing of the 22 TF-typed nodes absent from the Lambert census as genes
import sys as _sys
T22='APEX1 BMI1 BRCA1 CREBBP CTNNB1 EP300 EZH2 HDAC1 HDAC2 HDAC3 HDAC4 HDAC9 ILF3 MEN1 MKL1 MTA1 NCOR1 NF1 RB1 SIRT1 SUZ12 VHL'.split()
RETYPE=len(_sys.argv)>1 and _sys.argv[1]=='retype'
if RETYPE:
    assert all(nodes[x]=='TF' for x in T22); nodes.update({x:'Gene' for x in T22})
print('E5 sandbox: 22 nodes re-typed as Gene' if RETYPE else 'E5 sandbox: original node types')
E=set(); ET={}
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
    E.add((r['source'],r['target'])); ET[(r['source'],r['target'])]=r['edge_type']
out=collections.defaultdict(set)
for a,b in E: out[a].add(b)

sign={}
def load(path,s,t,m):
    if not os.path.exists(path): return 0
    n=0
    for r in csv.DictReader(open(path),delimiter='\t'):
        v={'Activation':+1,'Repression':-1}.get((r.get(m) or '').strip())
        if v is not None: sign[(r[s],r[t])]=v; n+=1
    return n
n1=load(f'{REV}/data/layer_TF_target.tsv','source','target','mode')
n2=load(f'{REV}/data/layer_TF_miRNA.tsv','source','target','mode')

def sgn(e):
    if e in sign: return sign[e]
    if ET.get(e)=='miRNA_target': return -1          # miRNA repression is mechanistically certain
    return None                                       # unannotated TF edge: not assumed

recip={(t,m) for (t,m) in E if nodes.get(t)=='TF' and nodes.get(m)=='miRNA' and (m,t) in E}
cores=[]
for (t,m) in recip:                                   # composite, counted once, apex = TF
    for g in (out[t]&out[m]):
        if nodes.get(g)!='miRNA' and g not in (t,m): cores.append((t,m,g,'Composite-FFL'))
for (t,m) in E:                                       # TF-FFL
    if nodes.get(t)=='TF' and nodes.get(m)=='miRNA' and (t,m) not in recip:
        for g in (out[t]&out[m]):
            if nodes.get(g)!='miRNA' and g not in (t,m): cores.append((t,m,g,'TF-FFL'))
for (m,t) in E:                                       # miRNA-FFL
    if nodes.get(m)=='miRNA' and nodes.get(t)=='TF' and (t,m) not in E:
        for g in (out[m]&out[t]):
            if nodes.get(g)!='miRNA' and g not in (m,t): cores.append((m,t,g,'miRNA-FFL'))

print("class counts:", dict(collections.Counter(c[3] for c in cores)), " total:", len(cores))
NAME={(1,1,1):'C1',(-1,-1,1):'C2',(-1,1,-1):'C3',(1,-1,-1):'C4',
      (1,-1,1):'I1',(-1,1,1):'I2',(-1,-1,-1):'I3',(1,1,-1):'I4'}
coh=collections.Counter(); byclass=collections.defaultdict(collections.Counter)
rows=[]
for R,M,T,cls in cores:
    s1,s2,s3=sgn((R,M)),sgn((M,T)),sgn((R,T))
    lab='unresolved' if None in (s1,s2,s3) else NAME[(s1,s2,s3)]
    coh[lab]+=1; byclass[cls][lab]+=1
    rows.append([R,M,T,cls,nodes.get(R),nodes.get(M),nodes.get(T),s1,s2,s3,lab])

res=sum(v for k,v in coh.items() if k!='unresolved')
C=sum(v for k,v in coh.items() if k.startswith('C')); I=sum(v for k,v in coh.items() if k.startswith('I'))
print(f"\nresolved {res} of {len(cores)}  (unresolved {coh['unresolved']}, TF edge unannotated in TRRUST)")
for k in ['C1','C2','C3','C4','I1','I2','I3','I4']:
    if coh[k]: print(f"   {k}  {coh[k]:>5}   ({'coherent' if k[0]=='C' else 'incoherent'})")
print(f"   coherent {C} ({100*C/res:.1f}%)   incoherent {I} ({100*I/res:.1f}%)")
print("\nby motif class:")
for cls in byclass:
    d=byclass[cls]; r=sum(v for k,v in d.items() if k!='unresolved')
    c=sum(v for k,v in d.items() if k.startswith('C'))
    print(f"   {cls:<15} resolved {r:>5}  coherent {c} ({100*c/r:.1f}%)" if r else f"   {cls:<15} none resolved")
print("\nThe canonical TF-driven circuit (TF -> miRNA, miRNA -| gene, TF -> gene) is type",
      NAME[(1,-1,1)], "- incoherent type 1, a pulse generator / response accelerator.")
import os as _os
with open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),'e5b_cores_'+('retyped' if RETYPE else 'original')+'.csv'),'w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['regulator','intermediate','target','class','type_R','type_M','type_T','sign_RM','sign_MT','sign_RT','coherence'])
    w.writerows(rows)
print("\nwrote e5b_cores_*.csv (E5 sandbox)")
