#!/usr/bin/env python3
# Sub-type numbering of Mangan and Alon (2003), signs (R->M, M->T, R->T); relabelled 2026-09-26 from an earlier non-standard numbering.
"""Independent cross-check of the 3-node FFL census and coherence typing.
Deliberately written separately from the census agent's implementation so the two
can be compared. Uses only the canonical network + TRRUST/TransmiR signs."""
import csv, collections, itertools, os, sys
REV='/path/to/revision'

nodes={r['name']:r['type'] for r in csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'),delimiter='\t')}
E=[]; ET={}
for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'),delimiter='\t'):
    E.append((r['source'],r['target'])); ET[(r['source'],r['target'])]=r['edge_type']

# --- signs from TRRUST / TransmiR where available
sign={}
def load_layer(path, key_s, key_t, key_mode):
    if not os.path.exists(path): return 0
    n=0
    for r in csv.DictReader(open(path),delimiter='\t'):
        m=(r.get(key_mode) or '').strip()
        s={'Activation':+1,'Repression':-1}.get(m)
        if s is not None:
            sign[(r[key_s],r[key_t])]=s; n+=1
    return n
n1=load_layer(f'{REV}/data/layer_TF_target.tsv','source','target','mode')
n2=load_layer(f'{REV}/data/layer_TF_miRNA.tsv','source','target','mode')
print(f"signed edges from TRRUST={n1} TransmiR={n2}")

def sgn(e):
    if e in sign: return sign[e], 'annotated'
    t=ET.get(e)
    if t=='miRNA_target': return -1,'inferred'      # miRNA repression is mechanistically certain
    return None,'unknown'                            # TF edges without annotation: not assumed

out=collections.defaultdict(set); inn=collections.defaultdict(set)
Eset=set(E)
for a,b in E: out[a].add(b); inn[b].add(a)

# --- enumerate 3-node cores: R->M, R->T, M->T
cores=[]
for R in list(out):
    for M in out[R]:
        common = out[R] & out[M]
        for T in common:
            if T==R or T==M: continue
            cores.append((R,M,T))
print("raw 3-node cores (any node types):", len(cores))

tR=lambda n: nodes.get(n,'Gene')
def classify(R,M,T):
    if tR(T)=='miRNA': return None                      # target must be a gene/TF
    if tR(R)=='miRNA' and tR(M)=='TF':
        return 'Composite-FFL' if (M,R) in Eset else 'miRNA-FFL'
    if tR(R)=='TF' and tR(M)=='miRNA':
        return 'Composite-FFL' if (M,R) in Eset else 'TF-FFL'
    return None

cls=collections.Counter(); typed=[]
for R,M,T in cores:
    c=classify(R,M,T)
    if c: cls[c]+=1; typed.append((R,M,T,c))
print("\nTyped 3-node FFL cores:")
for k,v in cls.most_common(): print(f"   {k:<16} {v}")
print(f"   {'TOTAL':<16} {sum(cls.values())}")
print("distinct nodes in cores:", len({x for t in typed for x in t[:3]}))

# --- coherence typing (Mangan & Alon)
coh=collections.Counter()
for R,M,T,c in typed:
    s1,k1=sgn((R,M)); s2,k2=sgn((M,T)); s3,k3=sgn((R,T))
    if None in (s1,s2,s3): coh['unresolved (TF edge unannotated)']+=1; continue
    direct=s3; indirect=s1*s2
    coherent = (direct==indirect)
    # Mangan-Alon numbering by the sign triple (X->Y, Y->Z, X->Z)
    key=(s1,s2,s3)
    names={(1,1,1):'C1',(-1,1,-1):'C2',(1,-1,-1):'C3',(-1,-1,1):'C4',
           (1,-1,1):'I1',(-1,-1,-1):'I2',(1,1,-1):'I3',(-1,1,1):'I4'}
    coh[f"{names[key]} ({'coherent' if coherent else 'incoherent'})"]+=1
print("\nCoherence types:")
for k,v in coh.most_common(): print(f"   {k:<34} {v}")

with open(f'{REV}/results/ffl_3node_cores_crosscheck.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['regulator','intermediate','target','class',
        'type_R','type_M','type_T','sign_RM','sign_MT','sign_RT','coherence'])
    for R,M,T,c in typed:
        s1,_=sgn((R,M)); s2,_=sgn((M,T)); s3,_=sgn((R,T))
        names={(1,1,1):'C1',(-1,1,-1):'C2',(1,-1,-1):'C3',(-1,-1,1):'C4',
           (1,-1,1):'I1',(-1,-1,-1):'I2',(1,1,-1):'I3',(-1,1,1):'I4'}
        lab = names.get((s1,s2,s3),'unresolved') if None not in (s1,s2,s3) else 'unresolved'
        w.writerow([R,M,T,c,tR(R),tR(M),tR(T),s1,s2,s3,lab])
print("\nwrote results/ffl_3node_cores_crosscheck.csv")
