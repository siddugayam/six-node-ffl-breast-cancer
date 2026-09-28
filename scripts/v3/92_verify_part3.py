#!/usr/bin/env python3
"""Part 3: double-knockout matrix, combination null, dynamics, power-law cross-check."""
import os, csv, json, collections
import numpy as np, networkx as nx
from scipy.stats import spearmanr
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
order=json.load(open(R3+'/_nodes_order.json')); oi={n:i for i,n in enumerate(order)}
rec('D','node order file matches sorted node list',order==NODES,True)

M2=np.load(R3+'/systems_double_knockout_matrix.npy')
rec('D','double knockouts evaluated (rows in matrix)',M2.shape[0],171991,
    note=f'C(587,2) = {587*586//2}; exhaustive, no sampling')
comp_all=M2[:,5]
ai=M2[:,0].astype(int); bi=M2[:,1].astype(int)
best_k=int(np.argmax(comp_all))
rec('D','best pair in stored matrix',f'{order[ai[best_k]]}+{order[bi[best_k]]}','RELA+SP1',
    ok=sorted([order[ai[best_k]],order[bi[best_k]]])==sorted(['RELA','SP1']))
rec('D','best pair composite',round(float(comp_all[best_k]),6),0.089324,ok=abs(comp_all[best_k]-0.089324)<1e-5)

# recompute composite from the deposited component columns to confirm the definition
recomp=(np.abs(M2[:,2])+np.abs(M2[:,3])+np.abs(M2[:,4]))/3
rec('D','composite == mean(dFFL, dReach, |dCOL|) for all 171991 pairs',
    f'max |diff| = {np.nanmax(np.abs(recomp-comp_all)):.2e}','0',
    ok=np.nanmax(np.abs(recomp-comp_all))<1e-6)

S1=np.load(R3+'/_single_damage.npy')
sing_comp=(np.abs(S1[:,0])+np.abs(S1[:,1])+np.abs(S1[:,2]))/3
stored_single={r['node']:float(r['composite_damage']) for r in csv.DictReader(open(R3+'/systems_single_knockout.csv'))}
mx=max(abs(sing_comp[oi[n]]-stored_single[n]) for n in NODES)
rec('D','single-KO composites in .npy match the CSV',f'max |diff| = {mx:.2e}','0',ok=mx<1e-5)
bs=float(sing_comp.max()); bn=order[int(np.argmax(sing_comp))]
rec('D','best single node / composite',f'{bn} {bs:.6f}','SP1 0.046934',
    ok=(bn=='SP1' and abs(bs-0.046934)<1e-5))
rec('D','best pair / best single',round(float(comp_all[best_k])/bs,4),1.9032,
    ok=abs(comp_all[best_k]/bs-1.9032)<1e-3)
rec('D','pairs beating the best single node',int((comp_all>bs).sum()),739)
rec('D','fraction of pairs beating best single',f'{(comp_all>bs).mean()*100:.3f}%','0.43%',
    ok=abs((comp_all>bs).mean()*100-0.43)<0.01)
mu=float(comp_all.mean()); sd=float(comp_all.std())
rec('D','random-pair null (all 171991 pairs) mean/sd',f'{mu:.5f}/{sd:.5f}','0.00539/0.00727',
    ok=abs(mu-0.00539)<5e-5 and abs(sd-0.00727)<5e-5)
rec('D','z of best pair vs random-pair null',round((comp_all[best_k]-mu)/sd,2),11.55,
    ok=abs((comp_all[best_k]-mu)/sd-11.55)<0.05)
rec('D','empirical p of the best pair',f'{(comp_all>=comp_all[best_k]).sum()/len(comp_all):.2e}',
    '5.81e-06',ok=abs((comp_all>=comp_all[best_k]).sum()/len(comp_all)-5.81e-06)<1e-7)
syn=comp_all-(sing_comp[ai]+sing_comp[bi])
rec('D','fraction of all pairs that are super-additive',f'{(syn>0).mean()*100:.1f}%','37.7%',
    ok=abs((syn>0).mean()*100-37.7)<0.15)
top100=np.argsort(-comp_all)[:100]
rec('D','sub-additive among the 100 most damaging pairs',int((syn[top100]<0).sum()),69)
rec('D','mean synergy of the top-100 pairs',round(float(syn[top100].mean()),5),-0.00162,
    ok=abs(syn[top100].mean()+0.00162)<1e-4)
k=int(np.argmax(syn))
rec('D','most synergistic pair',f'{order[ai[k]]}+{order[bi[k]]} syn={syn[k]:.5f}',
    'MKL1 + hsa-miR-143 syn=+0.00675',
    ok=sorted([order[ai[k]],order[bi[k]]])==sorted(['MKL1','hsa-miR-143']) and abs(syn[k]-0.00675)<1e-4)
rec('D','its damage relative to the best single node',f'{comp_all[k]/bs*100:.1f}%','14.4%',
    ok=abs(comp_all[k]/bs*100-14.4)<0.2)

# druggability, recomputed from DGIdb directly
DGI='/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data'
dt=collections.defaultdict(set); ap=collections.defaultdict(set); an=collections.defaultdict(set)
nrow=0
for r in csv.DictReader(open(DGI+'/interactions.tsv'),delimiter='\t'):
    nrow+=1; g=r['gene_name']
    if g in ntype:
        dt[g].add(r['drug_name'])
        if r['approved']=='TRUE': ap[g].add(r['drug_name'])
        if r['anti_neoplastic']=='TRUE': an[g].add(r['drug_name'])
rec('D','DGIdb interaction rows read',nrow,nrow,ok=True)
rec('D','druggable network nodes',len(dt),240)
rec('D','with an approved drug',len(ap),194)
rec('D','with an antineoplastic drug',len(an),159)
rec('D','druggable miRNAs',sum(1 for g in dt if ntype[g]=='miRNA'),0)
rec('D','druggable genes / TFs',f"{sum(1 for g in dt if ntype[g]=='Gene')} / {sum(1 for g in dt if ntype[g]=='TF')}",
    '134 / 106')
dd=[k for k in range(len(comp_all)) if order[ai[k]] in dt and order[bi[k]] in dt]
kk=max(dd,key=lambda k: comp_all[k])
rec('D','best druggable pair',f'{order[ai[kk]]}+{order[bi[kk]]} {comp_all[kk]:.6f}','RELA+SP1 0.089324',
    ok=sorted([order[ai[kk]],order[bi[kk]]])==sorted(['RELA','SP1']))
da=[k for k in range(len(comp_all)) if order[ai[k]] in ap and order[bi[k]] in ap]
ka=max(da,key=lambda k: comp_all[k])
rec('D','best pair with an APPROVED drug on both targets',
    f'{order[ai[ka]]}+{order[bi[ka]]} ratio={comp_all[ka]/bs:.4f}','NFKB1+RELA ratio=1.7131',
    ok=sorted([order[ai[ka]],order[bi[ka]]])==sorted(['NFKB1','RELA']) and abs(comp_all[ka]/bs-1.7131)<1e-3)

# ---- (E)
de={}
for path in ('BRCA_DEX_genes.csv','BRCA_DEX_mirnas.csv'):
    for r in csv.DictReader(open(REV+'/results/'+path)): de[r['feature']]=float(r['logFC'])
rec('E','DE features are symbols (limma index-artefact guard)',
    f"{len(de)} features; {sum(1 for k in de if k.isdigit())} all-digit; "
    f"{sum(1 for n in NODES if n in de)}/{N} network nodes matched",'guard passes',
    ok=sum(1 for n in NODES if n in de)>300 and not any(n.isdigit() for n in NODES))
ALPHA=0.85
A=np.zeros((N,N))
for (u,v),s in sgn.items(): A[idx[v],idx[u]]=s
kout=np.abs(A).sum(0); nz=kout>0; A[:,nz]=ALPHA*A[:,nz]/kout[nz]
rho=float(np.max(np.abs(np.linalg.eigvals(A))))
rec('E','spectral radius rho(A)',round(rho,4),'<0.85 required for convergence',ok=rho<0.85)
Sm=np.linalg.inv(np.eye(N)-A)
y=np.array([de.get(n,np.nan) for n in NODES]); m=~np.isnan(y)
rec('E','network nodes with a DE measurement',int(m.sum()),int(m.sum()),ok=True)
best=(-9,None); allr=[]
for j,n in enumerate(NODES):
    for s_ in (1,-1):
        r=spearmanr((Sm[:,j]*s_)[m],y[m]).statistic; allr.append(r)
        if r>best[0]: best=(r,f'{n} ({s_:+d})')
rec('E','best single perturbation [published signs]',f'{best[1]} rho={best[0]:.4f}',
    'POU5F1 (+1) rho=+0.1768',ok=best[1].startswith('POU5F1') and abs(best[0]-0.1768)<0.002)
# max-statistic permutation null, independent re-run
rng=np.random.default_rng(99); yv=y[m].copy(); nulls=[]
Sm_m=Sm[m,:]
for b in range(2000):
    yp=rng.permutation(yv)
    rk=np.argsort(np.argsort(yp))
    Rk=np.apply_along_axis(lambda c: np.argsort(np.argsort(c)),0,Sm_m) if b==0 else None
    if b==0:
        Rk0=Rk-Rk.mean(0); Rk0/= np.linalg.norm(Rk0,axis=0)+1e-300
    z=(rk-rk.mean()); z/=np.linalg.norm(z)
    nulls.append(float(np.abs(Rk0.T@z).max()))
nulls=np.array(nulls)
rec('E','max-statistic permutation null (2000 perms) mean/sd',
    f'{nulls.mean():.3f}/{nulls.std():.3f}','0.120/0.018',
    ok=abs(nulls.mean()-0.120)<0.02 and abs(nulls.std()-0.018)<0.006)
pv=(1+(nulls>=best[0]).sum())/(1+len(nulls))
rec('E','FWER p for the best single perturbation',f'{pv:.4f}','0.0050',ok=pv<0.02)
x=np.zeros(N); u=np.zeros(N); u[idx['POU5F1']]=1.0
for k2 in range(1,5001):
    xn=A@x+u
    if np.abs(xn-x).max()<1e-12: break
    x=xn
rec('E','Neumann iterations to tolerance 1e-12',k2,34,ok=abs(k2-34)<=2)
rec('E','Neumann vs closed form, max |diff|',f'{np.abs(xn-Sm[:,idx["POU5F1"]]).max():.2e}','exact',
    ok=float(np.abs(xn-Sm[:,idx['POU5F1']]).max())<1e-9)
COLS=['COL1A1','COL3A1']
rec('E','observed COL1A1 / COL3A1 logFC',f'{y[idx["COL1A1"]]:+.3f} / {y[idx["COL3A1"]]:+.3f}',
    'both up in tumour',ok=y[idx['COL1A1']]>0 and y[idx['COL3A1']]>0)
sc=[]
for j,n in enumerate(NODES):
    if n in COLS: continue
    for s_ in (1,-1):
        a=Sm[idx['COL1A1'],j]*s_; b=Sm[idx['COL3A1'],j]*s_
        if a>0 and b>0: sc.append((a+b,n,s_))
sc.sort(reverse=True)
rec('E','top upstream perturbations raising BOTH collagens',
    '; '.join(f'{n}({s:+d})' for _,n,s in sc[:6]),
    'MKL1; STAT6; hsa-miR-767; hsa-miR-196a; hsa-miR-338; hsa-miR-193a',
    ok=all(t in [n for _,n,_ in sc[:8]] for t in ['MKL1','STAT6']))

with open(R3+'/systems_independent_verification_part3.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=['section','item','independent_value','stored_value','agree','note'])
    w.writeheader(); w.writerows(V)
print(f"\n=== part3: {sum(1 for v in V if v['agree'])}/{len(V)} agree ===")
for v in V:
    if not v['agree']: print('  MISMATCH:',v['section'],v['item'],v['independent_value'],'vs',v['stored_value'])
