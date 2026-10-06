#!/usr/bin/env python3
"""PART 2H controls.  (1) Is the CAF-weighting of the miR-29 target set just the collagen
genes themselves?  (2) Where do the network's TFs and FFL hubs sit?  (3) How much of the
apparent CAF-weighting is circular, i.e. inherited from the way the anti-correlated set was
derived in bulk (miR-29 tracks stroma negatively)?"""
import numpy as np, pandas as pd
from scipy import stats
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
rng=np.random.default_rng(20260909)
def msg(*a): print(*a, flush=True)
CPM=pd.read_csv(f"{OUT}/sc_hbca_global_pseudobulk_by_celltype.csv", index_col=0)
cnt=pd.read_csv(f"{OUT}/sc_hbca_global_celltype_counts.csv")
keep=cnt[cnt.n_cells>=500].cell_type.tolist()
C=CPM[keep]; C=C[C.sum(axis=1)>0]
share=C.div(C.sum(axis=1), axis=0)
FIB=[c for c in keep if c.startswith('CAFs') or c=='Myofibroblast']
MAL=['malignant cell']
fibshare=share[FIB].sum(axis=1); malshare=share[MAL].sum(axis=1)
msg(f"genes {C.shape[0]}, fibroblast columns {FIB}")

sets={}
sets['mir29_anticorrelated']=pd.read_csv(f"{OUT}/screens_mir29_anticorrelated_genes.csv").gene.tolist()
a29=pd.read_csv(f"{OUT}/screens_mir29_target_set.csv")
sets['mir29_strong']=a29.query("tier=='STRONG_lowthroughput'").gene.tolist()
sets['mir29_targetscan_conserved']=a29.query("in_targetscan==True").gene.tolist()
sets['mir29_anticorr_noCOL']=[g for g in sets['mir29_anticorrelated'] if not g.startswith('COL')]
nodes=pd.read_csv(f"{REV}/data/canonical_nodes.tsv", sep="\t")
ctrl=pd.read_csv(f"{OUT}/systems_controllability_nodes.csv")
sets['network_TFs']=nodes.query("type=='TF'").name.tolist()
sets['network_genes']=nodes.query("type=='Gene'").name.tolist()
sets['FFL_hubs']=ctrl.query("is_FFL_hub==True").node.tolist()

tot=C.sum(axis=1); dec=pd.qcut(tot.rank(method='first'), 10, labels=False)
B=2000; rows=[]
for name,gl in sets.items():
    g=[x for x in dict.fromkeys(gl) if x in C.index]
    if len(g)<10: msg(f"skip {name} (n={len(g)})"); continue
    pool={d:C.index[(dec==d)&(~C.index.isin(set(gl)))].values for d in range(10)}
    gd=dec.loc[g].values
    of=float(fibshare.loc[g].mean()); om=float(malshare.loc[g].mean())
    nf=np.empty(B); nm=np.empty(B)
    for b in range(B):
        pick=[rng.choice(pool[d]) for d in gd]
        nf[b]=fibshare.loc[pick].mean(); nm[b]=malshare.loc[pick].mean()
    rows.append(dict(set=name, n_in_set=len(gl), n_detected=len(g),
        fibro_share_obs=of, fibro_share_null=nf.mean(), fibro_z=(of-nf.mean())/nf.std(ddof=1),
        fibro_emp_p=(1+np.sum(nf>=of))/(1+B),
        malig_share_obs=om, malig_share_null=nm.mean(), malig_z=(om-nm.mean())/nm.std(ddof=1),
        malig_emp_p=(1+np.sum(nm>=om))/(1+B)))
R=pd.DataFrame(rows)
R.to_csv(f"{OUT}/sc_target_set_compartment_summary.csv", index=False)
msg("\n=== fibroblast vs malignant expression share, gene sets vs expression-decile-matched null ===")
print(R.round(4).to_string(index=False))
msg("\nInterpretation guard: the anti-correlated sets were derived from BULK correlation with a")
msg("miRNA whose own level tracks stromal content, so their compartment bias is partly circular.")
msg("The database-derived sets (STRONG, TargetScan-conserved) carry no such circularity.")
