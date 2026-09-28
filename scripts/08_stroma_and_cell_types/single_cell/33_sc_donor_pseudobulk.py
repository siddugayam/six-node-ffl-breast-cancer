#!/usr/bin/env python3
"""PART 2J) Donor x cell-type pseudobulk from the global atlas, so that the
TF -> COL1A1 association can be tested WITHIN the CAF compartment alone."""
import h5py, numpy as np, pandas as pd, anndata as ad
from scipy.sparse import csr_matrix, coo_matrix
REV="/path/to/revision"; CA=f"{REV}/cache/v7"; OUT=f"{REV}/results/v3"
H5=f"{CA}/hbca_global.h5ad"
def msg(*a): print(*a, flush=True)
A=ad.read_h5ad(H5, backed='r'); obs=A.obs
genes=A.var['feature_name'].astype(str).values
ct=obs['cell_type'].astype(str).values; act=obs['author_cell_type'].astype(str).values
fine=np.where(pd.Series(act).str.startswith(('CAFs','Myofib','PCs','VSMC','Endo','Prolif Stromal')).values, act, ct)
comp=np.where(pd.Series(fine).str.startswith(('CAFs','Myofibroblast')).values, 'CAF',
      np.where(fine=='malignant cell','malignant',
      np.where(pd.Series(fine).str.startswith('Endo').values,'endothelial',
      np.where(np.isin(fine,['macrophage','cycling macrophage']),'macrophage',
      np.where(pd.Series(fine).str.contains('T cell').values,'Tcell',
      np.where(np.isin(fine,['PCs (ECM)','VSMC']),'perivascular','other'))))))
don=obs['donor_id'].astype(str).values; st=obs['batch'].astype(str).values
grp=np.array([f"{d}||{c}" for d,c in zip(don,comp)])
gl=sorted(set(grp)); gi={g:i for i,g in enumerate(gl)}; gidx=np.array([gi[g] for g in grp])
msg(f"donor x compartment groups: {len(gl)} (donors {len(set(don))}, compartments {sorted(set(comp))})")
f=h5py.File(H5,'r'); Xg=f['raw/X']; indptr=Xg['indptr'][:]; n_obs=A.n_obs; n_var=A.n_vars
# only accumulate genes we need (keeps the accumulator small)
want=set(['COL1A1','COL3A1','COL1A2','COL5A1','COL5A2','COL6A3','COL11A1','POSTN','FN1','ACTA2','TAGLN',
          'NFKB1','RELA','SP1','ETS1','MYC','STAT3','HIF1A','TP53','MKL1','SRF','TGFB1','TGFB2','TGFB3',
          'TGFBR1','TGFBR2','SMAD2','SMAD3','SMAD4','RUNX1','RUNX2','TWIST1','SNAI1','SNAI2','ZEB1',
          'JUN','FOS','EGR1','CEBPB','YAP1','WWTR1','SERPINE1','LOX','PDGFRA','PDGFRB','DCN','LUM','FAP',
          'EPCAM','KRT18','PTPRC','PECAM1','ESR1','VEGFA','MMP2','MMP11','MMP14','THY1','CXCL12','IL6','CFD'])
keepg=np.array([i for i,g in enumerate(genes) if g in want])
kg=genes[keepg]
msg(f"genes accumulated: {len(keepg)} ({len(set(kg))} unique symbols)")
acc=np.zeros((len(gl), len(keepg))); tot=np.zeros(len(gl)); ncell=np.bincount(gidx, minlength=len(gl))
CH=25000
for s0 in range(0,n_obs,CH):
    s1=min(s0+CH,n_obs); lo,hi=indptr[s0],indptr[s1]
    M=csr_matrix((Xg['data'][lo:hi],Xg['indices'][lo:hi],indptr[s0:s1+1]-lo),shape=(s1-s0,n_var))
    S=coo_matrix((np.ones(s1-s0),(gidx[s0:s1],np.arange(s1-s0))),shape=(len(gl),s1-s0)).tocsr()
    acc += np.asarray((S @ M[:,keepg]).todense())
    tot += np.asarray(S @ np.asarray(M.sum(axis=1))).ravel()
    if (s0//CH)%8==0: msg(f"  {s1}/{n_obs}")
f.close()
D=pd.DataFrame(acc, index=gl, columns=kg).T.groupby(level=0).sum().T
CPM=D.divide(np.where(tot==0,np.nan,tot), axis=0)*1e6
CPM['donor']=[g.split('||')[0] for g in gl]; CPM['compartment']=[g.split('||')[1] for g in gl]
CPM['n_cells']=ncell; CPM['total_umi']=tot.astype(np.int64)
d2s=dict(zip(don,st)); CPM['study']=[d2s[d] for d in CPM.donor]
CPM.to_csv(f"{OUT}/sc_hbca_donor_compartment_pseudobulk.csv", index=False)
msg(f"wrote donor x compartment pseudobulk: {CPM.shape}")
msg(CPM.groupby('compartment').agg(n_donors=('donor','nunique'), n_with_50cells=('n_cells',lambda x:(x>=50).sum())).to_string())
msg("DONE 33")
