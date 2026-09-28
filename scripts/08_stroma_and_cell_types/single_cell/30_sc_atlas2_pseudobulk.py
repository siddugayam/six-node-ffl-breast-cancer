#!/usr/bin/env python3
"""PART 2F) Second (and much larger) breast atlas: the Human Breast Cancer Single Cell Atlas
global object (621,200 cells, 138 donors, 8 constituent studies including Pal et al. 2021 and
Bassez et al. 2021).  This script streams the RAW counts once and writes pseudobulk CPM per
cell type and per cell type x study, so that the COL1A1/COL3A1 CAF enrichment can be
re-derived independently of Wu et al. 2021 (which is one of the eight studies and is
therefore excluded in the confirmatory analysis)."""
import h5py, numpy as np, pandas as pd, anndata as ad
from scipy.sparse import csr_matrix, coo_matrix

REV="/path/to/revision"
CA=f"{REV}/cache/v7"; OUT=f"{REV}/results/v3"
H5=f"{CA}/hbca_global.h5ad"
def msg(*a): print(*a, flush=True)

A=ad.read_h5ad(H5, backed='r')
obs=A.obs; var=A.var
genes=var['feature_name'].astype(str).values
msg(f"global atlas: {A.n_obs} cells x {A.n_vars} genes; donors {obs.donor_id.nunique()}; studies {obs.batch.nunique()}")
ct=obs['cell_type'].astype(str).values
act=obs['author_cell_type'].astype(str).values
st=obs['batch'].astype(str).values
# stromal cells keep their fine author label (CAF subtypes); everything else uses the ontology label
fine=np.where(pd.Series(act).str.startswith(('CAFs','Myofib','PCs','VSMC','Endo','Prolif Stromal')).values, act, ct)
grp=np.array([f"{a}||{b}" for a,b in zip(fine, st)])
gl=sorted(set(grp)); gi={g:i for i,g in enumerate(gl)}
gidx=np.array([gi[g] for g in grp], dtype=np.int64)
msg(f"groups (cell type x study): {len(gl)}")

f=h5py.File(H5,'r'); Xg=f['raw/X']
indptr=Xg['indptr'][:]; n_obs=A.n_obs; n_var=A.n_vars
msg(f"raw nnz = {indptr[-1]:,}")
acc=np.zeros((len(gl), n_var), dtype=np.float64)
ncell=np.bincount(gidx, minlength=len(gl))
CH=25000
for s0 in range(0, n_obs, CH):
    s1=min(s0+CH, n_obs); lo,hi=indptr[s0],indptr[s1]
    M=csr_matrix((Xg['data'][lo:hi], Xg['indices'][lo:hi], indptr[s0:s1+1]-lo), shape=(s1-s0, n_var))
    sub=gidx[s0:s1]
    S=coo_matrix((np.ones(s1-s0), (sub, np.arange(s1-s0))), shape=(len(gl), s1-s0)).tocsr()
    acc += np.asarray((S @ M).todense())
    if (s0//CH) % 5 == 0: msg(f"  {s1}/{n_obs}")
f.close()

sym=pd.Index(genes)
D=pd.DataFrame(acc.T, index=sym).groupby(level=0).sum()
tot=D.sum(axis=0).values
CPMs=D.divide(np.where(tot==0,np.nan,tot), axis=1)*1e6
CPMs.columns=gl
CPMs.round(4).to_csv(f"{OUT}/sc_hbca_global_pseudobulk_by_celltype_study.csv")
pd.DataFrame({'group':gl,'n_cells':ncell,'total_umi':tot.astype(np.int64)}).to_csv(
    f"{OUT}/sc_hbca_global_group_counts.csv", index=False)
msg(f"wrote pseudobulk: {CPMs.shape}")

# collapse to cell type (all studies pooled)
ctl=[g.split('||')[0] for g in gl]
Dt=D.copy(); Dt.columns=ctl; Dt=Dt.T.groupby(level=0).sum().T
tot2=Dt.sum(axis=0).values
CPM=Dt.divide(tot2, axis=1)*1e6
CPM.round(4).to_csv(f"{OUT}/sc_hbca_global_pseudobulk_by_celltype.csv")
pd.DataFrame({'cell_type':Dt.columns,'total_umi':tot2.astype(np.int64),
              'n_cells':[int(ncell[[i for i,g in enumerate(gl) if g.split('||')[0]==c]].sum()) for c in Dt.columns]}
             ).to_csv(f"{OUT}/sc_hbca_global_celltype_counts.csv", index=False)
msg("wrote per-cell-type pseudobulk", CPM.shape)
msg("cell types:", list(CPM.columns))
msg("DONE 30")
