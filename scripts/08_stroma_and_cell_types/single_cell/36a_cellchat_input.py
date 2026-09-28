#!/usr/bin/env python3
"""Export a subsampled, log-normalised expression matrix + cell labels for CellChat."""
import anndata as ad, numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, os
REV="/path/to/revision"; CA=f"{REV}/cache/v7"
rng=np.random.default_rng(20260909)
A=ad.read_h5ad(f"{CA}/hbca_global.h5ad", backed='r'); o=A.obs
act=o.author_cell_type.astype(str); ct=o.cell_type.astype(str)
lab=np.where(act.str.startswith('CAFs').values, act.values,
     np.where(act.str.startswith('Myofib').values,'Myofibroblast',
     np.where(ct.values=='malignant cell','Cancer cell',
     np.where(act.str.startswith('Endo').values,'Endothelial',
     np.where(np.isin(ct.values,['macrophage','cycling macrophage']),'Macrophage',
     np.where(ct.str.contains('T cell').values,'T cell',
     np.where(np.isin(act.values,['PCs (ECM)','VSMC']),'Perivascular','Other')))))))
idx=[]
for t in [x for x in set(lab) if x!='Other']:
    ii=np.where(lab==t)[0]
    idx.append(ii if len(ii)<=5000 else rng.choice(ii,5000,replace=False))
idx=np.sort(np.concatenate(idx))
S=A[idx].to_memory()
S.var_names=S.var['feature_name'].astype(str); S.var_names_make_unique()
X=sp.csc_matrix(S.X.T)     # genes x cells, log-normalised
print("export:", X.shape, "labels:", pd.Series(lab[idx]).value_counts().to_dict(), flush=True)
os.makedirs(f"{CA}/cc", exist_ok=True)
sio.mmwrite(f"{CA}/cc/X.mtx", X)
pd.Series(S.var_names).to_csv(f"{CA}/cc/genes.csv", index=False, header=False)
pd.Series(S.obs_names).to_csv(f"{CA}/cc/cells.csv", index=False, header=False)
pd.Series(lab[idx]).to_csv(f"{CA}/cc/labels.csv", index=False, header=False)
