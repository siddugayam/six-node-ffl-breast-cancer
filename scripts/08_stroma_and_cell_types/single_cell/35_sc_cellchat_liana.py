#!/usr/bin/env python3
"""PART 2I) Cell-cell communication on the Human Breast Cancer Single Cell Atlas,
focused on the CAF -> cancer-cell axis, collagen/ECM and TGF-beta.
Method: LIANA+ rank_aggregate over five scoring methods, run separately on the
CellChatDB and the LIANA consensus ligand-receptor resources."""
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, liana as li
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
rng=np.random.default_rng(20260909)
def msg(*a): print(*a, flush=True)

A=ad.read_h5ad(f"{REV}/cache/v7/hbca_global.h5ad", backed='r')
o=A.obs
act=o.author_cell_type.astype(str); ct=o.cell_type.astype(str)
lab=np.where(act.str.startswith('CAFs').values, act.values,
     np.where(act.str.startswith('Myofib').values, 'Myofibroblast',
     np.where(ct.values=='malignant cell','Cancer cell',
     np.where(act.str.startswith('Endo').values,'Endothelial',
     np.where(np.isin(ct.values,['macrophage','cycling macrophage']),'Macrophage',
     np.where(ct.str.contains('T cell').values,'T cell',
     np.where(np.isin(act.values,['PCs (ECM)','VSMC']),'Perivascular','Other')))))))
lab=pd.Series(lab, index=o.index)
keep_types=[t for t in lab.unique() if t!='Other']
msg("labels:", lab.value_counts().to_dict())

# subsample to keep the run tractable: at most 12,000 cells per label
idx=[]
for t in keep_types:
    ii=np.where(lab.values==t)[0]
    idx.append(ii if len(ii)<=12000 else rng.choice(ii, 12000, replace=False))
idx=np.sort(np.concatenate(idx))
msg(f"subsampled {len(idx)} cells from {A.n_obs}")
S=A[idx].to_memory()
S.obs['ctype']=pd.Categorical(lab.values[idx])
S.var_names=S.var['feature_name'].astype(str)
S.var_names_make_unique()
msg("subset:", S.shape, "| X is log-normalised (cellxgene convention)")
msg(S.obs.ctype.value_counts().to_string())

res={}
for resource in ['cellchatdb','consensus']:
    msg(f"\n=== LIANA rank_aggregate, resource = {resource} ===")
    li.mt.rank_aggregate(S, groupby='ctype', resource_name=resource, expr_prop=0.1,
                         use_raw=False, verbose=True, n_perms=100, seed=1)
    r=S.uns['liana_res'].copy(); r['resource']=resource
    res[resource]=r
    r.to_csv(f"{OUT}/sc_liana_{resource}_all.csv.gz", index=False, compression='gzip')
    msg(f"  {len(r)} ligand-receptor x source x target rows")

R=pd.concat(res.values())
caf_src=[t for t in keep_types if t.startswith('CAFs') or t=='Myofibroblast']
sub=R[(R.source.isin(caf_src)) & (R.target=='Cancer cell')].copy()
sub=sub.sort_values(['resource','magnitude_rank'])
sub.to_csv(f"{OUT}/sc_liana_CAF_to_cancer.csv", index=False)
msg(f"\nCAF -> cancer-cell interactions: {len(sub)} rows")
for resource,g in sub.groupby('resource'):
    msg(f"\n=== TOP 25 CAF -> cancer cell, {resource} (by magnitude rank) ===")
    cols=[c for c in ['source','target','ligand_complex','receptor_complex','magnitude_rank','specificity_rank','lr_means','cellphone_pvals'] if c in g.columns]
    print(g.drop_duplicates(['ligand_complex','receptor_complex','source']).head(25)[cols].to_string(index=False))

# focused: collagen / ECM and TGF-beta families
coll=sub[sub.ligand_complex.str.startswith(('COL','FN1','LAM','POSTN','THBS','TNC','SPP1','HSPG2','NID'))]
tgfb=sub[sub.ligand_complex.str.startswith('TGFB') | sub.receptor_complex.str.contains('TGFBR|ACVR|SMAD', na=False)]
coll.to_csv(f"{OUT}/sc_liana_CAF_to_cancer_ECM.csv", index=False)
tgfb.to_csv(f"{OUT}/sc_liana_CAF_to_cancer_TGFB.csv", index=False)
msg(f"\n=== CAF -> cancer, ECM/collagen ligands: {len(coll)} rows; top 25 ===")
cols=[c for c in ['resource','source','ligand_complex','receptor_complex','magnitude_rank','specificity_rank','lr_means'] if c in coll.columns]
print(coll.sort_values('magnitude_rank').head(25)[cols].to_string(index=False))
msg(f"\n=== CAF -> cancer, TGF-beta axis: {len(tgfb)} rows; top 15 ===")
print(tgfb.sort_values('magnitude_rank').head(15)[cols].to_string(index=False))

# which CAF subtype sends most of the collagen signal?
if len(coll):
    msg("\n=== ECM/collagen CAF->cancer interactions per CAF subtype (count in the top 5% by magnitude) ===")
    thr=R.magnitude_rank.quantile(0.05)
    print(coll[coll.magnitude_rank<=thr].groupby(['resource','source']).size().to_string())
msg("DONE 35")
