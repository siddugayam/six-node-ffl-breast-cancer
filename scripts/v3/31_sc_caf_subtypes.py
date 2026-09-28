#!/usr/bin/env python3
"""PART 2G) CAF subtypes in the Human Breast Cancer Single Cell Atlas stromal compartment.
Which CAF subset carries the collagen programme, and which carries the TFs?
Pseudobulk = summed RAW counts per (cell type) / total UMI * 1e6  (CPM), the same
construction the project already uses for the Wu et al. atlas.
"""
import h5py, numpy as np, pandas as pd, sys, os
from scipy.sparse import csr_matrix

REV="/path/to/revision"
CA=f"{REV}/cache/v7"; OUT=f"{REV}/results/v3"
H5=f"{CA}/hbca_stromal.h5ad"
def msg(*a): print(*a, flush=True)

import anndata as ad
A = ad.read_h5ad(H5, backed='r')
obs = A.obs.copy()
var = A.var.copy()
genes = var['feature_name'].astype(str).values
msg(f"stromal compartment: {A.n_obs} cells x {A.n_vars} genes")
msg("studies:", dict(obs.batch.value_counts()))
msg("donors:", obs.donor_id.nunique())

# duplicate-symbol guard (rule 4 in spirit: never let indices masquerade as symbols)
dupmask = pd.Series(genes).duplicated(keep=False).values
msg(f"duplicated gene symbols in var: {dupmask.sum()} rows "
    f"({pd.Series(genes[dupmask]).nunique()} symbols) -- these are summed per symbol")
assert (pd.Series(genes).str.match(r'^[A-Za-z]').mean() > 0.90), "var feature_name does not look like symbols"

ct = obs['author_cell_type'].astype(str).values
study = obs['batch'].astype(str).values
cats = sorted(set(ct))
cat_idx = {c:i for i,c in enumerate(cats)}
studies = sorted(set(study))
st_idx = {s:i for i,s in enumerate(studies)}

f = h5py.File(H5,'r')
Xg = f['raw/X']
indptr = Xg['indptr'][:]; n_obs = A.n_obs; n_var = A.n_vars
sums   = np.zeros((len(cats), n_var), dtype=np.float64)                 # by cell type
sums_s = np.zeros((len(cats)*len(studies), n_var), dtype=np.float64)    # by cell type x study
ncell  = np.zeros(len(cats), dtype=np.int64)
ncell_s= np.zeros(len(cats)*len(studies), dtype=np.int64)
CH=20000
for start in range(0, n_obs, CH):
    stop=min(start+CH, n_obs)
    lo, hi = indptr[start], indptr[stop]
    data = Xg['data'][lo:hi]; ind = Xg['indices'][lo:hi]
    ip = indptr[start:stop+1]-lo
    M = csr_matrix((data, ind, ip), shape=(stop-start, n_var))
    for c in cats:
        m = (ct[start:stop]==c)
        if m.any():
            sums[cat_idx[c]] += np.asarray(M[m].sum(axis=0)).ravel()
            ncell[cat_idx[c]] += int(m.sum())
    for c in cats:
        for s in studies:
            m = (ct[start:stop]==c) & (study[start:stop]==s)
            if m.any():
                k = cat_idx[c]*len(studies)+st_idx[s]
                sums_s[k] += np.asarray(M[m].sum(axis=0)).ravel()
                ncell_s[k] += int(m.sum())
    msg(f"  chunk {start}-{stop} done")
f.close()

# collapse duplicate symbols
sym = pd.Index(genes)
df = pd.DataFrame(sums.T, index=sym)
df = df.groupby(level=0).sum()
tot = df.sum(axis=0).values
cpm = df.divide(tot, axis=1) * 1e6
cpm.columns = cats
msg("pseudobulk CPM matrix:", cpm.shape, "| total UMI per cell type:", dict(zip(cats, tot.astype(int))))

dfs = pd.DataFrame(sums_s.T, index=sym).groupby(level=0).sum()
tots = dfs.sum(axis=0).values
cols_s = [f"{c}|{s}" for c in cats for s in studies]
cpms = dfs.divide(np.where(tots==0,np.nan,tots), axis=1)*1e6
cpms.columns = cols_s
ncs = pd.Series(ncell_s, index=cols_s)

cpm.round(4).to_csv(f"{OUT}/sc_hbca_stromal_pseudobulk_cpm.csv")
pd.DataFrame({'cell_type':cats,'n_cells':ncell,'total_umi':tot.astype(np.int64)}).to_csv(
    f"{OUT}/sc_hbca_stromal_celltype_counts.csv", index=False)

# ---------------- marker panels ----------------
panels = {
 'myCAF'   : ['ACTA2','TAGLN','MYL9','POSTN','COL11A1','MMP11','INHBA','CTHRC1','FN1','TPM2'],
 'iCAF'    : ['CXCL12','CFD','C3','C7','IL6','CXCL14','HAS1','PLA2G2A','APOD','PTGDS'],
 'apCAF'   : ['CD74','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','SLPI','CD24'],
 'collagen_programme': ['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','COL6A2','COL6A3'],
 'panfibro': ['PDGFRA','PDGFRB','DCN','LUM','FAP','THY1','VIM'],
}
rows=[]
for pname, gl in panels.items():
    have=[g for g in gl if g in cpm.index]
    miss=[g for g in gl if g not in cpm.index]
    if miss: msg(f"  panel {pname}: missing {miss}")
    # score = mean over genes of log2(CPM+1), then z across cell types
    L = np.log2(cpm.loc[have]+1)
    sc = L.mean(axis=0)
    z = (sc - sc.mean())/sc.std(ddof=0)
    for c in cats: rows.append(dict(panel=pname, cell_type=c, n_genes=len(have),
                                    mean_log2cpm=sc[c], z_across_celltypes=z[c]))
P = pd.DataFrame(rows)
P.to_csv(f"{OUT}/sc_hbca_caf_marker_panels.csv", index=False)
msg("\n=== marker panel scores (mean log2 CPM+1) by cell type ===")
print(P.pivot(index='cell_type', columns='panel', values='mean_log2cpm').round(3).to_string())

# ---------------- focus genes by CAF subtype ----------------
focus = ['COL1A1','COL3A1','COL1A2','COL5A1','COL11A1','POSTN','FN1','ACTA2','DCN','LUM','PDGFRA','PDGFRB','FAP',
         'NFKB1','RELA','SP1','ETS1','MYC','TP53','STAT3','HIF1A','TGFB1','TGFBR2','SMAD3','MKL1','SRF','TWIST1',
         'RUNX1','RUNX2','CREB1','JUN','FOS','ESR1','EPCAM','PTPRC','PECAM1']
have=[g for g in focus if g in cpm.index]
F = cpm.loc[have].round(3)
F.to_csv(f"{OUT}/sc_hbca_caf_focus_genes_cpm.csv")
msg("\n=== focus genes, CPM by stromal cell type ===")
print(F.to_string())

# which CAF subtype carries the collagen programme?
cafs=[c for c in cats if c.startswith('CAFs') or c=='Myofibroblast']
coll = np.log2(cpm.loc[[g for g in panels['collagen_programme'] if g in cpm.index]]+1).mean(axis=0)
msg("\ncollagen-programme score (mean log2 CPM+1), ranked:")
print(coll.sort_values(ascending=False).round(3).to_string())
msg("\nCOL1A1 CPM ranked:"); print(cpm.loc['COL1A1'].sort_values(ascending=False).round(1).to_string())
msg("\nCOL3A1 CPM ranked:"); print(cpm.loc['COL3A1'].sort_values(ascending=False).round(1).to_string())

tfs=['NFKB1','RELA','SP1','ETS1','MYC','STAT3','HIF1A','MKL1','SRF','RUNX1','TWIST1']
tfs=[t for t in tfs if t in cpm.index]
msg("\nTF CPM by stromal cell type (which subset carries the TFs):")
print(cpm.loc[tfs].round(2).to_string())

# ---------------- per-study reproducibility of the collagen ranking ----------------
rows=[]
for s in studies:
    sub = {c: cpms.get(f"{c}|{s}") for c in cats}
    for c in cats:
        col = sub[c]
        n = ncs.get(f"{c}|{s}",0)
        if col is None or n < 50 or not np.isfinite(col).any(): continue
        rows.append(dict(study=s, cell_type=c, n_cells=int(n),
                         COL1A1_cpm=float(col.get('COL1A1', np.nan)),
                         COL3A1_cpm=float(col.get('COL3A1', np.nan)),
                         collagen_score=float(np.log2(col[[g for g in panels['collagen_programme'] if g in col.index]]+1).mean()),
                         myCAF_score=float(np.log2(col[[g for g in panels['myCAF'] if g in col.index]]+1).mean()),
                         iCAF_score=float(np.log2(col[[g for g in panels['iCAF'] if g in col.index]]+1).mean())))
S = pd.DataFrame(rows)
S.to_csv(f"{OUT}/sc_hbca_caf_by_study.csv", index=False)
msg("\n=== per-study COL1A1 CPM by CAF subtype (n_cells>=50) ===")
piv = S.pivot(index='cell_type', columns='study', values='COL1A1_cpm')
print(piv.round(0).to_string())
msg("\n=== per-study collagen score ===")
print(S.pivot(index='cell_type', columns='study', values='collagen_score').round(2).to_string())
msg("DONE 31")
