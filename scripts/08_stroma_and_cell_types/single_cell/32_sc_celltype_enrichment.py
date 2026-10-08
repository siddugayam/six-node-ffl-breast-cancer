#!/usr/bin/env python3
"""PART 2F (confirmation) and 2H (in which cell types are the miR-29 target sets expressed?).
Uses the Human Breast Cancer Single Cell Atlas global pseudobulk written by script 30.
F: COL1A1/COL3A1 fibroblast-vs-carcinoma enrichment, computed per constituent study, with
   wu_natgen_2021 (the atlas the project already used) removed for the independent estimate.
H: per-cell-type module enrichment of the miR-29 anti-correlated and STRONG (low-throughput)
   target sets, against an expression-matched random-gene null.
"""
import numpy as np, pandas as pd
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
rng=np.random.default_rng(20260909)
def msg(*a): print(*a, flush=True)

CPM=pd.read_csv(f"{OUT}/sc_hbca_global_pseudobulk_by_celltype.csv", index_col=0)
CPMS=pd.read_csv(f"{OUT}/sc_hbca_global_pseudobulk_by_celltype_study.csv", index_col=0)
cnt=pd.read_csv(f"{OUT}/sc_hbca_global_celltype_counts.csv")
grp=pd.read_csv(f"{OUT}/sc_hbca_global_group_counts.csv")
msg("cell types:", CPM.shape[1], "| genes:", CPM.shape[0])
msg(cnt.sort_values('n_cells',ascending=False).head(40).to_string(index=False))

FIB=[c for c in CPM.columns if c.startswith('CAFs') or c in ('fibroblast','myofibroblast cell','Myofibroblast')]
MAL=[c for c in CPM.columns if c=='malignant cell']
msg(f"\nfibroblast/CAF columns: {FIB}\nmalignant columns: {MAL}")

# ---------------- F) COL1A1 / COL3A1 enrichment, per study ----------------
rows=[]
studies=sorted({g.split('||')[1] for g in CPMS.columns})
for g in ['COL1A1','COL3A1','COL1A2','COL5A1','DCN','LUM','PDGFRA','EPCAM','KRT18','PTPRC']:
    if g not in CPMS.index: continue
    for s in studies:
        fc=[f"{c}||{s}" for c in FIB if f"{c}||{s}" in CPMS.columns]
        mc=[f"{c}||{s}" for c in MAL if f"{c}||{s}" in CPMS.columns]
        if not fc or not mc: continue
        nf=grp.set_index('group').loc[fc,'n_cells'].sum(); nm=grp.set_index('group').loc[mc,'n_cells'].sum()
        if nf<100 or nm<100: continue
        # UMI-weighted pooled CPM within each compartment
        wf=grp.set_index('group').loc[fc,'total_umi'].values; wm=grp.set_index('group').loc[mc,'total_umi'].values
        vf=float(np.nansum(CPMS.loc[g,fc].values*wf)/wf.sum())
        vm=float(np.nansum(CPMS.loc[g,mc].values*wm)/wm.sum())
        rows.append(dict(gene=g, study=s, n_fibro_cells=int(nf), n_malignant_cells=int(nm),
                         cpm_fibroblast=vf, cpm_malignant=vm,
                         fold_fibro_over_malignant=(vf+0.01)/(vm+0.01)))
E=pd.DataFrame(rows); E.to_csv(f"{OUT}/sc_hbca_collagen_caf_enrichment_by_study.csv", index=False)
msg("\n=== COL1A1 / COL3A1 fibroblast-over-malignant fold enrichment, per constituent study ===")
piv=E[E.gene.isin(['COL1A1','COL3A1'])].pivot(index='study',columns='gene',values='fold_fibro_over_malignant')
print(piv.round(1).to_string())
ind=E[(E.study!='wu_natgen_2021')&(E.gene.isin(['COL1A1','COL3A1']))]
msg("\nEXCLUDING wu_natgen_2021 (independent studies only):")
print(ind.groupby('gene').fold_fibro_over_malignant.describe()[['count','min','50%','max']].round(1).to_string())
for g in ['COL1A1','COL3A1']:
    sub=E[(E.gene==g)]
    pooledf=np.average(sub[sub.study!='wu_natgen_2021'].cpm_fibroblast, weights=sub[sub.study!='wu_natgen_2021'].n_fibro_cells)
    pooledm=np.average(sub[sub.study!='wu_natgen_2021'].cpm_malignant, weights=sub[sub.study!='wu_natgen_2021'].n_malignant_cells)
    msg(f"  {g}: pooled (Wu excluded) fibroblast {pooledf:.1f} CPM vs malignant {pooledm:.1f} CPM -> {pooledf/max(pooledm,1e-9):.1f}x")

# ---------------- H) in which cell types are the miR-29 target sets expressed? ----------------
sets={}
sets['mir29_anticorrelated']=pd.read_csv(f"{OUT}/screens_mir29_anticorrelated_genes.csv").gene.tolist()
sets['mir29_strong']=pd.read_csv(f"{OUT}/screens_mir29_target_set.csv").query("tier=='STRONG_lowthroughput'").gene.tolist()

# restrict to cell types with enough cells, and to genes detected in the atlas
keep_ct=cnt[cnt.n_cells>=500].cell_type.tolist()
C=CPM[keep_ct]
det=C.sum(axis=1)
C=C[det>0]
msg(f"\ncell types kept (>=500 cells): {len(keep_ct)} | genes detected: {C.shape[0]}")
L=np.log2(C+1)
share=C.div(C.sum(axis=1), axis=0)                # fraction of a gene's summed CPM in each cell type
dom=share.idxmax(axis=1)
domfrac=share.max(axis=1)

# expression-decile-matched null over detected genes
tot=C.sum(axis=1)
dec=pd.qcut(tot.rank(method='first'), 10, labels=False)
rows=[]; B=2000
for name, gl in sets.items():
    g=[x for x in dict.fromkeys(gl) if x in C.index]
    msg(f"set {name}: {len(gl)} genes -> {len(g)} detected in atlas")
    obs=share.loc[g].mean(axis=0)
    # matched null
    pool_by_dec={d:C.index[(dec==d)&(~C.index.isin(set(gl)))].values for d in range(10)}
    gd=dec.loc[g].values
    nullmat=np.zeros((B, len(keep_ct)))
    for b in range(B):
        pick=[rng.choice(pool_by_dec[d]) for d in gd]
        nullmat[b]=share.loc[pick].mean(axis=0).values
    for i,c in enumerate(keep_ct):
        nl=nullmat[:,i]
        rows.append(dict(set=name, n_genes=len(g), cell_type=c,
                         obs_mean_share=float(obs[c]), null_mean_share=float(nl.mean()),
                         null_sd=float(nl.std(ddof=1)),
                         z=float((obs[c]-nl.mean())/max(nl.std(ddof=1),1e-12)),
                         emp_p_enriched=float((1+np.sum(nl>=obs[c]))/(1+B)),
                         n_genes_dominant_here=int((dom.loc[g]==c).sum()),
                         frac_genes_dominant_here=float((dom.loc[g]==c).mean())))
M=pd.DataFrame(rows); M.to_csv(f"{OUT}/sc_target_set_celltype_enrichment.csv", index=False)
msg("\n=== per-cell-type mean expression share, target sets vs expression-matched null (z) ===")
print(M.pivot(index='cell_type',columns='set',values='z').round(2).to_string())
msg("\n=== fraction of set genes whose dominant cell type is that one ===")
print(M.pivot(index='cell_type',columns='set',values='frac_genes_dominant_here').round(3).to_string())
msg("DONE 32")
