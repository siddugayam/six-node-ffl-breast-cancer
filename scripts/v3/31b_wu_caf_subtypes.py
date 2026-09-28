#!/usr/bin/env python3
"""PART 2G, cross-check in the atlas the project already uses (Wu et al. 2021, GSE176078),
which carries its own myCAF-like / MSC iCAF-like CAF labels."""
import pandas as pd, numpy as np
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
WU=f"{REV}/cache/celltype/Wu_etal_2021_BRCA_scRNASeq"; CA=f"{REV}/cache/v7/wu"
def msg(*a): print(*a, flush=True)
genes=[l.strip() for l in open(f"{WU}/count_matrix_genes.tsv")]
msg(f"Wu genes: {len(genes)} | duplicated symbols: {pd.Series(genes).duplicated().sum()}")
pb=pd.read_csv(f"{CA}/pb_minor.tsv", sep="\t", header=None, names=['gi','ct','umi'])
pb['gene']=[genes[i-1] for i in pb.gi]
M=pb.pivot_table(index='gene', columns='ct', values='umi', aggfunc='sum').fillna(0)
tot=M.sum(axis=0)
CPM=M.divide(tot, axis=1)*1e6
msg(f"Wu pseudobulk (celltype_minor): {CPM.shape}")
cells=pd.read_csv(f"{CA}/cellmap.tsv", sep="\t")
nc=cells.minor.value_counts()
msg("\ntotal UMI and cell counts per minor type:")
print(pd.DataFrame({'n_cells':nc, 'total_umi':tot.astype(np.int64)}).dropna().sort_values('n_cells',ascending=False).to_string())
CPM.round(3).to_csv(f"{OUT}/sc_wu_pseudobulk_celltype_minor_cpm.csv")

caf=[c for c in CPM.columns if c.startswith('CAFs')]
msg(f"\nWu CAF subsets: {caf}")
panels={'myCAF':['ACTA2','TAGLN','MYL9','POSTN','COL11A1','MMP11','INHBA','CTHRC1','FN1','TPM2'],
        'iCAF':['CXCL12','CFD','C3','C7','IL6','CXCL14','HAS1','PLA2G2A','APOD','PTGDS'],
        'apCAF':['CD74','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','SLPI'],
        'collagen':['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','COL6A2','COL6A3']}
rows=[]
for p,gl in panels.items():
    have=[g for g in gl if g in CPM.index]
    miss=[g for g in gl if g not in CPM.index]
    if miss: msg(f"  panel {p} missing: {miss}")
    sc=np.log2(CPM.loc[have]+1).mean(axis=0)
    for c in CPM.columns: rows.append(dict(panel=p, cell_type=c, n_genes=len(have), score=sc[c]))
P=pd.DataFrame(rows); P.to_csv(f"{OUT}/sc_wu_caf_marker_panels.csv", index=False)
msg("\n=== Wu et al.: marker-panel scores by CAF subset (and reference types) ===")
sel=caf+['Cancer_LumA_SC','Cancer_Basal_SC','Endothelial_ACKR1','PVL_Differentiated','Macrophage','Myoepithelial']
sel=[s for s in sel if s in CPM.columns]
print(P[P.cell_type.isin(sel)].pivot(index='cell_type',columns='panel',values='score').round(3).to_string())
focus=['COL1A1','COL3A1','COL1A2','COL11A1','POSTN','ACTA2','DCN','LUM','PDGFRA','FAP','CXCL12','CFD',
       'NFKB1','RELA','SP1','ETS1','MYC','STAT3','HIF1A','MKL1','SRF','RUNX1','RUNX2','TWIST1','TGFB1','TGFBR2']
have=[g for g in focus if g in CPM.index]
msg("\n=== Wu et al.: focus genes (CPM) across CAF subsets and reference types ===")
print(CPM.loc[have, sel].round(2).to_string())
CPM.loc[have].round(3).to_csv(f"{OUT}/sc_wu_caf_focus_genes_cpm.csv")

msg("\nCOL1A1 CPM in the two Wu CAF subsets:")
for c in caf: msg(f"   {c}: COL1A1 {CPM.loc['COL1A1',c]:.0f} | COL3A1 {CPM.loc['COL3A1',c]:.0f} | n_cells {int(nc[c])}")
if 'CAFs_myCAF-like' in CPM.columns and 'CAFs_MSC_iCAF-like' in CPM.columns:
    for g in ['COL1A1','COL3A1','COL1A2','COL11A1','POSTN','NFKB1','RELA','SP1','ETS1','MYC','RUNX1','TWIST1']:
        if g in CPM.index:
            a=CPM.loc[g,'CAFs_myCAF-like']; b=CPM.loc[g,'CAFs_MSC_iCAF-like']
            msg(f"   {g:8s} myCAF-like {a:10.2f} | iCAF-like {b:10.2f} | ratio {(a+0.01)/(b+0.01):6.2f}")
msg("DONE 31b")
