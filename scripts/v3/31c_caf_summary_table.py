#!/usr/bin/env python3
"""PART 2G summary: assign each CAF subtype to myCAF / iCAF / apCAF by marker score, and
state quantitatively which subtype carries the collagen programme and which carries the TFs."""
import pandas as pd, numpy as np
from scipy import stats
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
def msg(*a): print(*a, flush=True)
P=pd.read_csv(f"{OUT}/sc_hbca_caf_marker_panels.csv")
F=pd.read_csv(f"{OUT}/sc_hbca_caf_focus_genes_cpm.csv", index_col=0)
cnt=pd.read_csv(f"{OUT}/sc_hbca_stromal_celltype_counts.csv")
piv=P.pivot(index='cell_type', columns='panel', values='mean_log2cpm')
cafs=[c for c in piv.index if c.startswith('CAFs') or c=='Myofibroblast']
T=piv.loc[cafs, ['myCAF','iCAF','apCAF','collagen_programme','panfibro']].copy()
T['COL1A1_cpm']=F.loc['COL1A1', cafs].values
T['COL3A1_cpm']=F.loc['COL3A1', cafs].values
T['n_cells']=cnt.set_index('cell_type').loc[cafs,'n_cells'].values
T['myCAF_minus_iCAF']=T.myCAF-T.iCAF
T['assignment']=np.where(T.myCAF_minus_iCAF>0.5,'myCAF-like',
                 np.where(T.myCAF_minus_iCAF<-0.5,'iCAF-like','intermediate'))
T=T.sort_values('COL1A1_cpm', ascending=False)
T.round(3).to_csv(f"{OUT}/sc_caf_subtype_assignment.csv")
msg("=== CAF subtypes: marker assignment and collagen output (HBCA stromal compartment) ===")
print(T.round(2).to_string())
r=stats.spearmanr(T.myCAF_minus_iCAF, np.log2(T.COL1A1_cpm))
msg(f"\nacross the {len(T)} CAF subtypes: myCAF-minus-iCAF score vs log2 COL1A1 CPM: Spearman rho={r.statistic:.3f} p={r.pvalue:.3g}")
my=T[T.assignment=='myCAF-like']; ic=T[T.assignment=='iCAF-like']
if len(my) and len(ic):
    msg(f"myCAF-like subtypes {list(my.index)}: COL1A1 {my.COL1A1_cpm.mean():.0f} CPM, COL3A1 {my.COL3A1_cpm.mean():.0f}")
    msg(f"iCAF-like subtypes  {list(ic.index)}: COL1A1 {ic.COL1A1_cpm.mean():.0f} CPM, COL3A1 {ic.COL3A1_cpm.mean():.0f}")
    msg(f"ratio myCAF-like / iCAF-like: COL1A1 {my.COL1A1_cpm.mean()/ic.COL1A1_cpm.mean():.2f}x, "
        f"COL3A1 {my.COL3A1_cpm.mean()/ic.COL3A1_cpm.mean():.2f}x")
TF=['NFKB1','RELA','SP1','ETS1','MYC','STAT3','HIF1A','SRF','RUNX1','RUNX2','TWIST1','TGFB1','TGFBR2']
TF=[t for t in TF if t in F.index]
Z=F.loc[TF, cafs]
msg("\n=== TF CPM across the CAF subtypes, ordered by that subtype's COL1A1 output ===")
print(Z[T.index].round(2).to_string())
rows=[]
for t in TF:
    rr=stats.spearmanr(np.log2(F.loc[t,T.index].values+1), np.log2(T.COL1A1_cpm.values))
    rows.append(dict(TF=t, rho_vs_COL1A1_across_CAF_subtypes=rr.statistic, p=rr.pvalue, n_subtypes=len(T)))
R=pd.DataFrame(rows).sort_values('rho_vs_COL1A1_across_CAF_subtypes', ascending=False)
R.to_csv(f"{OUT}/sc_caf_subtype_TF_vs_collagen.csv", index=False)
msg("\n=== TF expression vs COL1A1 across CAF subtypes (n=6 subtypes; descriptive, not powered) ===")
print(R.round(3).to_string(index=False))
