#!/usr/bin/env python3
"""PART 2J, control: is the surviving within-CAF TF -> COL1A1 association just CAF-subtype
composition (donors with more myCAF-like cells have more collagen AND more RUNX2)?
Adds the per-donor CAF subtype fractions and re-tests with them partialled out."""
import anndata as ad, pandas as pd, numpy as np
from scipy import stats
import statsmodels.api as sm
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
def msg(*a): print(*a, flush=True)
A=ad.read_h5ad(f"{REV}/cache/v7/hbca_global.h5ad", backed='r'); o=A.obs
caf=o[o.author_cell_type.astype(str).str.startswith(('CAFs','Myofibroblast'))]
tab=pd.crosstab(caf.donor_id.astype(str), caf.author_cell_type.astype(str))
frac=tab.div(tab.sum(axis=1), axis=0)
frac['n_caf']=tab.sum(axis=1)
frac.to_csv(f"{OUT}/sc_hbca_donor_caf_subtype_fractions.csv")
msg("donors with CAFs:", len(frac)); msg(frac.drop(columns='n_caf').mean().round(3).to_string())

D=pd.read_csv(f"{OUT}/sc_hbca_donor_compartment_pseudobulk.csv")
C=D[(D.compartment=='CAF')&(D.n_cells>=50)].copy().set_index('donor')
F=frac.loc[frac.index.intersection(C.index)]
C=C.loc[F.index]
msg(f"donors used: {len(C)}")
# collagen-high myCAF-like fraction = CAFs (LAMP5) + CAFs (COL11A1+) (the two collagen-high subsets)
myc=[c for c in ['CAFs (LAMP5)','CAFs (COL11A1+)'] if c in F.columns]
C['myCAF_frac']=F[myc].sum(axis=1).values
msg(f"myCAF-like fraction: median {C.myCAF_frac.median():.3f}, IQR {C.myCAF_frac.quantile(.25):.3f}-{C.myCAF_frac.quantile(.75):.3f}")
lg=lambda v: np.log2(v+1)
y=lg(C.COL1A1.values)
ct=stats.spearmanr(C.myCAF_frac.values, y)
msg(f"myCAF-like fraction vs CAF COL1A1: rho={ct.statistic:.3f} p={ct.pvalue:.3g}")
rows=[]
TFS=['NFKB1','RELA','SP1','ETS1','MYC','STAT3','HIF1A','TP53','MKL1','SRF','RUNX1','RUNX2','TWIST1',
     'SNAI2','ZEB1','JUN','FOS','EGR1','CEBPB','YAP1','WWTR1','SMAD3','SMAD4','TGFB1','TGFBR2']
study=pd.get_dummies(C.study.astype(str), drop_first=True).astype(float)
for tf in TFS:
    if tf not in C.columns: continue
    x=lg(C[tf].values)
    r0=stats.spearmanr(x,y)
    # partial out study + myCAF fraction, on ranks
    X=np.column_stack([np.ones(len(C)), study.values, stats.rankdata(C.myCAF_frac.values)])
    rx=sm.OLS(stats.rankdata(x), X).fit().resid
    ry=sm.OLS(stats.rankdata(y), X).fit().resid
    r1=stats.pearsonr(rx,ry)
    # partial out study only
    X0=np.column_stack([np.ones(len(C)), study.values])
    rx0=sm.OLS(stats.rankdata(x), X0).fit().resid; ry0=sm.OLS(stats.rankdata(y), X0).fit().resid
    r2=stats.pearsonr(rx0,ry0)
    rows.append(dict(TF=tf, n=len(C), rho_raw=r0.statistic, p_raw=r0.pvalue,
                     rho_study_adj=r2.statistic, p_study_adj=r2.pvalue,
                     rho_study_and_myCAFfrac_adj=r1.statistic, p_study_and_myCAFfrac_adj=r1.pvalue))
R=pd.DataFrame(rows)
R['fdr_full_adj']=stats.false_discovery_control(R.p_study_and_myCAFfrac_adj)
R=R.sort_values('rho_study_and_myCAFfrac_adj', ascending=False)
R.to_csv(f"{OUT}/sc_tf_collagen_within_caf_composition_adjusted.csv", index=False)
msg("\n=== within-CAF TF vs COL1A1, before and after adjusting for CAF-subtype composition ===")
print(R.round(4).to_string(index=False))
msg(f"\nTFs still positive at FDR<0.05 after study + myCAF-fraction adjustment: "
    f"{sorted(R[(R.rho_study_and_myCAFfrac_adj>0)&(R.fdr_full_adj<0.05)].TF.tolist())}")
