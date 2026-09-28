#!/usr/bin/env python3
"""PART 2H, independent replication in the Wu et al. 2021 atlas (the project's original atlas),
using its own celltype_minor pseudobulk, to check that the compartment segregation of the two
target sets is not specific to the Human Breast Cancer Single Cell Atlas."""
import pandas as pd, numpy as np
from scipy import stats
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
rng=np.random.default_rng(20260909)
def msg(*a): print(*a, flush=True)
C=pd.read_csv(f"{OUT}/sc_wu_pseudobulk_celltype_minor_cpm.csv", index_col=0)
C=C[C.sum(axis=1)>0]
msg("Wu pseudobulk:", C.shape)
FIB=[c for c in C.columns if c.startswith('CAFs')]
MAL=[c for c in C.columns if c.startswith('Cancer_')]
msg("CAF columns:", FIB); msg("cancer columns:", MAL)
share=C.div(C.sum(axis=1), axis=0)
fib=share[FIB].sum(axis=1); mal=share[MAL].sum(axis=1)
sets={
 'mir29_anticorrelated': pd.read_csv(f"{OUT}/screens_mir29_anticorrelated_genes.csv").gene.tolist(),
 'mir130a_anticorrelated': pd.read_csv(f"{OUT}/mir130a_anticorrelated_validated_genes.csv").gene.tolist(),
 'mir29_strong': pd.read_csv(f"{OUT}/screens_mir29_target_set.csv").query("tier=='STRONG_lowthroughput'").gene.tolist(),
 'mir130a_strong': pd.read_csv(f"{OUT}/mir130a_target_set.csv").query("tier=='STRONG_lowthroughput'").gene.tolist()}
tot=C.sum(axis=1); dec=pd.qcut(tot.rank(method='first'), 10, labels=False)
B=2000; rows=[]
for name,gl in sets.items():
    g=[x for x in dict.fromkeys(gl) if x in C.index]
    pool={d:C.index[(dec==d)&(~C.index.isin(set(gl)))].values for d in range(10)}
    gd=dec.loc[g].values
    of=float(fib.loc[g].mean()); om=float(mal.loc[g].mean())
    nf=np.empty(B); nm=np.empty(B)
    for b in range(B):
        pick=[rng.choice(pool[d]) for d in gd]
        nf[b]=fib.loc[pick].mean(); nm[b]=mal.loc[pick].mean()
    rows.append(dict(set=name, n_detected=len(g),
        CAF_share_obs=of, CAF_share_null=nf.mean(), CAF_z=(of-nf.mean())/nf.std(ddof=1),
        CAF_emp_p=(1+np.sum(nf>=of))/(1+B),
        cancer_share_obs=om, cancer_share_null=nm.mean(), cancer_z=(om-nm.mean())/nm.std(ddof=1),
        cancer_emp_p=(1+np.sum(nm>=om))/(1+B)))
R=pd.DataFrame(rows); R.to_csv(f"{OUT}/sc_wu_target_set_compartment_summary.csv", index=False)
msg("\n=== Wu et al. 2021: CAF vs cancer-epithelial expression share of the two target sets ===")
print(R.round(4).to_string(index=False))
a=[x for x in set(sets['mir29_anticorrelated']) if x in C.index]
b=[x for x in set(sets['mir130a_anticorrelated']) if x in C.index]
u1=stats.mannwhitneyu(fib.loc[a], fib.loc[b]); u2=stats.mannwhitneyu(mal.loc[a], mal.loc[b])
msg(f"\nhead-to-head (Wu): CAF share median {fib.loc[a].median():.4f} vs {fib.loc[b].median():.4f}, MWU p={u1.pvalue:.3g}")
msg(f"                   cancer share median {mal.loc[a].median():.4f} vs {mal.loc[b].median():.4f}, MWU p={u2.pvalue:.3g}")
