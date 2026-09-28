#!/usr/bin/env python3
"""Independent re-derivation of the DEMETER2 and Project SCORE set-level decile-matched null tests."""
import numpy as np, pandas as pd, re, sys
REV="/path/to/revision"; OUT=f"{REV}/results/v3"; CA=f"{REV}/cache/v7"
rng=np.random.default_rng(777)

d2b=pd.read_csv("/path/to/scratch/d2_breast_means.csv",index_col=0).iloc[:,0]
A130=pd.read_csv(f"{OUT}/mir130a_target_set.csv"); R130=pd.read_csv(f"{OUT}/mir130a_tcga_target_correlations.csv")
A29=pd.read_csv(f"{OUT}/screens_mir29_target_set.csv"); R29=pd.read_csv(f"{OUT}/screens_mir29_tcga_correlations.csv")
edges=pd.read_csv(f"{REV}/data/canonical_edges.tsv",sep="\t"); nodes=pd.read_csv(f"{REV}/data/canonical_nodes.tsv",sep="\t")
ctrl=pd.read_csv(f"{OUT}/systems_controllability_nodes.csv"); hubs=list(ctrl.loc[ctrl.is_FFL_hub==True,"node"])
anti29=pd.read_csv(f"{OUT}/screens_mir29_anticorrelated_genes.csv")["gene"]
print("A130 tiers:",A130.tier.value_counts().to_dict())
print("A29 tiers:",A29.tier.value_counts().to_dict())
print("R29 columns:",list(R29.columns)[:12])

sets={
 "mir130a_STRONG":A130.loc[A130.tier=="STRONG_lowthroughput","gene"],
 "mir130a_VALIDATED_anticorr":R130.loc[R130.fdr_VALIDATED_any.notna()&(R130.fdr_VALIDATED_any<0.05)&(R130.rho_3p<0),"gene"],
 "mir130a_TARGETSCAN":A130.loc[A130.in_targetscan==True,"gene"],
 "mir130a_AUTHOR_NETWORK":edges.loc[(edges.source=="hsa-miR-130a")&(edges.edge_type=="miRNA_target"),"target"],
 "mir29_STRONG":A29.loc[A29.tier=="STRONG_lowthroughput","gene"],
 "mir29_anticorrelated":anti29,
 "mir29_TARGETSCAN":A29.loc[A29.in_targetscan==True,"gene"],
 "FFL_hubs":pd.Series(hubs),
 "network_TFs":nodes.loc[nodes.type=="TF","name"],
 "network_genes":nodes.loc[nodes.type=="Gene","name"]}

def run(score, label, B=2000, direction="less"):
    gs=set(score.dropna().index); universe=set(R130.gene)
    mexp=pd.Series(R130.mean_expr.values,index=R130.gene)
    pooluniv=sorted(gs&universe)
    q=np.quantile(mexp[pooluniv],np.arange(0,1.01,0.1)); q=np.unique(q)
    dec=pd.Series(np.clip(np.searchsorted(q,mexp[pooluniv],side="left")-1,0,len(q)-2),index=pooluniv)
    pools={d:np.array([g for g in pooluniv if dec[g]==d]) for d in dec.unique()}
    rows=[]
    for nm,gser in sets.items():
        g=sorted(set(gser.dropna())&gs&universe)
        if len(g)<5: continue
        excl=set(A130.gene) if nm.startswith("mir130a") else (set(A29.gene) if nm.startswith("mir29") else set(g))
        pb={d:np.array([x for x in pools[d] if x not in excl]) for d in pools}
        obs=float(np.nanmean(score[g]))
        nullm=np.empty(B)
        idx=[dec[x] for x in g]
        draws=np.empty((len(g),B),dtype=object)
        for i,d in enumerate(idx):
            draws[i,:]=rng.choice(pb[d],size=B,replace=True)
        vals=np.vectorize(lambda x: score.get(x,np.nan))(draws)
        nullm=np.nanmean(vals,axis=0)
        p_more=(1+np.sum(nullm<=obs))/(1+B)
        p_less=(1+np.sum(nullm>=obs))/(1+B)
        rows.append(dict(resource=label,set=nm,n=len(g),obs_mean=obs,null_mean=float(nullm.mean()),
                         null_sd=float(nullm.std(ddof=1)),emp_p_more_essential=p_more,emp_p_less_essential=p_less))
    return pd.DataFrame(rows)

r1=run(d2b,"DEMETER2")
print("\n=== DEMETER2 set-level (independent) ===")
print(r1.to_string(index=False,float_format=lambda x:"%.4f"%x))
dep=pd.read_csv(f"{OUT}/screens_demeter2_setlevel.csv")
m=r1.merge(dep[["set","n","obs_mean_D2_breast","null_mean_D2","emp_p_more_essential"]],on="set",suffixes=("_new","_dep"))
m["d_obs"]=m.obs_mean-m.obs_mean_D2_breast; m["d_null"]=m.null_mean-m.null_mean_D2
print("\n--- agreement with deposited ---")
print(m[["set","n_new","n_dep","obs_mean","obs_mean_D2_breast","d_obs","null_mean","null_mean_D2","d_null","emp_p_more_essential_new","emp_p_more_essential_dep"]].to_string(index=False,float_format=lambda x:"%.5f"%x))
r1.to_csv(f"{OUT}/screens_verify_demeter2_setlevel.csv",index=False)
