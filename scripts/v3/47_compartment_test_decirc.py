#!/usr/bin/env python3
"""Compartment test (Part 2H) run on RAW vs CAF-ADJUSTED anti-correlated target sets,
plus a formal two-set contrast for the database-derived (non-circular) sets."""
import numpy as np, pandas as pd, glob, os
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
SC="/path/to/scratch"
rng=np.random.default_rng(4242)
CPM=pd.read_csv(f"{OUT}/sc_hbca_global_pseudobulk_by_celltype.csv",index_col=0)
cnt=pd.read_csv(f"{OUT}/sc_hbca_global_celltype_counts.csv")
keep=cnt[cnt.n_cells>=500].cell_type.tolist()
C=CPM[keep]; C=C[C.sum(axis=1)>0]
share=C.div(C.sum(axis=1),axis=0)
FIB=[c for c in keep if c.startswith("CAFs") or c=="Myofibroblast"]
fib=share[FIB].sum(axis=1); mal=share["malignant cell"]
print("genes",C.shape[0],"| fibroblast columns:",FIB)
sets={}
for f in glob.glob(f"{SC}/set_*.csv"):
    sets[os.path.basename(f)[4:-4]]=pd.read_csv(f).gene.tolist()
a29=pd.read_csv(f"{OUT}/screens_mir29_target_set.csv"); a130=pd.read_csv(f"{OUT}/mir130a_target_set.csv")
sets["mir29_strong"]=a29.query("tier=='STRONG_lowthroughput'").gene.tolist()
sets["mir130a_strong"]=a130.query("tier=='STRONG_lowthroughput'").gene.tolist()
sets["mir29_targetscan"]=a29.query("in_targetscan==True").gene.tolist()
sets["mir130a_targetscan"]=a130.query("in_targetscan==True").gene.tolist()
tot=C.sum(axis=1); dec=pd.qcut(tot.rank(method="first"),10,labels=False)
B=2000; rows=[]
store={}
for name,gl in sets.items():
    g=[x for x in dict.fromkeys(gl) if x in C.index]
    if len(g)<10: continue
    pool={d:C.index[(dec==d)&(~C.index.isin(set(gl)))].values for d in range(10)}
    gd=dec.loc[g].values
    of=float(fib.loc[g].mean()); om=float(mal.loc[g].mean())
    nf=np.empty(B); nm=np.empty(B)
    for b in range(B):
        pick=[rng.choice(pool[d]) for d in gd]
        nf[b]=fib.loc[pick].mean(); nm[b]=mal.loc[pick].mean()
    store[name]=(g,gd,of,om)
    rows.append(dict(set=name,n=len(g),fib_obs=of,fib_null=nf.mean(),fib_z=(of-nf.mean())/nf.std(ddof=1),
                     fib_p=(1+np.sum(nf>=of))/(1+B),
                     mal_obs=om,mal_null=nm.mean(),mal_z=(om-nm.mean())/nm.std(ddof=1),
                     mal_p=(1+np.sum(nm>=om))/(1+B)))
R=pd.DataFrame(rows).sort_values("set")
print("\n=== compartment shares, RAW vs CAF-ADJUSTED sets (independent) ===")
print(R.to_string(index=False,float_format=lambda x:"%.4f"%x))
R.to_csv(f"{OUT}/sc_verify_compartment_decircularised.csv",index=False)
# formal two-set contrast: permutation on the pooled gene labels, decile-stratified
def contrast(a,b,B=5000):
    ga=[x for x in dict.fromkeys(sets[a]) if x in C.index]; gb=[x for x in dict.fromkeys(sets[b]) if x in C.index]
    obs=fib.loc[ga].mean()-fib.loc[gb].mean(); obsm=mal.loc[ga].mean()-mal.loc[gb].mean()
    pooled=np.array(ga+gb); na=len(ga); d=np.empty(B); dm=np.empty(B)
    for i in range(B):
        p=rng.permutation(pooled); d[i]=fib.loc[p[:na]].mean()-fib.loc[p[na:]].mean(); dm[i]=mal.loc[p[:na]].mean()-mal.loc[p[na:]].mean()
    return dict(pair=f"{a} - {b}",n_a=na,n_b=len(gb),d_fib=obs,p_fib=(1+np.sum(np.abs(d)>=abs(obs)))/(1+B),
                d_mal=obsm,p_mal=(1+np.sum(np.abs(dm)>=abs(obsm)))/(1+B))
cons=[contrast("mir29_strong","mir130a_strong"),contrast("mir29_targetscan","mir130a_targetscan"),
      contrast("mir29_anticorr_RAW","mir130a_anticorr_RAW"),contrast("mir29_anticorr_ADJ","mir130a_anticorr_ADJ")]
Cn=pd.DataFrame(cons)
print("\n=== FORMAL TWO-SET CONTRAST (miR-29 set minus miR-130a set), label permutation ===")
print(Cn.to_string(index=False,float_format=lambda x:"%.4f"%x))
Cn.to_csv(f"{OUT}/sc_verify_compartment_setcontrast.csv",index=False)
