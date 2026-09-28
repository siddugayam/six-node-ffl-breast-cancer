#!/usr/bin/env python3
"""(b)/(c) Do ETS1/NFKB1/RELA/SP1 motifs occur in accessible chromatin at
COL1A1/COL3A1, and at targets whose correlation survives compartment adjustment?

Null: GC-decile x accessibility-decile matched TCGA-BRCA ATAC peaks drawn from
outside a 1 Mb exclusion zone around the tested genes; 10,000 resamples.
Correction: Benjamini-Hochberg across all (gene set x motif) tests.
"""
import os, numpy as np, pandas as pd
ROOT="/path/to/revision"; C=f"{ROOT}/cache/v6/atac"; R=f"{ROOT}/results/v6"
rng=np.random.default_rng(20260910)
B=10000

hits=np.load(f"{C}/motif_hits.npy")
mmeta=pd.read_csv(f"{C}/motif_meta.tsv",sep="\t")
groups=pd.read_csv(f"{C}/motif_panel_groups.tsv",sep="\t").set_index("matrix_id")
meta=pd.read_csv(f"{C}/peak_stroma_correlation.tsv.gz",sep="\t")
tss=pd.read_csv(f"{C}/refseq_select_tss.tsv",sep="\t").drop_duplicates("symbol").set_index("symbol")
pc=pd.read_csv(f"{R}/atac_tf_target_compartment_adjusted.csv")
assert len(meta)==hits.shape[0]
mmeta["group"]=[groups.loc[m,"group"] if m in groups.index else "NA" for m in mmeta.matrix_id]
FOCAL={"ETS1":"MA0098.4","NFKB1":"MA0105.4","RELA":"MA0107.1","SP1":"MA0079.5"}

# strata: GC decile x accessibility decile
gcq=pd.qcut(meta.percentGC,10,labels=False,duplicates="drop")
acq=pd.qcut(meta.mean_acc,10,labels=False,duplicates="drop")
stratum=(gcq.astype(int)*10+acq.astype(int)).values
meta["stratum"]=stratum
chrom=meta.seqnames.values; mid=((meta.start+meta.end)//2).values

def peaks_for(genes,window):
    idx=[]
    for g in genes:
        if g not in tss.index: continue
        t=tss.loc[g]
        sel=np.where((chrom==t.chrom)&(np.abs(mid-t.tss)<=window))[0]
        idx.extend(sel.tolist())
    return np.array(sorted(set(idx)),dtype=int)

def exclusion(genes,pad=1_000_000):
    m=np.zeros(len(meta),dtype=bool)
    for g in genes:
        if g not in tss.index: continue
        t=tss.loc[g]
        m |= (chrom==t.chrom)&(np.abs(mid-t.tss)<=pad)
    return m

def test(focal_idx,excl_mask,mi,B=B):
    """binary: number of focal peaks with >=1 hit, vs matched background."""
    n=len(focal_idx)
    obs=int((hits[focal_idx,mi]>0).sum())
    bgmask=~excl_mask
    draws=np.zeros(B,dtype=np.int32)
    # per-stratum pools
    pools={}
    for s in np.unique(stratum[focal_idx]):
        pool=np.where((stratum==s)&bgmask)[0]
        if len(pool)==0: pool=np.where(bgmask)[0]
        pools[s]=(hits[pool,mi]>0).astype(np.int8)
    for s in np.unique(stratum[focal_idx]):
        k=int((stratum[focal_idx]==s).sum())
        p=pools[s]
        draws+=p[rng.integers(0,len(p),size=(B,k))].sum(1)
    exp=draws.mean()
    p_hi=(draws>=obs).mean(); p_lo=(draws<=obs).mean()
    p=min(1.0,2*min(p_hi,p_lo))
    return dict(n_peaks=n,obs_peaks_with_hit=obs,frac_obs=obs/n if n else np.nan,
                exp_peaks_with_hit=exp,frac_exp=exp/n if n else np.nan,
                ratio=(obs/exp) if exp>0 else np.nan,p_emp=p,p_greater=p_hi)

SETS=[]
SETS.append(("COL1A1_promoter_5kb",["COL1A1"],5000))
SETS.append(("COL1A1_locus_100kb",["COL1A1"],100000))
SETS.append(("COL3A1_promoter_5kb",["COL3A1"],5000))
SETS.append(("COL3A1_locus_100kb",["COL3A1"],100000))
SETS.append(("collagen_pair_locus_100kb",["COL1A1","COL3A1"],100000))
for tf in FOCAL:
    tg=sorted(pc.loc[(pc.TF==tf)&(pc.poscontrol),"target"].tolist())
    if tg:
        SETS.append((f"poscontrol_{tf}_promoter_5kb",tg,5000))
        SETS.append((f"poscontrol_{tf}_locus_100kb",tg,100000))
allpos=sorted(set(pc.loc[pc.poscontrol,"target"]))
SETS.append(("poscontrol_ALL_promoter_5kb",allpos,5000))
SETS.append(("poscontrol_ALL_locus_100kb",allpos,100000))

rows=[]
for name,genes,win in SETS:
    fi=peaks_for(genes,win); ex=exclusion(genes)
    if len(fi)==0:
        rows.append(dict(set=name,genes=";".join(genes)[:120],window=win,motif="ALL",motif_name="ALL",
                         group="NA",n_peaks=0)); continue
    for mi,mrow in mmeta.iterrows():
        r=test(fi,ex,mi)
        r.update(set=name,genes=";".join(genes)[:200],n_genes=len(genes),window=win,
                 motif=mrow.matrix_id,motif_name=mrow["name"],group=mrow.group)
        rows.append(r)
    print(f"{name}: {len(fi)} peaks over {len(genes)} genes",flush=True)
D=pd.DataFrame(rows)
p=D.p_emp.fillna(1).values; o=np.argsort(p); q=np.empty(len(p)); m=len(p)
q[o]=np.minimum.accumulate((p[o]*m/(np.arange(m)+1))[::-1])[::-1]
D["fdr_BH"]=np.clip(q,0,1)
D.to_csv(f"{R}/atac_motif_enrichment.csv",index=False)

pd.set_option("display.width",250)
foc=[v for v in FOCAL.values()]
print("\n=== focal TF motifs at the collagen loci ===")
print(D[(D.set.str.startswith("COL"))&(D.motif.isin(foc))][
    ["set","motif_name","n_peaks","obs_peaks_with_hit","frac_obs","frac_exp","ratio","p_emp","fdr_BH"]].to_string(index=False))
print("\n=== focal TF motifs at their own positive-control targets ===")
sub=D[D.set.str.startswith("poscontrol")]
sub=sub[[ (s.split("_")[1]==m) or s.startswith("poscontrol_ALL") for s,m in zip(sub.set,sub.motif_name)]]
print(sub[["set","motif_name","n_genes","n_peaks","obs_peaks_with_hit","frac_obs","frac_exp","ratio","p_emp","fdr_BH"]].to_string(index=False))
print("\n=== strongest enrichments at the collagen loci (any motif, FDR<0.1) ===")
cc=D[(D.set.str.startswith("COL")|D.set.str.startswith("collagen"))&(D.fdr_BH<0.1)].sort_values("p_emp")
print(cc[["set","motif_name","group","n_peaks","frac_obs","frac_exp","ratio","p_emp","fdr_BH"]].head(25).to_string(index=False))
print("DONE")
