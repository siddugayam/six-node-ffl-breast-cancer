#!/usr/bin/env python3
"""
(a) Matched PDX-vs-originator log2 fold change, corrected for the compositional
    renormalisation that loss of the stromal transcriptome causes (each pair's
    log2 ratios are centred on the median log2 ratio of the reference panel).
(b) CAF / cancer-epithelial expression ratio for the same genes from the
    Wu et al. 2021 breast atlas pseudobulk (cache/celltype/scrna_pseudobulk_all.tsv).
(c) Does single-cell stromality predict the PDX collapse?
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE="/path/to/revision"
CACHE=f"{BASE}/cache/v6/pdx"; RES=f"{BASE}/results/v6"

sym2ent=json.load(open(f"{CACHE}/gene_symbol_entrez.json"))
ent2sym={}
for s,e in sym2ent.items(): ent2sym.setdefault(e,s)
meta=json.load(open("/tmp/pdmr_breast_meta.json"))
long=pd.read_csv(f"{CACHE}/pdmr_breast_expr_long.tsv", sep="\t")
W=long.pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="mean")
W.columns=[ent2sym.get(c,str(c)) for c in W.columns]
W=W.loc[:,~pd.Index(W.columns).duplicated()]

ref=json.load(open("/tmp/ref_panel_entrez.json"))
TARGET=["COL1A1","COL3A1","FN1","PDGFRB","CXCL12","POSTN","MMP14","PLAU","MET","STAT5A",
        "CCND2","MYBL2","E2F1","E2F3","EZH2","DNMT1","BRCA1","GATA3","ESR1","JUN","EGR2","SREBF1",
        "NFKB1","ETS1","EPCAM","KRT8","KRT18","KRT19","ERBB2","PTPRC","PECAM1","ACTA2","DCN","LUM",
        "VIM","MKI67","COL1A2","THY1","FAP","PDGFRA"]
TARGET=[g for g in TARGET if g in W.columns]
REF=[g for g in sorted(set(ref)&set(W.columns)) if g not in TARGET]

pas=pd.Series({s:meta[s].get("PASSAGE") for s in W.index})
mod=pd.Series({s:meta[s].get("MODEL_ID") for s in W.index})
is_org=pas=="Originator"
is_pdx=pas.apply(lambda p: isinstance(p,str) and p.startswith("P") and p[1:].isdigit())

models=sorted(set(mod[is_org]) & set(mod[is_pdx]))
rows=[]
for m in models:
    o=W.loc[is_org & (mod==m)].mean(axis=0)
    x=W.loc[is_pdx & (mod==m)].mean(axis=0)
    lr=np.log2((x+1)/(o+1))
    off=np.nanmedian(lr[REF].values)          # compositional offset
    rows.append(pd.Series({**{g:lr[g] for g in TARGET}, "_offset":off}, name=m))
FC=pd.DataFrame(rows)
FC.to_csv(f"{RES}/pdx_matched_log2fc_permodel.csv")

res=[]
for g in TARGET:
    raw=FC[g].values; adj=FC[g].values-FC["_offset"].values
    w,pw=stats.wilcoxon(adj)
    res.append(dict(gene=g, n_models=len(FC),
        median_log2FC_raw=float(np.median(raw)),
        median_log2FC_panelnormalised=float(np.median(adj)),
        mean_log2FC_panelnormalised=float(np.mean(adj)),
        n_down=int((adj<0).sum()), n_up=int((adj>0).sum()), wilcoxon_p=pw))
R=pd.DataFrame(res)
R["wilcoxon_q"]=multipletests(R.wilcoxon_p.values, method="fdr_bh")[1]
R["median_compositional_offset_log2"]=float(np.median(FC["_offset"]))
R=R.sort_values("median_log2FC_panelnormalised")
R.to_csv(f"{RES}/pdx_matched_log2fc_summary.csv", index=False)
print(R.round(3).to_string(index=False))

# ---- (b) Wu 2021 atlas CAF vs cancer-epithelial ---------------------------
pb=pd.read_csv(f"{BASE}/cache/celltype/scrna_pseudobulk_all.tsv", sep="\t")
codes=pd.read_csv(f"{BASE}/cache/celltype/celltype_codes.tsv", sep="\t", header=None,
                  names=["ct_code","cell_type"])
pb=pb.merge(codes, on="ct_code")
tot=pb.groupby("cell_type").sum_umi.sum()
pb["cpm"]=pb.sum_umi/pb.cell_type.map(tot)*1e6
piv=pb.pivot_table(index="gene", columns="cell_type", values="cpm", aggfunc="sum").fillna(0)
piv["log2_CAF_over_CancerEpi"]=np.log2((piv["CAFs"]+1)/(piv["Cancer Epithelial"]+1))
piv[["CAFs","Cancer Epithelial","log2_CAF_over_CancerEpi"]].to_csv(f"{RES}/pdx_wu_atlas_caf_ratio_allgenes.csv")

# ---- (c) does stromality predict the PDX collapse? ------------------------
rk=pd.read_csv(f"{RES}/pdx_matched_pairs_rank_change.csv")
J=rk.merge(piv[["CAFs","Cancer Epithelial","log2_CAF_over_CancerEpi"]],
           left_on="gene", right_index=True, how="inner")
J=J.merge(R[["gene","median_log2FC_panelnormalised"]], on="gene", how="left")
J.to_csv(f"{RES}/pdx_stromality_vs_rankdrop.csv", index=False)
for yc in ["median_paired_delta","median_log2FC_panelnormalised"]:
    d=J.dropna(subset=["log2_CAF_over_CancerEpi",yc])
    r,p=stats.spearmanr(d.log2_CAF_over_CancerEpi, d[yc])
    rp,pp=stats.pearsonr(d.log2_CAF_over_CancerEpi, d[yc])
    print(f"{yc}: n={len(d)} spearman rho={r:.3f} p={p:.3g} | pearson r={rp:.3f} p={pp:.3g}")

# whole-panel version (all genes fetched, not just the 40)
rk_all=pd.read_csv(f"{RES}/pdx_percentile_ranks_persample.csv", index_col=0)
print("\ntop stromal-ratio genes among targets:")
print(J.sort_values("log2_CAF_over_CancerEpi", ascending=False)[
      ["gene","CAFs","Cancer Epithelial","log2_CAF_over_CancerEpi",
       "median_paired_delta","median_log2FC_panelnormalised"]].round(3).to_string(index=False))
print("DONE")
