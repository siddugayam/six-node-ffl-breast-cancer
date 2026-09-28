#!/usr/bin/env python3
"""Programme assignment and cross-cutting statistics for the compendium synthesis."""
import pandas as pd, numpy as np, os, itertools
from scipy import stats
ROOT="/path/to/revision"; V5=ROOT+"/results/v5"
T=pd.read_csv(V5+"/node_compendium_table.csv")
PROG={
"P1 Cell-cycle progression and genome integrity":["E2F1","E2F3","MYBL2","CCND2","BRCA1"],
"P2 Chromatin state and DNA methylation":["EZH2","DNMT1","hsa-miR-101","hsa-miR-195","hsa-miR-124"],
"P3 Luminal identity and hormone signalling":["ESR1","GATA3","STAT5A"],
"P4 Matrix remodelling and invasion":["COL1A1","FN1","MMP14","PLAU","MET"],
"P5 Stromal and inflammatory microenvironment":["PDGFRB","CXCL12","JUN","EGR2","hsa-miR-155","hsa-miR-21"],
"P6 Epithelial plasticity and TGF-beta signalling":["hsa-miR-429","hsa-miR-141","hsa-miR-204"],
"P7 Metabolic and stress-response control":["SREBF1","hsa-miR-383","hsa-miR-34a"],
}
assert sum(len(v) for v in PROG.values())==30
assert sorted(itertools.chain(*PROG.values()))==sorted(T.name)
T["programme"]=""
for k,v in PROG.items(): T.loc[T.name.isin(v),"programme"]=k

# ---- rule-based cis-mechanism call (identical to the rule used in 44_write_compendium.py) ----
def mech_call(r):
    if pd.isna(r.fdr_TCGA) or r.fdr_TCGA >= 0.05:
        return "not_DE"
    c = []; lf = r.logFC_TCGA
    if (pd.notna(r.meth_delta_beta) and abs(r.meth_delta_beta) >= 0.05 and r.meth_fdr < 0.05
            and np.sign(r.meth_delta_beta) == -np.sign(lf)
            and pd.notna(r.meth_expr_rho) and r.meth_expr_rho < -0.1):
        c.append("meth_supported")
    elif (pd.notna(r.meth_delta_beta) and abs(r.meth_delta_beta) >= 0.05 and r.meth_fdr < 0.05
          and np.sign(r.meth_delta_beta) == np.sign(lf)):
        c.append("meth_discordant")
    if (pd.notna(r.cn_frac_amp) and r.cn_frac_amp >= .30 and pd.notna(r.cn_expr_rho)
            and r.cn_expr_rho >= .20 and r.cn_expr_fdr < .05 and lf > 0):
        c.append("gain")
    if (pd.notna(r.cn_frac_del) and r.cn_frac_del >= .30 and pd.notna(r.cn_expr_rho)
            and r.cn_expr_rho >= .20 and r.cn_expr_fdr < .05 and lf < 0):
        c.append("loss")
    if pd.notna(r.mut_freq) and r.mut_freq >= 0.05:
        c.append("mutation")
    return ";".join(c) if c else "none"

T["mech"] = T.apply(mech_call, axis=1)
T[["name", "type", "programme", "mech"]].to_csv(V5 + "/node_mechanism_calls.csv", index=False)
T.to_csv(V5 + "/node_compendium_table.csv", index=False)
print("programme and mechanism columns added; node_compendium_table.csv rewritten\n")
print("\nmechanism calls:")
print(T.groupby("mech").name.apply(list).to_string())


g=T.groupby("programme")
S=pd.DataFrame(dict(
 n=g.size(), n_TF=g.apply(lambda d:(d.type=="TF").sum()), n_gene=g.apply(lambda d:(d.type=="Gene").sum()),
 n_miRNA=g.apply(lambda d:(d.type=="miRNA").sum()),
 median_logFC=g.logFC_TCGA.median().round(2), n_up=g.apply(lambda d:(d.logFC_TCGA>0).sum()),
 median_ffl=g.ffl_total.median(), median_degree=g.degree_total.median(),
 median_deltabeta=g.meth_delta_beta.median().round(3),
 mean_chronos=g.chronos_mean_breast.mean().round(3),
 n_OSp05=g.apply(lambda d:(d.tcga_OS_p<0.05).sum()),
 n_PFIp05=g.apply(lambda d:(d.tcga_PFI_p<0.05).sum()),
 n_MBq05=g.apply(lambda d:(d.metabric_OS_q<0.05).sum()),
 median_frac_strong=g.frac_strong.median().round(3)))
S.to_csv(V5+"/node_compendium_programme_summary.csv")
print(S.to_string())

# within/between programme edges among the 30
P=pd.read_csv(V5+"/node_compendium_partners_full.csv")
prog=dict(zip(T.name,T.programme))
pp=P[(P.role=="regulates")&(P.partner.isin(set(T.name)))].copy()
pp["prog_src"]=pp.node.map(prog); pp["prog_tgt"]=pp.partner.map(prog)
pp["within"]=pp.prog_src==pp.prog_tgt
print("\nprioritised->prioritised edges:",len(pp),"within-programme:",int(pp.within.sum()))
print(pp[["node","partner","edge_type","evidence_tier","rho","rho_caf_adj","prog_src","prog_tgt"]]
      .sort_values("rho").round(3).to_string(index=False))
