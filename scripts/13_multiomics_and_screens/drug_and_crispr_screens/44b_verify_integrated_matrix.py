#!/usr/bin/env python3
"""Independent verification of PART 1E (integrated essentiality matrix)."""
import pandas as pd, numpy as np
from scipy.stats import spearmanr
REV="/path/to/revision"; OUT=f"{REV}/results/v3"
m=pd.read_csv(f"{OUT}/screens_integrated_essentiality_matrix.csv")
print("matrix:",m.shape)
print("\nconsensus counts:"); print(m.consensus.value_counts().to_string())
print("\nn_modalities:"); print(m.n_modalities.value_counts().to_string())
print("\nessential in EVERY modality tested:",m.loc[(m.n_modalities>0)&(m.n_essential==m.n_modalities),"node"].tolist())
print("essential in SOME:",int(((m.n_essential>0)&(m.n_essential<m.n_modalities)).sum()))
print("non-essential in every modality:",int(((m.n_modalities>0)&(m.n_essential==0)).sum()))
print("not screenable:",int((m.n_modalities==0).sum()),"| of which miRNA:",int(((m.n_modalities==0)&(m.type=='miRNA')).sum()))
h=m[m.is_FFL_hub==True]
print("\nFFL hubs:",len(h),"| screened:",int((h.n_modalities>0).sum()),"| essential in any:",h.loc[h.n_essential>0,"node"].tolist())
s=m.dropna(subset=["chronos_breast","d2_breast","score_breast"])
print("\nnodes with all three modalities:",len(s))
print("CRISPR vs RNAi  Spearman %+.4f"%spearmanr(s.chronos_breast,s.d2_breast).statistic)
print("CRISPR vs SCORE Spearman %+.4f"%spearmanr(s.chronos_breast,s.score_breast).statistic)
print("RNAi   vs SCORE Spearman %+.4f"%spearmanr(s.d2_breast,s.score_breast).statistic)
foc=["NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","MYC","HIF1A","STAT3","TP53","ESR1","MKL1","VEGFA"]
out=m[m.node.isin(foc)][["node","type","is_FFL_hub","chronos_breast","d2_breast","score_breast","consensus","orcs_breast_screens_tested","orcs_breast_screens_hit"]]
print("\nfocus nodes:"); print(out.to_string(index=False,float_format=lambda x:"%.4f"%x))
out.to_csv(f"{OUT}/screens_verify_integrated_focus.csv",index=False)
