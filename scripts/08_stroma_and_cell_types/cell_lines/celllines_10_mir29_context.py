#!/usr/bin/env python3
"""
celllines_10_mir29_context.py

Where does miR-29 sit relative to the stromal / mesenchymal axis in each system?
 - TCGA-BRCA bulk: miR-29a/b/c vs the 4-marker CAF score and vs COL1A1/COL3A1
 - CCLE breast lines: miR-29a/b/c vs the mesenchymal score and vs COL1A1/COL3A1
This is what decides whether the bulk miR-29 -> collagen correlation can be read as
cell-intrinsic repression.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")

zc = lambda df: (df - df.mean()) / df.std(ddof=1)
rows = []

# TCGA
G = pd.read_csv(os.path.join(RES, "_tcga_targets.tsv"), sep="\t", index_col=0).T
M = pd.read_csv(os.path.join(RES, "_tcga_allmirs.tsv"), sep="\t", index_col=0).T
G = G.loc[M.index]
caf = zc(G[["DCN", "LUM", "FAP", "THY1"]]).mean(axis=1)
for m in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-miR-200c", "hsa-miR-21"]:
    if m not in M.columns:
        continue
    r, p = stats.spearmanr(M[m].values, caf.values)
    rows.append(dict(system="TCGA-BRCA paired tumours", n=len(caf), miRNA=m,
                     partner="CAF score (DCN/LUM/FAP/THY1)", rho=float(r), p=float(p)))
    for c in ["COL1A1", "COL3A1"]:
        r, p = stats.spearmanr(M[m].values, G[c].values)
        rows.append(dict(system="TCGA-BRCA paired tumours", n=len(caf), miRNA=m,
                         partner=c, rho=float(r), p=float(p)))

# cell lines
E = pd.read_csv(os.path.join(RES, "celllines_expr_breast.tsv.gz"), sep="\t", index_col=0)
mir = np.log2(pd.read_csv(os.path.join(RES, "celllines_mirna_breast.tsv"),
                          sep="\t", index_col=0) + 1)
ids = [i for i in mir.index if i in E.index]
E, mir = E.loc[ids], mir.loc[ids]
MES = [g for g in ["VIM", "ZEB1", "ZEB2", "SNAI2", "TWIST1", "CDH2", "SPARC", "THY1"]
       if g in E.columns]
mes = zc(E[MES]).mean(axis=1)
for m in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-miR-200c", "hsa-miR-21"]:
    if m not in mir.columns:
        continue
    r, p = stats.spearmanr(mir[m].values, mes.values)
    rows.append(dict(system="CCLE breast cancer lines", n=len(ids), miRNA=m,
                     partner="mesenchymal score (%s)" % "/".join(MES),
                     rho=float(r), p=float(p)))
    for c in ["COL1A1", "COL3A1"]:
        r, p = stats.spearmanr(mir[m].values, E[c].values)
        rows.append(dict(system="CCLE breast cancer lines", n=len(ids), miRNA=m,
                         partner=c, rho=float(r), p=float(p)))

d = pd.DataFrame(rows)
d.to_csv(os.path.join(RES, "celllines_F_mir29_context.csv"), index=False)
print(d.round(4).to_string(index=False))
print("DONE")
