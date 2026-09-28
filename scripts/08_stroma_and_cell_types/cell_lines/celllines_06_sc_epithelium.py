#!/usr/bin/env python3
"""
celllines_06_sc_epithelium.py

Control for part C: are breast cancer CELL LINES collagen-negative because they are
cultured, or because tumour epithelial cells are collagen-negative in vivo too?

Pseudobulk CPM of COL1A1 / COL3A1 (+ controls) per cell type in the Wu et al. 2021
breast cancer single-cell atlas (GSE176078, 100,064 cells, 29,733 genes), computed here
from the raw sparse count matrix.  CPM = sum(gene counts in cell type) /
sum(nCount_RNA in cell type) * 1e6.
"""
import os
import numpy as np
import pandas as pd

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
SC = os.path.join(ROOT, "cache", "celltype", "Wu_etal_2021_BRCA_scRNASeq")

GENE_IDX = {17587: "COL1A1", 3199: "COL3A1", 8041: "COL1A2", 2443: "EPCAM",
            1727: "PTPRC", 12041: "VIM", 12401: "ACTA2", 3347: "FN1", 6234: "SPARC",
            5937: "LOX", 11906: "ETS1", 5127: "NFKB1", 11385: "RELA", 13233: "SP1",
            9869: "MYC", 13216: "KRT18"}

genes = [l.rstrip("\n") for l in open(os.path.join(SC, "count_matrix_genes.tsv"))]
bcs = [l.rstrip("\n") for l in open(os.path.join(SC, "count_matrix_barcodes.tsv"))]
print("genes", len(genes), "barcodes", len(bcs))
for i, g in GENE_IDX.items():
    assert genes[i - 1] == g, (i, genes[i - 1], g)
print("gene index assertions passed")

meta = pd.read_csv(os.path.join(SC, "metadata.csv"), index_col=0)
print("metadata", meta.shape, list(meta.columns))
meta = meta.loc[bcs]                     # enforce mtx column order
assert list(meta.index) == bcs
ct = meta["celltype_major"].values
ctm = meta["celltype_minor"].values
tot = meta["nCount_RNA"].values.astype(float)

ent = pd.read_csv(os.path.join(RES, "_sc_selected_entries.tsv"), sep="\t",
                  header=None, names=["gi", "ci", "count"])
print("selected non-zero entries:", len(ent))
ent["gene"] = ent["gi"].map(GENE_IDX)
ent["celltype_major"] = ct[ent["ci"].values - 1]
ent["celltype_minor"] = ctm[ent["ci"].values - 1]

for key in ["celltype_major", "celltype_minor"]:
    tot_by = pd.Series(tot).groupby(pd.Series(meta[key].values)).sum()
    n_by = pd.Series(meta[key].values).value_counts()
    s = ent.groupby([key, "gene"])["count"].sum().unstack(fill_value=0)
    ncells_pos = ent.groupby([key, "gene"])["ci"].nunique().unstack(fill_value=0)
    cpm = s.div(tot_by.reindex(s.index), axis=0) * 1e6
    pct = ncells_pos.div(n_by.reindex(ncells_pos.index), axis=0) * 100
    out = cpm.round(3).add_suffix("_CPM").join(pct.round(2).add_suffix("_pct_cells_pos"))
    out.insert(0, "n_cells", n_by.reindex(out.index).values)
    out.insert(1, "total_UMI", tot_by.reindex(out.index).values)
    out = out.sort_values("COL1A1_CPM", ascending=False)
    out.to_csv(os.path.join(RES, f"celllines_C_scRNA_{key}.csv"))
    print(f"\n===== {key} =====")
    cols = ["n_cells", "COL1A1_CPM", "COL3A1_CPM", "COL1A2_CPM", "EPCAM_CPM",
            "COL1A1_pct_cells_pos", "COL3A1_pct_cells_pos"]
    print(out[cols].to_string())

# headline ratios
maj = pd.read_csv(os.path.join(RES, "celllines_C_scRNA_celltype_major.csv"), index_col=0)
if "CAFs" in maj.index and "Cancer Epithelial" in maj.index:
    for g in ["COL1A1", "COL3A1"]:
        a = maj.loc["CAFs", g + "_CPM"]; b = maj.loc["Cancer Epithelial", g + "_CPM"]
        print(f"{g}: CAF {a:.1f} CPM vs Cancer Epithelial {b:.1f} CPM -> {a/b:.1f}x")
print("DONE")
