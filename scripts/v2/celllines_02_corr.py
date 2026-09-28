#!/usr/bin/env python3
"""
celllines_02_corr.py -- the stroma-free test, parts B and C.

B) In breast cancer cell lines (no fibroblasts present):
     COL1A1 ~ COL3A1
     ETS1 / NFKB1 / RELA / SP1 / MYB / MKL1 / STAT6 / TFAP2A -> COL1A1, COL3A1
     miR-29a/b/c, let-7b/e, miR-143, miR-218, miR-133a/b -> COL1A1, COL3A1
   with the matched TCGA-BRCA bulk values recomputed here in the same run.

C) Absolute expression of COL1A1 / COL3A1 in breast lines vs other lineages and
   vs bulk tumours (rank-based, because the two platforms are not on a common scale).

Outputs: results/v2/celllines_B_*.csv, celllines_C_*.csv
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")

TFS = ["ETS1", "NFKB1", "RELA", "SP1", "MYB", "MKL1", "STAT6", "TFAP2A"]
MIRS = ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-let-7b", "hsa-let-7e",
        "hsa-miR-143", "hsa-miR-218", "hsa-miR-133a", "hsa-miR-133b", "hsa-miR-101",
        "hsa-miR-21", "hsa-miR-130a"]
COLS = ["COL1A1", "COL3A1"]


def spear(x, y, min_n=8):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < min_n:
        return dict(n=n, rho=np.nan, p=np.nan, r_pearson=np.nan, p_pearson=np.nan)
    rho, p = stats.spearmanr(x[m], y[m])
    r, pp = stats.pearsonr(x[m], y[m])
    return dict(n=n, rho=float(rho), p=float(p), r_pearson=float(r), p_pearson=float(pp))


# ------------------------------------------------------------------ cell lines
expr = pd.read_csv(os.path.join(RES, "celllines_expr_breast.tsv.gz"),
                   sep="\t", index_col=0)
models = pd.read_csv(os.path.join(RES, "celllines_breast_models.csv"))
models = models.set_index("ModelID")
cancer_ids = [i for i in expr.index
              if models.loc[i, "OncotreePrimaryDisease"] != "Non-Cancerous"]
print(f"breast models with expression: {expr.shape[0]}  (cancer only: {len(cancer_ids)})")

mir = pd.read_csv(os.path.join(RES, "celllines_mirna_breast.tsv"), sep="\t", index_col=0)
print(f"breast models with miRNA: {mir.shape[0]}")
print(f"CCLE miRNA raw value range {np.nanmin(mir.values):.3f} .. {np.nanmax(mir.values):.1f}"
      " -> log2 transformed (Spearman is unaffected; Pearson uses the log scale)")
mir = np.log2(mir + 1.0)

rows = []


def add(dataset, pair_type, a, b, n_note, d, extra=None):
    r = dict(dataset=dataset, pair_type=pair_type, regulator=a, target=b, subset=n_note)
    r.update(d)
    if extra:
        r.update(extra)
    rows.append(r)


for subset_name, ids in [("all_breast_lines", list(expr.index)),
                         ("breast_cancer_lines_only", cancer_ids)]:
    E = expr.loc[ids]
    add("CCLE/DepMap_24Q4", "gene_gene", "COL1A1", "COL3A1", subset_name,
        spear(E["COL1A1"].values, E["COL3A1"].values))
    for tf in TFS:
        if tf not in E.columns:
            print("MISSING TF", tf)
            continue
        for c in COLS:
            add("CCLE/DepMap_24Q4", "TF_target", tf, c, subset_name,
                spear(E[tf].values, E[c].values))

# miRNA, matched on ModelID
both = [i for i in mir.index if i in expr.index]
both_ca = [i for i in both if i in set(cancer_ids)]
print(f"breast models with miRNA AND expression: {len(both)}  (cancer only {len(both_ca)})")
for subset_name, ids in [("miRNA_matched_all", both), ("miRNA_matched_cancer_only", both_ca)]:
    M = mir.loc[ids]
    E = expr.loc[ids]
    for m in MIRS:
        if m not in M.columns:
            print("MISSING miRNA probe", m)
            continue
        for c in COLS:
            add("CCLE_nanoString_miRNA + DepMap_24Q4", "miRNA_target", m, c, subset_name,
                spear(M[m].values, E[c].values))
    # positive control: a miRNA-target pair that should work anywhere
    for m, c in [("hsa-miR-101", "EZH2")]:
        if m in M.columns and c in E.columns:
            add("CCLE_nanoString_miRNA + DepMap_24Q4", "miRNA_target(control)", m, c,
                subset_name, spear(M[m].values, E[c].values))

# ------------------------------------------------------------------ TCGA bulk (recomputed here)
import subprocess
rcode = r'''
g <- readRDS("%s/brca_gene_expr.rds")
m <- readRDS("%s/brca_mirna_expr_canonical.rds")
p <- readRDS("%s/brca_pheno.rds")
tum <- p$sample[p$sample_type == "Primary Tumor"]
gt <- intersect(colnames(g), tum); mt <- intersect(colnames(m), tum)
both <- intersect(gt, mt)
cat("N_gene_tumours", length(gt), "\n"); cat("N_paired", length(both), "\n")
genes <- c("COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1","MYB","MKL1","STAT6","TFAP2A",
           "EZH2","DCN","LUM","FAP","THY1")
genes <- genes[genes %%in%% rownames(g)]
write.table(g[genes, gt], "%s/_tcga_genes_tumour.tsv", sep="\t", quote=FALSE)
mirs <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-let-7b","hsa-let-7e",
          "hsa-miR-143","hsa-miR-218","hsa-miR-133a","hsa-miR-133b","hsa-miR-101",
          "hsa-miR-21","hsa-miR-130a")
mirs <- mirs[mirs %%in%% rownames(m)]
write.table(m[mirs, both], "%s/_tcga_mirs_tumour.tsv", sep="\t", quote=FALSE)
write.table(g[genes, both], "%s/_tcga_genes_paired.tsv", sep="\t", quote=FALSE)
# within-sample percentile rank of the collagens across ALL genes
r <- apply(g[, gt], 2, function(v) rank(v)/length(v))
write.table(t(r[c("COL1A1","COL3A1"), ]), "%s/_tcga_collagen_pct.tsv", sep="\t", quote=FALSE)
cat("N_genes_matrix", nrow(g), "\n")
''' % ((DATA,) * 3 + (RES,) * 4)
tmp_r = os.path.join(ROOT, "scripts", "v2", "_celllines_tcga_pull.R")
open(tmp_r, "w").write(rcode)
out = subprocess.run(["/usr/bin/Rscript", tmp_r], capture_output=True, text=True)
print(out.stdout, out.stderr)
assert out.returncode == 0, out.stderr

tg = pd.read_csv(os.path.join(RES, "_tcga_genes_tumour.tsv"), sep="\t", index_col=0).T
tgp = pd.read_csv(os.path.join(RES, "_tcga_genes_paired.tsv"), sep="\t", index_col=0).T
tm = pd.read_csv(os.path.join(RES, "_tcga_mirs_tumour.tsv"), sep="\t", index_col=0).T
print("TCGA gene tumours", tg.shape, "paired", tgp.shape, tm.shape)

add("TCGA-BRCA_bulk", "gene_gene", "COL1A1", "COL3A1", "primary_tumours",
    spear(tg["COL1A1"].values, tg["COL3A1"].values))
for tf in TFS:
    if tf not in tg.columns:
        continue
    for c in COLS:
        add("TCGA-BRCA_bulk", "TF_target", tf, c, "primary_tumours",
            spear(tg[tf].values, tg[c].values))
for m in MIRS:
    if m not in tm.columns:
        continue
    for c in COLS:
        add("TCGA-BRCA_bulk", "miRNA_target", m, c, "paired_tumours",
            spear(tm[m].values, tgp[c].values))
if "hsa-miR-101" in tm.columns:
    add("TCGA-BRCA_bulk", "miRNA_target(control)", "hsa-miR-101", "EZH2", "paired_tumours",
        spear(tm["hsa-miR-101"].values, tgp["EZH2"].values))

# TCGA partial correlation given a 4-marker CAF score (DCN, LUM, FAP, THY1)
caf_genes = [g for g in ["DCN", "LUM", "FAP", "THY1"] if g in tg.columns]
print("CAF marker genes used:", caf_genes)


def partial_spear(x, y, z):
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(z), axis=1)
    x, y, z = x[m], y[m], z[m]
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    rz = np.column_stack([stats.rankdata(z[:, j]) for j in range(z.shape[1])])
    rz = np.column_stack([np.ones(len(rz)), rz])
    bx = np.linalg.lstsq(rz, rx, rcond=None)[0]
    by = np.linalg.lstsq(rz, ry, rcond=None)[0]
    ex = rx - rz @ bx
    ey = ry - rz @ by
    r, p = stats.pearsonr(ex, ey)
    return dict(n=int(m.sum()), rho=float(r), p=float(p), r_pearson=np.nan, p_pearson=np.nan)


if caf_genes:
    Zt = tg[caf_genes].values
    Ztp = tgp[caf_genes].values
    for tf in TFS:
        if tf not in tg.columns:
            continue
        for c in COLS:
            add("TCGA-BRCA_bulk_CAFadj", "TF_target", tf, c, "primary_tumours,partial|DCN+LUM+FAP+THY1",
                partial_spear(tg[tf].values, tg[c].values, Zt))
    for m in MIRS:
        if m not in tm.columns:
            continue
        for c in COLS:
            add("TCGA-BRCA_bulk_CAFadj", "miRNA_target", m, c, "paired_tumours,partial|DCN+LUM+FAP+THY1",
                partial_spear(tm[m].values, tgp[c].values, Ztp))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(RES, "celllines_B_correlations.csv"), index=False)
print(df.head(20).to_string())

# ------------------------------------------------------------------ C: absolute levels
lin = pd.read_csv(os.path.join(RES, "celllines_C_lineage_collagen.csv"))
pct_all = pd.read_csv(os.path.join(RES, "celllines_C_percentile_all_models.csv"), index_col=0)
tp = pd.read_csv(os.path.join(RES, "_tcga_collagen_pct.tsv"), sep="\t", index_col=0)

crows = []
for g in COLS:
    v = expr.loc[cancer_ids, g].values
    tpm = 2 ** v - 1
    crows.append(dict(
        gene=g, group="DepMap breast cancer cell lines", n=len(v),
        median_log2TPM1=float(np.median(v)), mean_log2TPM1=float(np.mean(v)),
        q25_log2TPM1=float(np.percentile(v, 25)), q75_log2TPM1=float(np.percentile(v, 75)),
        max_log2TPM1=float(np.max(v)),
        median_TPM=float(np.median(tpm)), max_TPM=float(np.max(tpm)),
        frac_TPM_gt1=float(np.mean(tpm > 1)), frac_TPM_gt10=float(np.mean(tpm > 10)),
        frac_TPM_gt100=float(np.mean(tpm > 100)),
        median_within_sample_pctile=float(np.median(pct_all.loc[cancer_ids, g]))))
    tv = tp[g].values
    crows.append(dict(
        gene=g, group="TCGA-BRCA primary tumours (bulk)", n=len(tv),
        median_log2TPM1=np.nan, mean_log2TPM1=np.nan, q25_log2TPM1=np.nan,
        q75_log2TPM1=np.nan, max_log2TPM1=np.nan, median_TPM=np.nan, max_TPM=np.nan,
        frac_TPM_gt1=np.nan, frac_TPM_gt10=np.nan, frac_TPM_gt100=np.nan,
        median_within_sample_pctile=float(np.median(tv))))
    # other lineages for scale
    for L in ["Fibroblast", "Soft Tissue", "Bone", "Skin", "Lung", "Pancreas"]:
        sel = pct_all.index[pct_all["lineage"] == L]
        if len(sel) < 3:
            continue
        crows.append(dict(gene=g, group=f"DepMap {L} lines", n=len(sel),
                          median_log2TPM1=np.nan, mean_log2TPM1=np.nan, q25_log2TPM1=np.nan,
                          q75_log2TPM1=np.nan, max_log2TPM1=np.nan, median_TPM=np.nan,
                          max_TPM=np.nan, frac_TPM_gt1=np.nan, frac_TPM_gt10=np.nan,
                          frac_TPM_gt100=np.nan,
                          median_within_sample_pctile=float(np.median(pct_all.loc[sel, g]))))
pd.DataFrame(crows).to_csv(os.path.join(RES, "celllines_C_absolute_expression.csv"), index=False)
print(pd.DataFrame(crows).to_string())

tab = expr.loc[cancer_ids, ["COL1A1", "COL3A1", "COL1A2", "FN1", "SPARC", "VIM", "EPCAM", "CDH1"]].copy()
tab.insert(0, "cell_line", models.loc[tab.index, "StrippedCellLineName"])
tab["COL1A1_TPM"] = 2 ** tab["COL1A1"] - 1
tab["COL3A1_TPM"] = 2 ** tab["COL3A1"] - 1
tab["COL1A1_pctile"] = pct_all.loc[tab.index, "COL1A1"]
tab["COL3A1_pctile"] = pct_all.loc[tab.index, "COL3A1"]
tab.sort_values("COL1A1", ascending=False).to_csv(
    os.path.join(RES, "celllines_C_per_line_collagen.csv"))
print("DONE")
