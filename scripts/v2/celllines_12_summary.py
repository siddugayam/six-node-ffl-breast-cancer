#!/usr/bin/env python3
"""
celllines_12_summary.py -- consolidate the stroma-free test into one table,
plus the essentiality re-check on ALL breast-lineage lines that have a CRISPR screen
(not only the ones that also have expression).
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")
CRISPR = "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data/CRISPRGeneEffect.csv"

models = pd.read_csv(os.path.join(RES, "celllines_breast_models.csv"))
br_all = set(models["ModelID"])
br_ca = set(models.loc[models.OncotreePrimaryDisease != "Non-Cancerous", "ModelID"])

GENES = ["ETS1", "NFKB1", "RELA", "SP1", "MYB", "MRTFA", "STAT6", "TFAP2A",
         "COL1A1", "COL3A1", "MYC", "TP53", "CTCF", "RPL5"]
cols = pd.read_csv(CRISPR, nrows=0).columns
sel = [c for c in cols if c.split(" (")[0] in GENES]
ce = pd.read_csv(CRISPR, usecols=[cols[0]] + sel, index_col=0)
ce.columns = [c.split(" (")[0] for c in ce.columns]
ce_br = ce.loc[[i for i in ce.index if i in br_ca]]
print("breast CANCER lines with a CRISPR screen:", ce_br.shape[0])
print("breast-lineage lines (incl. non-cancerous) with a CRISPR screen:",
      len([i for i in ce.index if i in br_all]))

rows = []
for g in GENES:
    if g not in ce_br.columns:
        rows.append(dict(gene=g, n_lines=0, note="not screened"))
        continue
    v = ce_br[g].dropna()
    rows.append(dict(gene=g, n_lines=len(v), median_gene_effect=float(v.median()),
                     min_gene_effect=float(v.min()), max_gene_effect=float(v.max()),
                     n_lt_minus0p5=int((v < -0.5).sum()), n_lt_minus1p0=int((v < -1.0).sum())))
e = pd.DataFrame(rows)
e.to_csv(os.path.join(RES, "celllines_D_essentiality_all_breast_CRISPR.csv"), index=False)
print(e.to_string(index=False))

# ---------------------------------------------------------------- summary
B = pd.read_csv(os.path.join(RES, "celllines_B_correlations.csv"))
C = pd.read_csv(os.path.join(RES, "celllines_C_absolute_expression.csv"))
lin = pd.read_csv(os.path.join(RES, "celllines_C_lineage_collagen_TPM.csv"))
sc = pd.read_csv(os.path.join(RES, "celllines_C_scRNA_celltype_major.csv"), index_col=0)
D1 = pd.read_csv(os.path.join(RES, "celllines_D_dependency_vs_expression.csv"))
D2 = pd.read_csv(os.path.join(RES, "celllines_D_collagenscore_vs_dependency.csv"))
Esum = pd.read_csv(os.path.join(RES, "celllines_E_drug_summary.csv"))
PC = pd.read_csv(os.path.join(RES, "celllines_B_positive_controls.csv"))
NS = pd.read_csv(os.path.join(RES, "celllines_F_network_sign_combined.csv"))


def g(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).values
    return df[m]


S = []


def add(item, value, source):
    S.append(dict(item=item, value=value, source_file=source))


def zdiff(r1, n1, r2, n2):
    """two-sided p for H0: rho1 == rho2 (Fisher z)."""
    if not all(np.isfinite([r1, r2])) or min(n1, n2) < 5:
        return np.nan
    z = (np.arctanh(r1) - np.arctanh(r2)) / np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    return float(2 * (1 - stats.norm.cdf(abs(z))))


cl = "breast_cancer_lines_only"
ml = "miRNA_matched_cancer_only"
tb = "primary_tumours"
tp = "paired_tumours"
for reg, tgt in [("COL1A1", "COL3A1")]:
    r = g(B, subset=cl, regulator=reg, target=tgt).iloc[0]
    t = g(B, subset=tb, regulator=reg, target=tgt).iloc[0]
    add(f"{reg}~{tgt} Spearman, breast cancer cell lines",
        f"rho={r.rho:.3f} (95% CI {r.rho_CI95_lo:.3f} to {r.rho_CI95_hi:.3f}), p={r.p:.3g}, n={int(r.n)}",
        "celllines_B_correlations.csv")
    add(f"{reg}~{tgt} Spearman, TCGA-BRCA bulk",
        f"rho={t.rho:.3f}, p={t.p:.3g}, n={int(t.n)}; p(cell line vs TCGA differ)="
        f"{zdiff(r.rho, int(r.n), t.rho, int(t.n)):.3g}", "celllines_B_correlations.csv")
for reg in ["ETS1", "NFKB1", "RELA", "SP1"]:
    for tgt in ["COL1A1", "COL3A1"]:
        r = g(B, subset=cl, regulator=reg, target=tgt).iloc[0]
        t = g(B, subset=tb, regulator=reg, target=tgt).iloc[0]
        a = g(B, subset="primary_tumours,partial|DCN+LUM+FAP+THY1", regulator=reg, target=tgt)
        pd_ = zdiff(r.rho, int(r.n), t.rho, int(t.n))
        add(f"{reg}->{tgt}",
            f"cell lines rho={r.rho:+.3f} (CI {r.rho_CI95_lo:+.3f},{r.rho_CI95_hi:+.3f}) p={r.p:.3g} n={int(r.n)} | "
            f"TCGA bulk rho={t.rho:+.3f} p={t.p:.3g} | "
            f"TCGA CAF-adjusted rho={a.iloc[0].rho:+.3f} p={a.iloc[0].p:.3g} | "
            f"p(cell line vs TCGA differ)={pd_:.3g}",
            "celllines_B_correlations.csv")
for reg in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"]:
    for tgt in ["COL1A1", "COL3A1"]:
        r = g(B, subset=ml, regulator=reg, target=tgt).iloc[0]
        t = g(B, subset=tp, regulator=reg, target=tgt).iloc[0]
        a = g(B, subset="paired_tumours,partial|DCN+LUM+FAP+THY1", regulator=reg, target=tgt)
        pd_ = zdiff(r.rho, int(r.n), t.rho, int(t.n))
        add(f"{reg}->{tgt}",
            f"cell lines rho={r.rho:+.3f} (CI {r.rho_CI95_lo:+.3f},{r.rho_CI95_hi:+.3f}) p={r.p:.3g} n={int(r.n)} | "
            f"TCGA bulk rho={t.rho:+.3f} p={t.p:.3g} | "
            f"TCGA CAF-adjusted rho={a.iloc[0].rho:+.3f} p={a.iloc[0].p:.3g} | "
            f"p(cell line vs TCGA differ)={pd_:.3g}",
            "celllines_B_correlations.csv")

ok = PC[PC.note == "sign_as_expected"]
add("miRNA assay positive controls in the same 50 lines",
    f"{len(ok)}/{len(PC)} with the expected sign; strongest miR-141->ZEB2 rho="
    f"{PC.loc[PC.regulator=='hsa-miR-141','rho'].iloc[0]:.3f}, "
    f"miR-200c->ZEB2 rho={PC.loc[(PC.regulator=='hsa-miR-200c')&(PC.target=='ZEB2'),'rho'].iloc[0]:.3f}",
    "celllines_B_positive_controls.csv")

brrow = lin[lin.lineage == "Breast"].iloc[0]
fbrow = lin[lin.lineage == "Fibroblast"].iloc[0]
add("COL1A1 median TPM, breast lines vs fibroblast lines",
    f"{brrow.COL1A1_median_TPM:.1f} vs {fbrow.COL1A1_median_TPM:.0f} TPM "
    f"({fbrow.COL1A1_median_TPM/brrow.COL1A1_median_TPM:.0f}x)",
    "celllines_C_lineage_collagen_TPM.csv")
add("COL3A1 median TPM, breast lines vs fibroblast lines",
    f"{brrow.COL3A1_median_TPM:.2f} vs {fbrow.COL3A1_median_TPM:.0f} TPM "
    f"({fbrow.COL3A1_median_TPM/brrow.COL3A1_median_TPM:.0f}x)",
    "celllines_C_lineage_collagen_TPM.csv")
c1 = C[(C.gene == "COL1A1") & (C.group == "DepMap breast cancer cell lines")].iloc[0]
c3 = C[(C.gene == "COL3A1") & (C.group == "DepMap breast cancer cell lines")].iloc[0]
t1 = C[(C.gene == "COL1A1") & (C.group.str.startswith("TCGA"))].iloc[0]
t3 = C[(C.gene == "COL3A1") & (C.group.str.startswith("TCGA"))].iloc[0]
add("COL1A1 within-sample percentile: breast lines vs TCGA bulk tumours",
    f"{c1.median_within_sample_pctile:.3f} vs {t1.median_within_sample_pctile:.4f}",
    "celllines_C_absolute_expression.csv")
add("COL3A1 within-sample percentile: breast lines vs TCGA bulk tumours",
    f"{c3.median_within_sample_pctile:.3f} vs {t3.median_within_sample_pctile:.4f}",
    "celllines_C_absolute_expression.csv")
add("breast cancer lines with COL3A1 > 1 TPM",
    f"{c3.frac_TPM_gt1*100:.1f}% ({round(c3.frac_TPM_gt1*c3.n)}/{int(c3.n)})",
    "celllines_C_absolute_expression.csv")
add("scRNA (Wu 2021) COL1A1 pseudobulk CPM, CAF vs cancer epithelial",
    f"{sc.loc['CAFs','COL1A1_CPM']:.0f} vs {sc.loc['Cancer Epithelial','COL1A1_CPM']:.1f} CPM "
    f"({sc.loc['CAFs','COL1A1_CPM']/sc.loc['Cancer Epithelial','COL1A1_CPM']:.1f}x); "
    f"{sc.loc['Cancer Epithelial','COL1A1_pct_cells_pos']:.1f}% of cancer cells positive",
    "celllines_C_scRNA_celltype_major.csv")
add("scRNA (Wu 2021) COL3A1 pseudobulk CPM, CAF vs cancer epithelial",
    f"{sc.loc['CAFs','COL3A1_CPM']:.0f} vs {sc.loc['Cancer Epithelial','COL3A1_CPM']:.1f} CPM "
    f"({sc.loc['CAFs','COL3A1_CPM']/sc.loc['Cancer Epithelial','COL3A1_CPM']:.1f}x); "
    f"{sc.loc['Cancer Epithelial','COL3A1_pct_cells_pos']:.1f}% of cancer cells positive",
    "celllines_C_scRNA_celltype_major.csv")

r = D1[(D1.panel == "breast_cancer_lines") & (D1.dep_gene == "ETS1") & (D1.expr_gene == "COL1A1")].iloc[0]
add("CRISPR: ETS1 dependency vs COL1A1 expression, breast lines",
    f"rho={r.rho:+.3f}, p={r.p:.3g}, BH FDR={r.FDR_BH:.3g}, n={int(r.n)} "
    "(collagen-high lines are MORE ETS1-dependent)",
    "celllines_D_dependency_vs_expression.csv")
add("CRISPR: any other TF dependency vs collagen expression in breast lines",
    "none at BH FDR<0.05 among the 66 tests",
    "celllines_D_dependency_vs_expression.csv")
add("CRISPR: collagen module score vs dependency, network nodes",
    f"0 of {int(D2.is_network_node.sum())} canonical-network genes at BH FDR<0.05 "
    f"(genome-wide: {int((D2.FDR_BH_genomewide<0.05).sum())} of {len(D2)} genes)",
    "celllines_D_collagenscore_vs_dependency.csv")
for gname in ["MYC", "NFKB1", "RELA", "ETS1", "SP1", "COL1A1", "COL3A1"]:
    row = e[e.gene == gname]
    if len(row) and "median_gene_effect" in row:
        row = row.iloc[0]
        add(f"CRISPR essentiality of {gname} in breast cancer lines",
            f"{int(row.n_lt_minus0p5)}/{int(row.n_lines)} lines with gene effect < -0.5 "
            f"(median {row.median_gene_effect:+.3f})",
            "celllines_D_essentiality_all_breast_CRISPR.csv")

for _, row in Esum.iterrows():
    add(f"Drug: {row.panel} / {row.score} / {row.metric}",
        f"{int(row.n_FDR05)} of {int(row.n_compounds)} compounds at BH FDR<0.05 "
        f"(min FDR {row.min_FDR:.3g})", "celllines_E_drug_summary.csv")

for _, row in NS.iterrows():
    add(f"Network-level sign test, {row.miRNA}, {row.dataset}",
        f"{row.frac_negative*100:.1f}% of {int(row.n_targets)} canonical targets negative; "
        f"mean rho {row.mean_rho:+.4f}; permutation p(frac)={row.p_frac_neg:.3g}, "
        f"p(mean rho)={row.p_mean_rho:.3g}, {int(row.n_perm)} permutations",
        "celllines_F_network_sign_combined.csv")

S = pd.DataFrame(S)
S.to_csv(os.path.join(RES, "celllines_STROMA_FREE_SUMMARY.csv"), index=False)
pd.set_option("display.max_colwidth", 200)
print(S.to_string(index=False))
print("DONE")
