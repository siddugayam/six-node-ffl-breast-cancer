#!/usr/bin/env python3
"""
celllines_04_crispr.py -- part D of the stroma-free test.

D1  Is there a dependency-expression relationship?  Spearman( CRISPR gene effect of
    ETS1 / NFKB1 / RELA / SP1 / MYB / MRTFA / STAT6 / TFAP2A , expression of COL1A1 or
    COL3A1 ) across breast cancer lines, and pan-cancer for power.
D2  Are collagen-high lines differentially dependent on any network hub?  Spearman of a
    collagen module score against the dependency of every gene, restricted to the
    canonical-network nodes (primary) and genome-wide (exploratory).  BH-FDR.
D3  Essentiality re-check: how many breast lines is each network TF / collagen essential in.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")
CRISPR = "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data/CRISPRGeneEffect.csv"

expr_br = pd.read_csv(os.path.join(RES, "celllines_expr_breast.tsv.gz"), sep="\t", index_col=0)
models = pd.read_csv(os.path.join(RES, "celllines_breast_models.csv")).set_index("ModelID")
nodes = pd.read_csv(os.path.join(DATA, "canonical_nodes.tsv"), sep="\t")
print("canonical nodes:", nodes.shape, list(nodes.columns))

ce = pd.read_csv(CRISPR, index_col=0, low_memory=False)
ce.columns = [c.split(" (")[0] for c in ce.columns]
print("CRISPRGeneEffect:", ce.shape)
assert not pd.Index(ce.columns).has_duplicates, "duplicated gene symbols in CRISPR matrix"

br_ids = [i for i in expr_br.index if i in ce.index]
print("breast lines with BOTH CRISPR and expression:", len(br_ids))

# full DepMap expression only for the pan-cancer arm
full_expr = pd.read_csv(os.path.join(DATA, "depmap",
                                     "OmicsExpressionProteinCodingGenesTPMLogp1_24Q4.csv"),
                        index_col=0, usecols=None, low_memory=False)
full_expr.columns = [c.split(" (")[0] for c in full_expr.columns]
pan_ids = [i for i in full_expr.index if i in ce.index]
print("pan-cancer lines with BOTH:", len(pan_ids))

TFS = ["ETS1", "NFKB1", "RELA", "SP1", "MYB", "MRTFA", "STAT6", "TFAP2A", "MYC", "COL1A1", "COL3A1"]
COLS = ["COL1A1", "COL3A1"]


def sp(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


# ------------------------------------------------------------------ D1
d1 = []
for g in TFS:
    if g not in ce.columns:
        d1.append(dict(panel="breast", dep_gene=g, expr_gene="-", n=np.nan, rho=np.nan,
                       p=np.nan, note="gene not screened"))
        continue
    for c in COLS:
        r, p, n = sp(ce.loc[br_ids, g].values, expr_br.loc[br_ids, c].values)
        d1.append(dict(panel="breast_cancer_lines", dep_gene=g, expr_gene=c, n=n, rho=r, p=p))
        r2, p2, n2 = sp(ce.loc[pan_ids, g].values, full_expr.loc[pan_ids, c].values)
        d1.append(dict(panel="pan_cancer_lines", dep_gene=g, expr_gene=c, n=n2, rho=r2, p=p2))
    # self-dependency control: does dependency on the gene track its own expression?
    r, p, n = sp(ce.loc[br_ids, g].values, expr_br.loc[br_ids, g].values)
    d1.append(dict(panel="breast_cancer_lines", dep_gene=g, expr_gene=g + "(self)", n=n, rho=r, p=p))
    r2, p2, n2 = sp(ce.loc[pan_ids, g].values, full_expr.loc[pan_ids, g].values)
    d1.append(dict(panel="pan_cancer_lines", dep_gene=g, expr_gene=g + "(self)", n=n2, rho=r2, p=p2))
d1 = pd.DataFrame(d1)
d1["FDR_BH"] = np.nan
msk = d1["p"].notna()
d1.loc[msk, "FDR_BH"] = stats.false_discovery_control(d1.loc[msk, "p"].values, method="bh")
d1.to_csv(os.path.join(RES, "celllines_D_dependency_vs_expression.csv"), index=False)
print(d1.to_string())

# ------------------------------------------------------------------ D2
COLLAGEN_MODULE = ["COL1A1", "COL3A1", "COL1A2", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
                   "COL6A3", "COL4A1", "COL11A1", "FN1", "SPARC", "LOX"]
cm = [g for g in COLLAGEN_MODULE if g in expr_br.columns]
print("collagen module genes used:", cm)
Z = (expr_br[cm] - expr_br[cm].mean()) / expr_br[cm].std(ddof=1)
coll_score = Z.mean(axis=1)
coll_score.to_frame("collagen_module_score").join(
    models[["StrippedCellLineName"]]).to_csv(
        os.path.join(RES, "celllines_D_collagen_module_score.csv"))

X = coll_score.loc[br_ids].values
CE = ce.loc[br_ids]
ok = CE.notna().sum(axis=0) >= 25
CE = CE.loc[:, ok]
print("genes with >=25 non-missing dependency values in breast lines:", CE.shape[1])

rk_x = stats.rankdata(X)
res = []
Y = CE.values
for j, g in enumerate(CE.columns):
    y = Y[:, j]
    m = np.isfinite(y)
    if m.sum() < 25:
        continue
    r, p = stats.spearmanr(X[m], y[m])
    res.append((g, int(m.sum()), float(r), float(p)))
d2 = pd.DataFrame(res, columns=["gene", "n", "rho", "p"])
d2["FDR_BH_genomewide"] = stats.false_discovery_control(d2["p"].values, method="bh")
net_nodes = set(nodes["name"]) if "name" in nodes.columns else set(nodes.iloc[:, 0])
d2["is_network_node"] = d2["gene"].isin(net_nodes)
sub = d2[d2.is_network_node].copy()
sub["FDR_BH_networkonly"] = stats.false_discovery_control(sub["p"].values, method="bh")
d2 = d2.merge(sub[["gene", "FDR_BH_networkonly"]], on="gene", how="left")
d2 = d2.sort_values("p")
d2.to_csv(os.path.join(RES, "celllines_D_collagenscore_vs_dependency.csv"), index=False)
print("n network nodes screened:", int(d2.is_network_node.sum()))
print("genome-wide FDR<0.05:", int((d2.FDR_BH_genomewide < 0.05).sum()),
      "| network-only FDR<0.05:", int((d2.FDR_BH_networkonly < 0.05).sum()))
print(d2.head(25).to_string())
print("--- network nodes, top 20 ---")
print(d2[d2.is_network_node].head(20).to_string())

# ------------------------------------------------------------------ D3 essentiality
d3 = []
for g in sorted(set(list(nodes["name"]) if "name" in nodes.columns else nodes.iloc[:, 0])):
    if g not in ce.columns:
        continue
    v = ce.loc[br_ids, g].dropna()
    if len(v) < 25:
        continue
    d3.append(dict(gene=g, n_lines=len(v), median_gene_effect=float(v.median()),
                   n_essential_lt_minus0p5=int((v < -0.5).sum()),
                   n_essential_lt_minus1p0=int((v < -1.0).sum()),
                   frac_essential_lt_minus0p5=float((v < -0.5).mean())))
d3 = pd.DataFrame(d3).sort_values("median_gene_effect")
d3.to_csv(os.path.join(RES, "celllines_D_network_essentiality_breast.csv"), index=False)
print(d3.head(15).to_string())
print(d3[d3.gene.isin(["NFKB1", "RELA", "SP1", "ETS1", "MYB", "STAT6", "TFAP2A",
                       "COL1A1", "COL3A1", "MYC", "TP53"])].to_string())
print("DONE")
