#!/usr/bin/env python3
"""
celllines_03_controls.py

Everything needed to interpret celllines_02_corr.py honestly:

 1. Assay-quality positive controls for the CCLE nanoString miRNA panel in the same
    50 breast lines (miR-200/141/205 -> ZEB1/ZEB2, let-7 -> HMGA2, miR-21 -> PDCD4).
    If these work and miR-29 -> collagen does not, the miR-29 null is a real null.
 2. Dynamic range of miR-29a/b/c and of COL1A1/COL3A1 across the 50 lines.
 3. Network-level test: over ALL canonical targets of each miRNA, what fraction have
    the expected negative sign in cell lines vs in TCGA bulk, against a permutation null.
 4. Epithelial-mesenchymal confounding INSIDE the cell-line panel: partial Spearman of
    TF -> collagen and miR-29 -> collagen given a mesenchymal score, plus the same
    correlations restricted to lines with detectable collagen.
 5. Absolute TPM by lineage (fibroblast lines as the positive control for collagen).
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")
RNG = np.random.default_rng(20260909)

expr = pd.read_csv(os.path.join(RES, "celllines_expr_breast.tsv.gz"), sep="\t", index_col=0)
mir = pd.read_csv(os.path.join(RES, "celllines_mirna_breast.tsv"), sep="\t", index_col=0)
mir = np.log2(mir + 1.0)
models = pd.read_csv(os.path.join(RES, "celllines_breast_models.csv")).set_index("ModelID")
ids = [i for i in mir.index if i in expr.index]
E = expr.loc[ids]
M = mir.loc[ids]
print("matched breast lines:", len(ids))


def sp(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def partial_spear(x, y, Z):
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(Z), axis=1)
    x, y, Z = x[m], y[m], Z[m]
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rz = np.column_stack([np.ones(m.sum())] + [stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    ex = rx - rz @ np.linalg.lstsq(rz, rx, rcond=None)[0]
    ey = ry - rz @ np.linalg.lstsq(rz, ry, rcond=None)[0]
    r, p = stats.pearsonr(ex, ey)
    return float(r), float(p), int(m.sum())


rows = []

# ---------------------------------------------------------------- 1. positive controls
CTRL = [("hsa-miR-200c", "ZEB1", -1), ("hsa-miR-200c", "ZEB2", -1),
        ("hsa-miR-200b", "ZEB1", -1), ("hsa-miR-200b", "ZEB2", -1),
        ("hsa-miR-200a", "ZEB1", -1), ("hsa-miR-141", "ZEB2", -1),
        ("hsa-miR-205", "ZEB1", -1), ("hsa-miR-203", "ZEB2", -1),
        ("hsa-let-7a", "HMGA2", -1), ("hsa-let-7g", "HMGA2", -1),
        ("hsa-miR-21", "PDCD4", -1), ("hsa-miR-200c", "CDH1", +1),
        ("hsa-miR-200c", "VIM", -1)]
for m, g, exp_sign in CTRL:
    if m not in M.columns or g not in E.columns:
        rows.append(dict(block="1_positive_control", regulator=m, target=g,
                         expected_sign=exp_sign, n=np.nan, rho=np.nan, p=np.nan,
                         note="probe or gene missing"))
        continue
    r, p, n = sp(M[m].values, E[g].values)
    rows.append(dict(block="1_positive_control", regulator=m, target=g,
                     expected_sign=exp_sign, n=n, rho=r, p=p,
                     note="sign_as_expected" if np.sign(r) == exp_sign else "sign_WRONG"))

# ---------------------------------------------------------------- 2. dynamic range
rng_rows = []
for m in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-miR-200c", "hsa-miR-21",
          "hsa-let-7b", "hsa-miR-101"]:
    if m not in M.columns:
        continue
    v = M[m].values
    rng_rows.append(dict(feature=m, platform="CCLE nanoString log2(count+1)", n=len(v),
                         min=float(v.min()), q25=float(np.percentile(v, 25)),
                         median=float(np.median(v)), q75=float(np.percentile(v, 75)),
                         max=float(v.max()), IQR=float(np.percentile(v, 75) - np.percentile(v, 25)),
                         sd=float(v.std(ddof=1))))
for g in ["COL1A1", "COL3A1", "ZEB1", "ZEB2", "CDH1", "VIM", "EZH2", "PDCD4", "HMGA2"]:
    if g not in E.columns:
        continue
    v = E[g].values
    rng_rows.append(dict(feature=g, platform="DepMap log2(TPM+1)", n=len(v),
                         min=float(v.min()), q25=float(np.percentile(v, 25)),
                         median=float(np.median(v)), q75=float(np.percentile(v, 75)),
                         max=float(v.max()), IQR=float(np.percentile(v, 75) - np.percentile(v, 25)),
                         sd=float(v.std(ddof=1))))
pd.DataFrame(rng_rows).to_csv(os.path.join(RES, "celllines_B_dynamic_range.csv"), index=False)
print(pd.DataFrame(rng_rows).to_string())

# ---------------------------------------------------------------- 3. network-level sign test
edges = pd.read_csv(os.path.join(DATA, "canonical_edges.tsv"), sep="\t")
mt = edges[edges.edge_type == "miRNA_target"]
focus = ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-miR-200c", "hsa-miR-21",
         "hsa-let-7b", "hsa-miR-101", "hsa-miR-34a", "hsa-miR-124"]
avail_mirs = [m for m in M.columns]
net_rows = []
N_PERM = 2000
for m in focus:
    if m not in M.columns:
        continue
    tg = sorted(set(mt.loc[mt.source == m, "target"]) & set(E.columns))
    if len(tg) < 5:
        continue
    rhos = np.array([sp(M[m].values, E[g].values)[0] for g in tg])
    frac_neg = float(np.mean(rhos < 0))
    mean_rho = float(np.nanmean(rhos))
    # permutation null: same target set, but the miRNA vector is replaced by a random
    # OTHER miRNA measured on the same panel (preserves target-set gene properties)
    others = [x for x in avail_mirs if x != m]
    null_frac = np.empty(N_PERM)
    null_mean = np.empty(N_PERM)
    for i in range(N_PERM):
        mm = others[RNG.integers(len(others))]
        rr = np.array([sp(M[mm].values, E[g].values)[0] for g in tg])
        null_frac[i] = np.mean(rr < 0)
        null_mean[i] = np.nanmean(rr)
    p_frac = (np.sum(null_frac >= frac_neg) + 1) / (N_PERM + 1)
    p_mean = (np.sum(null_mean <= mean_rho) + 1) / (N_PERM + 1)
    net_rows.append(dict(dataset="cell lines (n=%d)" % len(ids), miRNA=m, n_targets=len(tg),
                         frac_negative=frac_neg, mean_rho=mean_rho,
                         null_frac_negative=float(null_frac.mean()),
                         null_mean_rho=float(null_mean.mean()),
                         p_frac_neg=float(p_frac), p_mean_rho=float(p_mean),
                         n_perm=N_PERM))
    print("net", m, len(tg), round(frac_neg, 3), round(mean_rho, 4), p_frac, p_mean, flush=True)
pd.DataFrame(net_rows).to_csv(os.path.join(RES, "celllines_B_network_sign_test.csv"), index=False)

# ---------------------------------------------------------------- 4. EMT confounding in lines
MES = [g for g in ["VIM", "ZEB1", "ZEB2", "SNAI2", "TWIST1", "CDH2", "SPARC", "THY1"] if g in E.columns]
EPI = [g for g in ["CDH1", "EPCAM", "KRT8", "KRT18", "ESRP1"] if g in E.columns]
z = lambda df: (df - df.mean()) / df.std(ddof=1)
mes_score = z(E[MES]).mean(axis=1)
epi_score = z(E[EPI]).mean(axis=1)
print("mesenchymal score genes:", MES, " epithelial:", EPI)
r_mc1, p_mc1, _ = sp(mes_score.values, E["COL1A1"].values)
r_mc3, p_mc3, _ = sp(mes_score.values, E["COL3A1"].values)
print(f"mes_score ~ COL1A1 rho={r_mc1:.3f} p={p_mc1:.2g}; ~ COL3A1 rho={r_mc3:.3f} p={p_mc3:.2g}")

Eall = expr  # all 70 breast models for the TF analysis
mes_all = z(Eall[MES]).mean(axis=1)
TFS = ["ETS1", "NFKB1", "RELA", "SP1", "MYB", "MRTFA", "STAT6", "TFAP2A"]
conf_rows = []
for tf in TFS:
    if tf not in Eall.columns:
        conf_rows.append(dict(block="4_EMT_adjust", regulator=tf, target="-", n=np.nan,
                              rho_raw=np.nan, p_raw=np.nan, rho_partial_MES=np.nan,
                              p_partial_MES=np.nan, note="gene symbol absent"))
        continue
    for c in ["COL1A1", "COL3A1"]:
        r0, p0, n0 = sp(Eall[tf].values, Eall[c].values)
        r1, p1, n1 = partial_spear(Eall[tf].values, Eall[c].values,
                                   mes_all.values.reshape(-1, 1))
        rtm, ptm, _ = sp(Eall[tf].values, mes_all.values)
        conf_rows.append(dict(block="4_EMT_adjust", regulator=tf, target=c, n=n0,
                              rho_raw=r0, p_raw=p0, rho_partial_MES=r1, p_partial_MES=p1,
                              rho_TF_vs_MESscore=rtm, p_TF_vs_MESscore=ptm))
for m in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"]:
    for c in ["COL1A1", "COL3A1"]:
        r0, p0, n0 = sp(M[m].values, E[c].values)
        r1, p1, _ = partial_spear(M[m].values, E[c].values, mes_score.values.reshape(-1, 1))
        rtm, ptm, _ = sp(M[m].values, mes_score.values)
        conf_rows.append(dict(block="4_EMT_adjust", regulator=m, target=c, n=n0,
                              rho_raw=r0, p_raw=p0, rho_partial_MES=r1, p_partial_MES=p1,
                              rho_TF_vs_MESscore=rtm, p_TF_vs_MESscore=ptm))
pd.DataFrame(conf_rows).to_csv(os.path.join(RES, "celllines_B_EMT_adjusted.csv"), index=False)
print(pd.DataFrame(conf_rows).to_string())

# restricted to lines where the collagen is actually detected
det_rows = []
for c, thr in [("COL1A1", 1.0), ("COL3A1", 1.0)]:
    tpm = 2 ** Eall[c] - 1
    keep = Eall.index[tpm > thr]
    for tf in TFS:
        if tf not in Eall.columns:
            continue
        r, p, n = sp(Eall.loc[keep, tf].values, Eall.loc[keep, c].values)
        det_rows.append(dict(block="4b_detected_only", regulator=tf, target=c,
                             threshold_TPM=thr, n=n, rho=r, p=p))
    keep_m = [i for i in keep if i in M.index]
    for m in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"]:
        r, p, n = sp(M.loc[keep_m, m].values, Eall.loc[keep_m, c].values)
        det_rows.append(dict(block="4b_detected_only", regulator=m, target=c,
                             threshold_TPM=thr, n=n, rho=r, p=p))
pd.DataFrame(det_rows).to_csv(os.path.join(RES, "celllines_B_detected_only.csv"), index=False)
print(pd.DataFrame(det_rows).to_string())

pd.DataFrame(rows).to_csv(os.path.join(RES, "celllines_B_positive_controls.csv"), index=False)
print(pd.DataFrame(rows).to_string())

# ---------------------------------------------------------------- 5. absolute TPM by lineage
lin = pd.read_csv(os.path.join(RES, "celllines_C_lineage_collagen.csv"))
for c in ["COL1A1", "COL3A1"]:
    lin[c + "_median_TPM"] = 2 ** lin[c + "_median"] - 1
    lin[c + "_max_TPM"] = 2 ** lin[c + "_max"] - 1
lin.to_csv(os.path.join(RES, "celllines_C_lineage_collagen_TPM.csv"), index=False)
br = lin[lin.lineage == "Breast"].iloc[0]
fb = lin[lin.lineage == "Fibroblast"].iloc[0]
print(f"COL1A1 median TPM: Breast {br.COL1A1_median_TPM:.2f}  Fibroblast {fb.COL1A1_median_TPM:.1f}"
      f"  ratio {fb.COL1A1_median_TPM/br.COL1A1_median_TPM:.0f}x")
print(f"COL3A1 median TPM: Breast {br.COL3A1_median_TPM:.3f} Fibroblast {fb.COL3A1_median_TPM:.1f}"
      f"  ratio {fb.COL3A1_median_TPM/br.COL3A1_median_TPM:.0f}x")
print("DONE")
