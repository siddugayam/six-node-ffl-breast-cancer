#!/usr/bin/env python3
"""
Spatial transcriptomics test of the compartment argument.
Datasets: 10x Genomics public Visium human breast cancer sections
          + Wu et al. 2021 (Zenodo 4739739) Visium sections if present.
For each section: assign spots to stromal / epithelial compartments using
marker scores that contain NO collagen or ECM outcome genes, then compute
Spearman correlation of COL1A1/COL3A1/FN1 with each TF (a) across all spots
and (b) within each compartment separately.
"""
import os, json, glob, warnings
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats, sparse
warnings.filterwarnings("ignore")

BASE = "/path/to/revision"
CACHE = f"{BASE}/cache/v6"
RES   = f"{BASE}/results/v6"
FIG   = f"{BASE}/figures/v6"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# --- marker panels -----------------------------------------------------------
# Epithelial / carcinoma markers (no ECM genes)
EPI = ["EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","ELF3","CLDN4","CLDN3",
       "MUC1","KRT5","KRT14","SLPI","AGR2","TACSTD2"]
# Fibroblast / stromal markers - deliberately EXCLUDES every collagen, FN1,
# POSTN, PDGFRB and CXCL12 so that compartment calling cannot be circular.
STR = ["DCN","LUM","FBLN1","MMP2","THY1","PDGFRA","VCAN","SFRP2","C1S","C1R",
       "CFD","SERPINF1","FBN1","MFAP4","ISLR","PDPN","CALD1"]
IMM = ["PTPRC","CD3D","CD2","CD68","CD14","MS4A1","LYZ","CXCR4"]

OUTCOMES = ["COL1A1","COL3A1","FN1","POSTN","PDGFRB","CXCL12"]
TFS = ["NFKB1","ETS1","E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1",
       "SREBF1","DNMT1","E2F3"]

def load_section(path, libid, count_file="filtered_feature_bc_matrix.h5"):
    a = sc.read_visium(path, count_file=count_file, library_id=libid)
    a.var_names_make_unique()
    a.obs["section"] = libid
    return a

def prep(a):
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    a = a[a.obs.total_counts >= 500, :].copy()
    sc.pp.filter_genes(a, min_cells=3)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a

def score(a, genes, name):
    g = [x for x in genes if x in a.var_names]
    sc.tl.score_genes(a, g, score_name=name, random_state=0)
    return g

def dense(a, g):
    x = a[:, g].X
    return np.asarray(x.todense()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()

def spear(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30: return np.nan, np.nan, int(ok.sum())
    r, p = stats.spearmanr(x[ok], y[ok])
    return r, p, int(ok.sum())

def partial_spear(x, y, z):
    """Spearman partial correlation of x,y given z (rank-based residualisation)."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 30: return np.nan, np.nan
    rx, ry, rz = (stats.rankdata(v[ok]) for v in (x, y, z))
    def resid(v, w):
        b = np.polyfit(w, v, 1); return v - np.polyval(b, w)
    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return r, p

def fisher_z_diff(r1, n1, r2, n2):
    if not all(np.isfinite([r1, r2])) or n1 < 4 or n2 < 4: return np.nan, np.nan
    r1 = np.clip(r1, -0.999999, 0.999999); r2 = np.clip(r2, -0.999999, 0.999999)
    z1, z2 = np.arctanh(r1), np.arctanh(r2)
    se = np.sqrt(1/(n1-3) + 1/(n2-3))
    z = (z1 - z2)/se
    return z, 2*stats.norm.sf(abs(z))

sections = []
# 10x public
for s in ["V1_Breast_Cancer_Block_A_Section_1","V1_Breast_Cancer_Block_A_Section_2",
          "Visium_FFPE_Human_Breast_Cancer"]:
    p = f"{CACHE}/spatial/{s}"
    if os.path.isdir(p): sections.append(("10x_public", s, p, "filtered_feature_bc_matrix.h5"))
# Wu et al 2021
wu = sorted(glob.glob(f"{CACHE}/spatial_wu/*/"))
for p in wu:
    lib = os.path.basename(p.rstrip("/"))
    if os.path.isdir(os.path.join(p, "spatial")):
        cf = [f for f in os.listdir(p) if f.endswith(".h5")]
        if cf: sections.append(("Wu2021", lib, p.rstrip("/"), cf[0]))

print("sections:", [(a,b) for a,b,_,_ in sections], flush=True)

rows_marker, rows_corr, spotframes, qc = [], [], [], []
adatas = {}

for src, lib, path, cf in sections:
    try:
        a = load_section(path, lib, cf)
    except Exception as e:
        print("SKIP", lib, e, flush=True); continue
    n0 = a.n_obs
    a = prep(a)
    ge = score(a, EPI, "epi_score"); gs = score(a, STR, "stroma_score"); gi = score(a, IMM, "imm_score")
    a.obs["delta"] = a.obs["stroma_score"] - a.obs["epi_score"]
    # primary classifier: sign of the difference in marker scores
    a.obs["compartment"] = np.where(a.obs["delta"] > 0, "stromal", "epithelial")
    # tertile classifier (robustness)
    q1, q2 = np.quantile(a.obs["delta"], [1/3, 2/3])
    a.obs["compartment_tertile"] = pd.cut(a.obs["delta"], [-np.inf, q1, q2, np.inf],
                                          labels=["epithelial","intermediate","stromal"])
    qc.append(dict(source=src, section=lib, spots_raw=n0, spots_used=a.n_obs,
                   genes_used=a.n_vars, median_counts=float(np.median(a.obs.total_counts)),
                   median_genes=float(np.median(a.obs.n_genes_by_counts)),
                   n_epi_markers=len(ge), n_stroma_markers=len(gs),
                   n_stromal_spots=int((a.obs.compartment=="stromal").sum()),
                   n_epithelial_spots=int((a.obs.compartment=="epithelial").sum())))

    genes = [g for g in OUTCOMES + TFS + ["EPCAM","PTPRC","ACTA2","DCN"] if g in a.var_names]
    E = pd.DataFrame({g: dense(a, g) for g in genes}, index=a.obs_names)
    E["compartment"] = a.obs["compartment"].values
    E["compartment_tertile"] = a.obs["compartment_tertile"].values
    E["stroma_score"] = a.obs["stroma_score"].values
    E["epi_score"] = a.obs["epi_score"].values
    E["section"] = lib; E["source"] = src
    E["x"] = a.obsm["spatial"][:,0]; E["y"] = a.obsm["spatial"][:,1]
    spotframes.append(E)
    adatas[lib] = a

    st = E.compartment == "stromal"; ep = E.compartment == "epithelial"
    # compartment marker enrichment
    for g in genes:
        v = E[g].values
        d = float(np.mean(v[st]) - np.mean(v[ep]))
        try: u, pu = stats.mannwhitneyu(v[st], v[ep], alternative="two-sided")
        except Exception: u, pu = np.nan, np.nan
        auc = u/(st.sum()*ep.sum()) if np.isfinite(u) else np.nan
        rows_marker.append(dict(source=src, section=lib, gene=g,
            mean_stromal=float(np.mean(v[st])), mean_epithelial=float(np.mean(v[ep])),
            diff_log_stromal_minus_epithelial=d, auc_stromal_vs_epithelial=auc, p_mwu=pu,
            pct_spots_detected=float(np.mean(v > 0))))
    # correlations
    for out in [o for o in OUTCOMES if o in E.columns]:
        for tf in [t for t in TFS if t in E.columns]:
            ra, pa, na = spear(E[out].values, E[tf].values)
            rs, ps, ns = spear(E.loc[st, out].values, E.loc[st, tf].values)
            re_, pe, ne = spear(E.loc[ep, out].values, E.loc[ep, tf].values)
            pr, pp = partial_spear(E[out].values, E[tf].values, E["stroma_score"].values)
            zs, pzs = fisher_z_diff(rs, ns, ra, na)
            ze, pze = fisher_z_diff(re_, ne, ra, na)
            rows_corr.append(dict(source=src, section=lib, outcome=out, tf=tf,
                rho_all=ra, p_all=pa, n_all=na,
                rho_stromal=rs, p_stromal=ps, n_stromal=ns,
                rho_epithelial=re_, p_epithelial=pe, n_epithelial=ne,
                rho_partial_given_stromascore=pr, p_partial=pp,
                z_stromal_vs_all=zs, p_z_stromal_vs_all=pzs,
                z_epithelial_vs_all=ze, p_z_epithelial_vs_all=pze,
                delta_rho_stromal_minus_all=rs-ra if np.isfinite(rs) and np.isfinite(ra) else np.nan,
                delta_rho_epithelial_minus_all=re_-ra if np.isfinite(re_) and np.isfinite(ra) else np.nan))
    print("done", lib, a.shape, flush=True)

pd.DataFrame(qc).to_csv(f"{RES}/spatial_section_qc.csv", index=False)
pd.DataFrame(rows_marker).to_csv(f"{RES}/spatial_compartment_marker_enrichment.csv", index=False)
corr = pd.DataFrame(rows_corr); corr.to_csv(f"{RES}/spatial_tf_collagen_correlation_bysection.csv", index=False)
spots = pd.concat(spotframes); spots.to_csv(f"{RES}/spatial_spot_level_values.csv.gz", index=True, compression="gzip")

# meta-analysis across sections (Fisher z, weighted by n-3)
def meta(g, rcol, ncol):
    d = g[[rcol, ncol]].dropna()
    d = d[d[ncol] > 3]
    if len(d) == 0: return pd.Series({rcol+"_meta": np.nan, rcol+"_meta_p": np.nan})
    z = np.arctanh(np.clip(d[rcol], -0.999999, 0.999999)); w = d[ncol]-3
    zm = np.sum(z*w)/np.sum(w); se = 1/np.sqrt(np.sum(w))
    return pd.Series({rcol+"_meta": np.tanh(zm), rcol+"_meta_p": 2*stats.norm.sf(abs(zm/se))})

mrows = []
for (out, tf), g in corr.groupby(["outcome","tf"]):
    r = dict(outcome=out, tf=tf, n_sections=len(g))
    for rc, nc in [("rho_all","n_all"),("rho_stromal","n_stromal"),("rho_epithelial","n_epithelial")]:
        r.update(meta(g, rc, nc).to_dict())
    r["delta_stromal_minus_all"] = r["rho_stromal_meta"] - r["rho_all_meta"]
    r["delta_epithelial_minus_all"] = r["rho_epithelial_meta"] - r["rho_all_meta"]
    mrows.append(r)
pd.DataFrame(mrows).sort_values(["outcome","tf"]).to_csv(
    f"{RES}/spatial_tf_collagen_correlation_meta.csv", index=False)

mk = pd.DataFrame(rows_marker)
mkm = mk.groupby("gene").agg(n_sections=("section","size"),
        mean_diff=("diff_log_stromal_minus_epithelial","mean"),
        min_diff=("diff_log_stromal_minus_epithelial","min"),
        max_diff=("diff_log_stromal_minus_epithelial","max"),
        mean_auc=("auc_stromal_vs_epithelial","mean"),
        mean_pct_detected=("pct_spots_detected","mean")).reset_index().sort_values("mean_diff", ascending=False)
mkm.to_csv(f"{RES}/spatial_compartment_marker_enrichment_meta.csv", index=False)
print(mkm.to_string(index=False))
json.dump({k: str(v) for k, v in dict(sections_used=[s[1] for s in sections]).items()},
          open(f"{RES}/spatial_sections_used.json","w"))
print("WROTE results")
