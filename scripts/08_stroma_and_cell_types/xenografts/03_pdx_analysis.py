#!/usr/bin/env python3
"""
PDX natural experiment.
NCI Patient-Derived Models Repository (PDMR), cBioPortal study pancan_pdmr_2025,
profile pancan_pdmr_2025_rna_seq_v2_mrna (RSEM, human genome).
Breast-cancer samples only. PASSAGE == 'Originator' are the patient specimens
(human tumour WITH human stroma); PASSAGE == P0..Pn are PDX (human carcinoma,
mouse stroma; mouse reads do not map to the human transcriptome).
TCGA-BRCA primary tumours are the external bulk reference.
Standardisation: within-sample percentile rank of each gene among a fixed
reference panel of genes, which is invariant to the compositional
renormalisation that removing the stromal transcriptome causes.
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats

BASE="/path/to/revision"
CACHE=f"{BASE}/cache/v6/pdx"; RES=f"{BASE}/results/v6"
os.makedirs(RES, exist_ok=True)

TARGET_ORDER = ["COL1A1","COL3A1","FN1","PDGFRB","CXCL12","POSTN","MMP14","PLAU",
                "MET","STAT5A","CCND2","MYBL2","E2F1","E2F3","EZH2","DNMT1","BRCA1",
                "GATA3","ESR1","JUN","EGR2","SREBF1"]
EXTRA = ["NFKB1","ETS1","EPCAM","KRT8","KRT18","KRT19","ERBB2","PTPRC","PECAM1",
         "ACTA2","DCN","LUM","VIM","MKI67","COL1A2","THY1","FAP","PDGFRA"]

sym2ent = json.load(open(f"{CACHE}/gene_symbol_entrez.json"))
ent2sym = {}
for s,e in sym2ent.items(): ent2sym.setdefault(e, s)
meta = json.load(open("/path/to/scratch/pdmr_breast_meta.json"))

long = pd.read_csv(f"{CACHE}/pdmr_breast_expr_long.tsv", sep="\t")
W = long.pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="mean")
W.columns = [ent2sym.get(c, str(c)) for c in W.columns]
W = W.loc[:, ~pd.Index(W.columns).duplicated()]
print("PDMR matrix", W.shape)

ref_panel = json.load(open("/path/to/scratch/ref_panel_entrez.json"))
REF = sorted(set(ref_panel.keys()) & set(W.columns))
TARGETS = [g for g in TARGET_ORDER+EXTRA if g in W.columns]
REF = [g for g in REF if g not in TARGETS]
print("reference panel genes:", len(REF), "targets:", len(TARGETS))

def pct_rank(df, ref, targets):
    """within-sample percentile of each target among the reference panel."""
    R = df[ref].to_numpy(dtype=float)
    out = {}
    for g in targets:
        v = df[g].to_numpy(dtype=float)
        pr = np.full(len(v), np.nan)
        for i in range(len(v)):
            r = R[i]; ok = np.isfinite(r)
            if np.isfinite(v[i]) and ok.sum() > 100:
                pr[i] = 100.0*(np.mean(r[ok] < v[i]) + 0.5*np.mean(r[ok] == v[i]))
        out[g] = pr
    return pd.DataFrame(out, index=df.index)

P = pct_rank(W, REF, TARGETS)
info = pd.DataFrame({s: dict(passage=meta[s].get("PASSAGE"), model=meta[s].get("MODEL_ID"),
                             culture=meta[s].get("CULTURE_ORIGIN"),
                             human_read=meta[s].get("HUMAN_READ"),
                             onco=meta[s].get("ONCOTREE_CODE"))
                     for s in P.index}).T
def grp(p):
    if p == "Originator": return "PDMR_patient_originator"
    if isinstance(p,str) and p.startswith("P") and p[1:].isdigit(): return "PDMR_PDX"
    if p in ("PDOrg","PDC"): return "PDMR_invitro"
    return "other"
info["group"] = info.passage.map(grp)
P = P.join(info)
P.to_csv(f"{RES}/pdx_percentile_ranks_persample.csv")
print(P.group.value_counts().to_dict())

# ---- TCGA reference -------------------------------------------------------
tc = pd.read_csv(f"{CACHE}/tcga_primary_expr_panel.tsv", sep="\t", index_col=0)
print("TCGA matrix", tc.shape)
REF_T = [g for g in REF if g in tc.columns]
TG_T  = [g for g in TARGETS if g in tc.columns]
PT = pct_rank(tc, REF_T, TG_T)
PT["group"] = "TCGA_primary_tumour"
PT.to_csv(f"{RES}/pdx_percentile_ranks_tcga.csv")

# ---- cohort-level comparison ---------------------------------------------
rows=[]
pdx = P[P.group=="PDMR_PDX"]; org = P[P.group=="PDMR_patient_originator"]; vit = P[P.group=="PDMR_invitro"]
for g in TARGET_ORDER+EXTRA:
    if g not in P.columns: continue
    a = pdx[g].dropna().values; b = PT[g].dropna().values if g in PT.columns else np.array([])
    o = org[g].dropna().values; v = vit[g].dropna().values
    def mw(x,y):
        if len(x)<3 or len(y)<3: return np.nan, np.nan
        u,p = stats.mannwhitneyu(x,y,alternative="two-sided")
        return u/(len(x)*len(y)), p
    auc_t,p_t = mw(a,b); auc_o,p_o = mw(a,o)
    rows.append(dict(gene=g,
        n_pdx=len(a), n_tcga=len(b), n_originator=len(o), n_invitro=len(v),
        pctrank_TCGA_primary=np.median(b) if len(b) else np.nan,
        pctrank_PDMR_originator=np.median(o) if len(o) else np.nan,
        pctrank_PDMR_PDX=np.median(a) if len(a) else np.nan,
        pctrank_PDMR_invitro=np.median(v) if len(v) else np.nan,
        delta_PDX_minus_TCGA=(np.median(a)-np.median(b)) if len(a) and len(b) else np.nan,
        delta_PDX_minus_originator=(np.median(a)-np.median(o)) if len(a) and len(o) else np.nan,
        auc_PDX_vs_TCGA=auc_t, p_PDX_vs_TCGA=p_t,
        auc_PDX_vs_originator=auc_o, p_PDX_vs_originator=p_o))
D = pd.DataFrame(rows)
for c in ["p_PDX_vs_TCGA","p_PDX_vs_originator"]:
    v = D[c].values; ok=np.isfinite(v); q=np.full(len(v),np.nan)
    if ok.sum():
        from statsmodels.stats.multitest import multipletests
        q[ok]=multipletests(v[ok],method="fdr_bh")[1]
    D[c.replace("p_","q_")]=q
D = D.sort_values("delta_PDX_minus_TCGA")
D.to_csv(f"{RES}/pdx_cohort_rank_change.csv", index=False)
print(D[["gene","pctrank_TCGA_primary","pctrank_PDMR_originator","pctrank_PDMR_PDX",
         "delta_PDX_minus_TCGA","auc_PDX_vs_TCGA","q_PDX_vs_TCGA"]].round(3).to_string(index=False))

# ---- matched model-level analysis -----------------------------------------
mrows=[]
models = sorted(set(org.model) & set(pdx.model))
print("matched models:", len(models), models)
for g in TARGET_ORDER+EXTRA:
    if g not in P.columns: continue
    d=[]
    for m in models:
        o = org[(org.model==m)][g].dropna(); a = pdx[(pdx.model==m)][g].dropna()
        if len(o) and len(a): d.append((m, float(o.mean()), float(a.mean())))
    if len(d) < 5: continue
    dd = pd.DataFrame(d, columns=["model","originator","pdx"])
    diff = dd.pdx - dd.originator
    w,pw = stats.wilcoxon(dd.pdx, dd.originator)
    t,pt = stats.ttest_rel(dd.pdx, dd.originator)
    mrows.append(dict(gene=g, n_models=len(dd),
        median_pctrank_originator=float(dd.originator.median()),
        median_pctrank_pdx=float(dd.pdx.median()),
        median_paired_delta=float(diff.median()), mean_paired_delta=float(diff.mean()),
        n_models_down=int((diff<0).sum()), n_models_up=int((diff>0).sum()),
        wilcoxon_p=pw, ttest_p=pt))
M = pd.DataFrame(mrows).sort_values("median_paired_delta")
from statsmodels.stats.multitest import multipletests
M["wilcoxon_q"]=multipletests(M.wilcoxon_p.values, method="fdr_bh")[1]
M.to_csv(f"{RES}/pdx_matched_pairs_rank_change.csv", index=False)
print(M.round(4).to_string(index=False))

# per-model detail for the collagen genes
det=[]
for m in models:
    for g in ["COL1A1","COL3A1","FN1","POSTN","PDGFRB","CXCL12","EPCAM","E2F1","EZH2","MYBL2","CCND2"]:
        if g not in P.columns: continue
        o = org[org.model==m][g].dropna(); a = pdx[pdx.model==m][g].dropna()
        if len(o) and len(a):
            det.append(dict(model=m, gene=g, n_pdx_samples=len(a),
                            pctrank_originator=float(o.mean()), pctrank_pdx=float(a.mean()),
                            delta=float(a.mean()-o.mean())))
pd.DataFrame(det).to_csv(f"{RES}/pdx_matched_permodel_detail.csv", index=False)

# human read fraction QC
hr = info.dropna(subset=["human_read"]).copy()
hr["human_read"]=hr.human_read.astype(float)
hr.groupby("group").human_read.describe().to_csv(f"{RES}/pdx_human_read_fraction_qc.csv")
print(hr.groupby("group").human_read.describe())
print("DONE")
