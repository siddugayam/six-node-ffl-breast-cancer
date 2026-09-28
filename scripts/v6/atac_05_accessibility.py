#!/usr/bin/env python3
"""(a) Are COL1A1/COL3A1 regulatory elements accessible in TCGA-BRCA tumours, and
does that accessibility track stromal content?

ATAC : TCGA-BRCA ATAC-seq, Corces et al. Science 2018 (GDC ATAC-seq AWG),
       BRCA-specific peak set (215,920 peaks x 74 patients after collapsing
       technical replicates and the one duplicated tissue fragment).
Stroma: ssGSEA/mean-z stromal scores computed from the SAME patients' TCGA-BRCA
       RNA-seq. ESTIMATE stromal signature with every COL* gene removed
       (ESTIMATE_STROMAL_noCOL) and a scRNA-derived CAF-50 signature that is
       collagen-free and network-node-free by construction -> non-circular.
"""
import os, numpy as np, pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
C = os.path.join(ROOT, "cache/v6/atac"); R = os.path.join(ROOT, "results/v6")

meta = pd.read_csv(f"{C}/brca_peak_meta_full.tsv.gz", sep="\t")
S = np.load(f"{C}/brca_sample_matrix.npy")
samp = pd.read_csv(f"{C}/brca_samples.tsv", sep="\t")
tss = pd.read_csv(f"{C}/refseq_select_tss.tsv", sep="\t")
sc = pd.read_csv(f"{ROOT}/results/v5/farmer_tcga_persample_scores.csv")

# collapse the duplicated patient (2 tissue fragments) -> one column per patient
pat = samp.tcga_patient.values
upat = pd.unique(pat)
P = np.zeros((S.shape[0], len(upat)), dtype=np.float32)
for j, p in enumerate(upat):
    P[:, j] = S[:, pat == p].mean(1)
print("peaks x patients:", P.shape)

sc["patient"] = sc["sample"].str.slice(0, 12)
scm = sc[sc.patient.isin(upat)].drop_duplicates("patient").set_index("patient").reindex(upat)
print("patients with matched RNA stromal score:", scm["meanZ_ESTIMATE_STROMAL_noCOL"].notna().sum())

STR = {"ESTIMATE_STROMAL_noCOL": scm["meanZ_ESTIMATE_STROMAL_noCOL"].values,
       "CAF_scRNA_50":           scm["meanZ_CAF_scRNA_50"].values,
       "ESTIMATE_IMMUNE":        scm["meanZ_ESTIMATE_IMMUNE"].values}
ok = ~np.isnan(STR["ESTIMATE_STROMAL_noCOL"])
print("usable patients:", ok.sum())

def spearman_rows(M, v):
    """Spearman of every row of M against vector v (no ties in v)."""
    rv = stats.rankdata(v); rv = (rv - rv.mean()) / rv.std()
    rM = stats.rankdata(M, axis=1)
    rM = (rM - rM.mean(1, keepdims=True)) / rM.std(1, keepdims=True)
    rho = (rM @ rv) / len(rv)
    n = len(rv); df = n - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = rho * np.sqrt(df / (1 - rho ** 2))
    p = 2 * stats.t.sf(np.abs(t), df)
    return rho, p

def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); r = np.empty(n)
    r[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(r, 0, 1)

Pk = P[:, ok]
res = meta.copy()
res["mean_acc"] = Pk.mean(1)
res["acc_pctile"] = stats.rankdata(res.mean_acc) / len(res) * 100
for key, v in STR.items():
    rho, p = spearman_rows(Pk, v[ok])
    res[f"rho_{key}"] = rho
    res[f"p_{key}"] = p
    res[f"fdr_{key}"] = bh(p)
    res[f"rho_pctile_{key}"] = stats.rankdata(rho) / len(rho) * 100
    print(f"{key}: peaks FDR<0.05 = {(res[f'fdr_{key}']<0.05).sum()} "
          f"({(res[f'fdr_{key}']<0.05).mean()*100:.1f}%), median rho {np.median(rho):.3f}")

res.to_csv(f"{C}/peak_stroma_correlation.tsv.gz", sep="\t", index=False)

# ---- locus tables -----------------------------------------------------------
def locus(gene, win=100000):
    t = tss[tss.symbol == gene]
    if not len(t): return None
    t = t.iloc[0]
    sub = res[(res.seqnames == t.chrom) & (res.end > t.tss - win) & (res.start < t.tss + win)].copy()
    mid = (sub.start + sub.end) // 2
    d = mid - t.tss
    sub["dist_tss"] = d if t.strand == "+" else -d
    sub["gene"] = gene
    sub["is_promoter"] = (sub.start <= t.tss + 1000) & (sub.end >= t.tss - 1000)
    return sub

genes = ["COL1A1", "COL3A1", "COL1A2", "EPCAM", "KRT8", "ESR1", "GATA3", "PTPRC",
         "DCN", "LUM", "FAP", "POSTN", "ACTA2", "PDGFRB", "GAPDH", "ACTB",
         "ETS1", "NFKB1", "RELA", "SP1"]
L = pd.concat([locus(g) for g in genes if locus(g) is not None])
L.to_csv(f"{R}/atac_locus_peaks.csv", index=False)

summ = []
for g in genes:
    sub = L[L.gene == g]
    if not len(sub): continue
    pr = sub[sub.is_promoter]
    summ.append(dict(gene=g, n_peaks_100kb=len(sub), n_promoter_peaks=len(pr),
        prom_peak=";".join(pr.name) if len(pr) else "NONE",
        prom_mean_acc=pr.mean_acc.max() if len(pr) else np.nan,
        prom_acc_pctile=pr.acc_pctile.max() if len(pr) else np.nan,
        prom_rho_stroma=pr.loc[pr.mean_acc.idxmax(), "rho_ESTIMATE_STROMAL_noCOL"] if len(pr) else np.nan,
        prom_fdr_stroma=pr.loc[pr.mean_acc.idxmax(), "fdr_ESTIMATE_STROMAL_noCOL"] if len(pr) else np.nan,
        prom_rho_pctile=pr.loc[pr.mean_acc.idxmax(), "rho_pctile_ESTIMATE_STROMAL_noCOL"] if len(pr) else np.nan,
        locus_median_rho=sub["rho_ESTIMATE_STROMAL_noCOL"].median(),
        locus_n_fdr05_pos=int(((sub.fdr_ESTIMATE_STROMAL_noCOL < 0.05) & (sub.rho_ESTIMATE_STROMAL_noCOL > 0)).sum()),
        locus_n_fdr05_neg=int(((sub.fdr_ESTIMATE_STROMAL_noCOL < 0.05) & (sub.rho_ESTIMATE_STROMAL_noCOL < 0)).sum())))
Su = pd.DataFrame(summ)
Su.to_csv(f"{R}/atac_locus_summary.csv", index=False)
pd.set_option("display.width", 250)
print(Su.to_string(index=False))

# ---- ATAC vs matched RNA of the same gene ----------------------------------
rna = pd.read_csv(f"{C}/rna_matched_genes.csv")
rna["patient"] = rna["sample"].str.slice(0, 12)
rna = rna.drop_duplicates("patient").set_index("patient").reindex(upat)
rows = []
for g in ["COL1A1", "COL3A1", "EPCAM", "ESR1", "DCN", "LUM", "GAPDH"]:
    if g not in rna.columns: continue
    sub = L[(L.gene == g) & L.is_promoter]
    if not len(sub): 
        sub = L[(L.gene == g)].nlargest(1, "mean_acc")
        tag = "best_locus_peak"
    else:
        sub = sub.nlargest(1, "mean_acc"); tag = "promoter_peak"
    i = res.index[res.name == sub.name.values[0]][0]
    y = rna[g].values[ok]; x = P[i, ok]
    r, p = stats.spearmanr(x, y)
    rows.append(dict(gene=g, peak=sub.name.values[0], peak_type=tag, n=int(ok.sum()),
                     rho_atac_vs_rna=r, p=p))
AR = pd.DataFrame(rows); AR["fdr"] = bh(AR.p.values)
AR.to_csv(f"{R}/atac_vs_rna_same_patients.csv", index=False)
print("\n== promoter accessibility vs matched RNA (same patients) ==")
print(AR.to_string(index=False))
print("DONE")
