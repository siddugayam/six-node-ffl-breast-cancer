#!/usr/bin/env python3
"""
celllines_09_network_sign_tcga.py

The same network-level sign test that celllines_03_controls.py ran in 50 breast cancer
cell lines, run in 1,066 paired TCGA-BRCA tumours, so the two are directly comparable:
over ALL canonical targets of a miRNA, what fraction correlate negatively with it, and
what is the mean rho, against a permutation null built from the other measured miRNAs.
Also reports the CAF-low tertile version.
"""
import os, subprocess
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")
RNG = np.random.default_rng(20260909)
N_PERM = 2000

FOCUS = ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-miR-200c", "hsa-miR-21",
         "hsa-let-7b", "hsa-miR-101", "hsa-miR-34a", "hsa-miR-124"]

edges = pd.read_csv(os.path.join(DATA, "canonical_edges.tsv"), sep="\t")
mt = edges[edges.edge_type == "miRNA_target"]
need_genes = sorted(set(mt.loc[mt.source.isin(FOCUS), "target"]) |
                    {"DCN", "LUM", "FAP", "THY1"})

rcode = f'''
g <- readRDS("{DATA}/brca_gene_expr.rds"); m <- readRDS("{DATA}/brca_mirna_expr_canonical.rds")
p <- readRDS("{DATA}/brca_pheno.rds")
tum <- p$sample[p$sample_type == "Primary Tumor"]
both <- intersect(intersect(colnames(g), tum), intersect(colnames(m), tum))
cat("N_paired", length(both), "\\n")
gn <- readLines("{RES}/_need_genes.txt"); gn <- gn[gn %in% rownames(g)]
cat("N_genes_found", length(gn), "\\n")
write.table(g[gn, both], "{RES}/_tcga_targets.tsv", sep="\\t", quote=FALSE)
write.table(m[, both], "{RES}/_tcga_allmirs.tsv", sep="\\t", quote=FALSE)
'''
open(os.path.join(RES, "_need_genes.txt"), "w").write("\n".join(need_genes) + "\n")
rp = os.path.join(ROOT, "scripts", "08_stroma_and_cell_types", "cell_lines", "_celllines_tcga_targets.R")
open(rp, "w").write(rcode)
o = subprocess.run(["/usr/bin/Rscript", rp], capture_output=True, text=True)
print(o.stdout, o.stderr)
assert o.returncode == 0, o.stderr

G = pd.read_csv(os.path.join(RES, "_tcga_targets.tsv"), sep="\t", index_col=0).T
M = pd.read_csv(os.path.join(RES, "_tcga_allmirs.tsv"), sep="\t", index_col=0).T
G = G.loc[M.index]
print("TCGA target matrix", G.shape, "miRNA matrix", M.shape)

zc = lambda df: (df - df.mean()) / df.std(ddof=1)
caf = zc(G[["DCN", "LUM", "FAP", "THY1"]]).mean(axis=1)
q = np.quantile(caf.values, [1 / 3, 2 / 3])
tert = np.digitize(caf.values, q)


def rank_corr_matrix(x, Y):
    """Spearman of vector x against every column of Y, vectorised."""
    rx = stats.rankdata(x)
    rx = (rx - rx.mean()) / rx.std()
    RY = np.apply_along_axis(stats.rankdata, 0, Y)
    RY = (RY - RY.mean(0)) / RY.std(0)
    return (rx @ RY) / len(rx)


rows = []
allm = [c for c in M.columns]
for subset, mask in [("all_paired_tumours", np.ones(len(G), bool)),
                     ("CAF_low_tertile", tert == 0)]:
    Gs, Ms = G.values[mask], M.values[mask]
    Gcols, Mcols = list(G.columns), list(M.columns)
    for mi in FOCUS:
        if mi not in Mcols:
            continue
        tg = sorted(set(mt.loc[mt.source == mi, "target"]) & set(Gcols))
        if len(tg) < 5:
            continue
        idx = [Gcols.index(t) for t in tg]
        Y = Gs[:, idx]
        keep = Y.std(0) > 0
        Y = Y[:, keep]
        x = Ms[:, Mcols.index(mi)]
        r = rank_corr_matrix(x, Y)
        frac = float(np.mean(r < 0)); mr = float(np.mean(r))
        others = [c for c in Mcols if c != mi]
        nf = np.empty(N_PERM); nm = np.empty(N_PERM)
        for i in range(N_PERM):
            xo = Ms[:, Mcols.index(others[RNG.integers(len(others))])]
            if xo.std() == 0:
                nf[i], nm[i] = np.nan, np.nan
                continue
            rr = rank_corr_matrix(xo, Y)
            nf[i] = np.mean(rr < 0); nm[i] = np.mean(rr)
        p_f = (np.nansum(nf >= frac) + 1) / (np.sum(np.isfinite(nf)) + 1)
        p_m = (np.nansum(nm <= mr) + 1) / (np.sum(np.isfinite(nm)) + 1)
        rows.append(dict(dataset=f"TCGA-BRCA {subset} (n={int(mask.sum())})", miRNA=mi,
                         n_targets=int(keep.sum()), frac_negative=frac, mean_rho=mr,
                         null_frac_negative=float(np.nanmean(nf)),
                         null_mean_rho=float(np.nanmean(nm)),
                         p_frac_neg=float(p_f), p_mean_rho=float(p_m), n_perm=N_PERM))
        print(rows[-1], flush=True)

t = pd.DataFrame(rows)
t.to_csv(os.path.join(RES, "celllines_F_network_sign_tcga.csv"), index=False)

cl = pd.read_csv(os.path.join(RES, "celllines_B_network_sign_test.csv"))
comb = pd.concat([cl, t], ignore_index=True)
comb.to_csv(os.path.join(RES, "celllines_F_network_sign_combined.csv"), index=False)
print(comb.round(4).to_string(index=False))
print("DONE")
