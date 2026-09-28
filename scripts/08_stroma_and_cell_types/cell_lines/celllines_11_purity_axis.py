#!/usr/bin/env python3
"""
celllines_11_purity_axis.py

Decisive test of whether the bulk miR-29 -> collagen correlation is anything more than a
tumour-composition (purity) effect.

For EVERY miRNA measured in TCGA-BRCA (n = 496), compute
    a = Spearman(miRNA, CAF score)          (how much it tracks stromal content)
    b = Spearman(miRNA, COL1A1)  and COL3A1 (the "regulatory" correlation the paper uses)
If b is simply a linear function of a, then a negative miRNA->collagen correlation carries
no regulatory information at all.  We then ask whether miR-29a/b/c are outliers from that
relationship (residual + its rank among all 496 miRNAs), and do the same for all 20,530
genes as putative TF regulators.
"""
import os, subprocess
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")

rcode = f'''
g <- readRDS("{DATA}/brca_gene_expr.rds"); m <- readRDS("{DATA}/brca_mirna_expr_canonical.rds")
p <- readRDS("{DATA}/brca_pheno.rds")
tum <- p$sample[p$sample_type == "Primary Tumor"]
both <- intersect(intersect(colnames(g), tum), intersect(colnames(m), tum))
cat("N_paired", length(both), "\\n")
saveRDS(list(g = g[, both], m = m[, both]), "{RES}/_tcga_full_paired.rds")
'''
rp = os.path.join(ROOT, "scripts", "08_stroma_and_cell_types", "cell_lines", "_celllines_full_pull.R")
open(rp, "w").write(rcode)
o = subprocess.run(["/usr/bin/Rscript", rp], capture_output=True, text=True)
print(o.stdout, o.stderr)
assert o.returncode == 0

rcode2 = f'''
x <- readRDS("{RES}/_tcga_full_paired.rds")
write.table(x$g, gzfile("{RES}/_tcga_g.tsv.gz"), sep="\\t", quote=FALSE)
write.table(x$m, "{RES}/_tcga_m.tsv", sep="\\t", quote=FALSE)
'''
rp2 = os.path.join(ROOT, "scripts", "08_stroma_and_cell_types", "cell_lines", "_celllines_full_write.R")
open(rp2, "w").write(rcode2)
o = subprocess.run(["/usr/bin/Rscript", rp2], capture_output=True, text=True)
assert o.returncode == 0, o.stderr

Gg = pd.read_csv(os.path.join(RES, "_tcga_g.tsv.gz"), sep="\t", index_col=0).T
Mm = pd.read_csv(os.path.join(RES, "_tcga_m.tsv"), sep="\t", index_col=0).T
Gg = Gg.loc[Mm.index]
print("genes", Gg.shape, "miRNAs", Mm.shape)
zc = lambda df: (df - df.mean()) / df.std(ddof=1)
caf = zc(Gg[["DCN", "LUM", "FAP", "THY1"]]).mean(axis=1).values


def rc(x, Y):
    rx = stats.rankdata(x); rx = (rx - rx.mean()) / rx.std()
    RY = np.apply_along_axis(stats.rankdata, 0, Y)
    sd = RY.std(0); sd[sd == 0] = np.nan
    RY = (RY - RY.mean(0)) / sd
    return (rx @ RY) / len(rx)


def build(X, names, label):
    a = rc(caf, X.values)                       # here a_i = corr(CAF, feature_i)
    # need corr(feature_i, COL1A1): same thing, symmetric
    b1 = rc(Gg["COL1A1"].values, X.values)
    b3 = rc(Gg["COL3A1"].values, X.values)
    d = pd.DataFrame(dict(feature=names, rho_vs_CAFscore=a,
                          rho_vs_COL1A1=b1, rho_vs_COL3A1=b3))
    d = d.dropna()
    for c in ["COL1A1", "COL3A1"]:
        y = d[f"rho_vs_{c}"].values
        x = d["rho_vs_CAFscore"].values
        sl, ic, r, p, se = stats.linregress(x, y)
        d[f"resid_{c}"] = y - (sl * x + ic)
        print(f"[{label}] rho_vs_{c} ~ rho_vs_CAFscore : slope {sl:.3f}, "
              f"R = {r:.4f}, R2 = {r**2:.4f}, n = {len(d)}")
        d[f"resid_{c}_pctile"] = stats.rankdata(d[f"resid_{c}"]) / len(d)
    d["panel"] = label
    return d


mi = build(Mm, list(Mm.columns), "miRNA (n=%d)" % Mm.shape[1])
ge = build(Gg, list(Gg.columns), "gene (n=%d)" % Gg.shape[1])
out = pd.concat([mi, ge], ignore_index=True)
out.to_csv(os.path.join(RES, "celllines_F_purity_axis.csv"), index=False)

show = ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-miR-200c", "hsa-miR-21",
        "hsa-miR-101", "hsa-let-7b", "ETS1", "NFKB1", "RELA", "SP1", "MKL1", "STAT6",
        "TFAP2A", "MYB", "COL1A2", "FN1"]
print("\n", out[out.feature.isin(show)][
    ["panel", "feature", "rho_vs_CAFscore", "rho_vs_COL1A1", "resid_COL1A1",
     "resid_COL1A1_pctile", "rho_vs_COL3A1", "resid_COL3A1", "resid_COL3A1_pctile"]]
    .round(4).to_string(index=False))

# how extreme is miR-29 among all miRNAs for raw correlation vs for residual?
for f in ["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"]:
    row = mi[mi.feature == f]
    if len(row) == 0:
        continue
    raw_rank = (mi.rho_vs_COL1A1 < row.rho_vs_COL1A1.iloc[0]).mean()
    res_rank = row.resid_COL1A1_pctile.iloc[0]
    print(f"{f}: raw rho vs COL1A1 = {row.rho_vs_COL1A1.iloc[0]:.3f} "
          f"(percentile {raw_rank:.3f} among {len(mi)} miRNAs); "
          f"residual percentile after removing the CAF axis = {res_rank:.3f}")
print("DONE")
