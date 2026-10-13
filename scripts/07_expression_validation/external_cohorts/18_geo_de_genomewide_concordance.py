#!/usr/bin/env python3
"""
Genome-wide direction concordance of the TCGA-BRCA differential expression in the three GEO tumour-vs-normal series.

For each series (GSE42568, GSE45827, GSE10780) the genes measured in both TCGA-BRCA and the series are taken from
results/external_GEO_DE_vs_TCGA.csv (written by 16_external_geo_de.R), and the genes that are differentially expressed in
TCGA-BRCA (Benjamini-Hochberg q < 0.05 and |log2 fold change| > 1, the criterion of Section 2.4) are counted together with
the number whose log2 fold change has the same sign in the series. The looser criterion q < 0.05 alone is also tabulated, and each count is given for all genes and for the protein-coding
network nodes alone.

Input : results/external_GEO_DE_vs_TCGA.csv
Output: results/external_GEO_DE_genomewide_concordance.csv
"""
import csv, os
from collections import defaultdict

BASE = "/path/to/six-node-ffl-breast-cancer"
rows = list(csv.DictReader(open(os.path.join(BASE, "results", "external_GEO_DE_vs_TCGA.csv"), encoding="utf-8")))

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

by = defaultdict(list)
for r in rows:
    by[r["gse"]].append(r)

out = []
SUBSETS = (("all genes", lambda r: True), ("network nodes", lambda r: r["is_network_node"] == "TRUE"))
CRITERIA = (("q<0.05 and |log2FC|>1", 1.0), ("q<0.05", 0.0))
for gse in ("GSE42568", "GSE45827", "GSE10780"):
    for subset, keep in SUBSETS:
        rs = [r for r in by[gse] if keep(r)]
        for label, lfc_min in CRITERIA:
            sel = [r for r in rs
                   if num(r["q_tcga"]) is not None and num(r["q_tcga"]) < 0.05
                   and num(r["logFC_tcga"]) is not None and abs(num(r["logFC_tcga"])) > lfc_min
                   and num(r["logFC_ext"]) is not None]
            same = sum(1 for r in sel if num(r["logFC_tcga"]) * num(r["logFC_ext"]) > 0)
            out.append({"series": gse, "gene_set": subset, "tcga_criterion": label, "genes_measured_in_both": len(rs),
                        "tcga_significant_and_measured": len(sel), "same_direction": same,
                        "pct_same_direction": round(100 * same / len(sel), 1)})

path = os.path.join(BASE, "results", "external_GEO_DE_genomewide_concordance.csv")
with open(path, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)
for o in out:
    print(o)
