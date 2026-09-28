#!/usr/bin/env python3
"""Human Protein Atlas single-cell RNA: cell-type source of network genes.
Inputs (already downloaded, cache/celltype/):
  rna_single_cell_type_tissue.tsv  (HPA v23, per tissue x cluster x cell type, nTPM)
  rna_single_cell_type.tsv         (HPA v25, pooled per cell type, nCPM)
Outputs: cache/celltype/hpa_breast_celltype.csv, hpa_pantissue_celltype.csv
"""
import csv, sys, collections

CACHE = "/path/to/revision/cache/celltype/"
genes = set(l.strip() for l in open(CACHE + "network_protein_nodes.txt") if l.strip())
# markers used to sanity-check the cell-type assignment
markers = {"PTPRC","EPCAM","KRT18","KRT19","CDH1","PECAM1","VWF","DCN","LUM","FAP",
           "PDGFRA","PDGFRB","THY1","ACTA2","CD68","CD3E","MS4A1","COL1A2","COL5A1",
           "COL6A3","POSTN","SPARC","FN1","MKI67","ESR1","KRT14","TP63","MYLK"}
want = genes | markers
print("genes wanted:", len(want), flush=True)

# ---------- breast, tissue-resolved ----------
acc = collections.defaultdict(lambda: [0.0, 0.0, 0])  # (gene, celltype) -> [sum nTPM, sum reads, n clusters]
n = 0
with open(CACHE + "rna_single_cell_type_tissue.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t")
    hdr = next(rd)
    assert hdr[:7] == ['Gene','Gene name','Tissue','Cluster','Cell type','Read count','nTPM'], hdr
    for row in rd:
        n += 1
        if row[2] != "breast":
            continue
        g = row[1]
        if g not in want:
            continue
        k = (g, row[4])
        a = acc[k]
        a[0] += float(row[6]); a[1] += float(row[5]); a[2] += 1
print("rows scanned:", n, " breast gene-cluster rows kept:", sum(v[2] for v in acc.values()), flush=True)

with open(CACHE + "hpa_breast_celltype.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene","cell_type","mean_nTPM","n_clusters","total_read_count"])
    for (g, ct), v in sorted(acc.items()):
        w.writerow([g, ct, round(v[0]/v[2], 4), v[2], int(v[1])])

# ---------- pan-tissue pooled ----------
acc2 = collections.defaultdict(float)
n = 0
with open(CACHE + "rna_single_cell_type.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t")
    hdr = next(rd)
    assert hdr[:4] == ['Gene','Gene name','Cell type','nCPM'], hdr
    for row in rd:
        n += 1
        if row[1] in want:
            acc2[(row[1], row[2])] = float(row[3])
print("pan-tissue rows scanned:", n, " kept:", len(acc2), flush=True)
with open(CACHE + "hpa_pantissue_celltype.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene","cell_type","nCPM"])
    for (g, ct), v in sorted(acc2.items()):
        w.writerow([g, ct, v])
print("done")
