#!/usr/bin/env python3
"""Per-cell-type expression of the network's protein-coding nodes in 26 primary breast
tumours (GSE176078, Wu et al. 2021, 100,064 cells), combined with the Human Protein Atlas
normal-breast single-cell reference.  Writes results/multiomics/celltype_source.csv."""
import csv, collections, math, statistics

CACHE = "/path/to/revision/cache/celltype/"
RES   = "/path/to/revision/results/multiomics/"

# ---- cell metadata ----
cell_ct, cell_pat = {}, {}
with open(CACHE + "cell_index_celltype.tsv") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        i = int(r["cell_index"]); cell_ct[i] = r["celltype_major"]; cell_pat[i] = r["patient"]
n_cells_ct = collections.Counter(cell_ct.values())
print("cells:", len(cell_ct), "| per type:", dict(n_cells_ct))

# ---- per-cell library sizes ----
cell_tot = {}
with open(CACHE + "scrna_cell_totals.tsv") as fh:
    for line in fh:
        a, b = line.split()
        cell_tot[int(a)] = int(b)
umi_ct = collections.Counter()
for i, t in cell_tot.items():
    umi_ct[cell_ct[i]] += t
print("total UMI per cell type:", dict(umi_ct))

# ---- gene counts ----
gsum   = collections.defaultdict(int)      # (gene, ct) -> summed UMI
gcells = collections.defaultdict(int)      # (gene, ct) -> n cells expressing
gcp10k = collections.defaultdict(float)    # (gene, ct) -> sum of log1p(CP10K)
n = 0
with open(CACHE + "scrna_gene_counts.tsv") as fh:
    for line in fh:
        g, c, v = line.rstrip("\n").split("\t")
        c = int(c); v = int(v); n += 1
        ct = cell_ct[c]
        gsum[(g, ct)] += v
        gcells[(g, ct)] += 1
        gcp10k[(g, ct)] += math.log1p(v / cell_tot[c] * 1e4)
print("gene-cell entries read:", n)

genes = sorted({k[0] for k in gsum})
cts = sorted(n_cells_ct)
rows = []
for g in genes:
    tot_g = sum(gsum[(g, ct)] for ct in cts)
    for ct in cts:
        s = gsum[(g, ct)]
        rows.append(dict(
            gene=g, cell_type=ct,
            pseudobulk_CPM=round(s / umi_ct[ct] * 1e6, 3),
            pct_cells_expressing=round(100.0 * gcells[(g, ct)] / n_cells_ct[ct], 2),
            mean_log1p_CP10K=round(gcp10k[(g, ct)] / n_cells_ct[ct], 4),
            frac_of_gene_total_UMI=round(s / tot_g, 4) if tot_g else None,
            n_cells=n_cells_ct[ct], total_UMI_in_celltype=umi_ct[ct]))
with open(CACHE + "scrna_celltype_expression.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print("scrna per-celltype rows:", len(rows), "genes:", len(genes))

# ---- summarise per gene: dominant cell type, CAF vs cancer-epithelial ratio ----
byg = collections.defaultdict(dict)
for r in rows:
    byg[r["gene"]][r["cell_type"]] = r
hpa = collections.defaultdict(dict)
with open(CACHE + "hpa_breast_celltype.csv") as fh:
    for r in csv.DictReader(fh):
        hpa[r["gene"]][r["cell_type"]] = float(r["mean_nTPM"])

out = []
for g in genes:
    d = byg[g]
    tot = sum(d[ct]["pseudobulk_CPM"] for ct in cts)
    top = max(cts, key=lambda c: d[c]["pseudobulk_CPM"])
    caf = d["CAFs"]["pseudobulk_CPM"]; can = d["Cancer Epithelial"]["pseudobulk_CPM"]
    h = hpa.get(g, {})
    hfib = h.get("fibroblasts"); hgl = h.get("breast glandular cells")
    htop = max(h, key=lambda c: h[c]) if h else None
    out.append(dict(
        gene=g,
        scrna_top_celltype=top,
        scrna_top_CPM=d[top]["pseudobulk_CPM"],
        scrna_CAF_CPM=caf, scrna_CancerEpi_CPM=can,
        scrna_CAF_over_CancerEpi=round(caf / can, 3) if can > 0 else None,
        scrna_log2_CAF_over_CancerEpi=round(math.log2((caf + 0.1) / (can + 0.1)), 3),
        scrna_frac_UMI_from_CAF=d["CAFs"]["frac_of_gene_total_UMI"],
        scrna_frac_UMI_from_CancerEpi=d["Cancer Epithelial"]["frac_of_gene_total_UMI"],
        scrna_pct_CAFs_expressing=d["CAFs"]["pct_cells_expressing"],
        scrna_pct_CancerEpi_expressing=d["Cancer Epithelial"]["pct_cells_expressing"],
        scrna_pct_Myeloid_expressing=d["Myeloid"]["pct_cells_expressing"],
        scrna_pct_Tcells_expressing=d["T-cells"]["pct_cells_expressing"],
        scrna_CPM_Endothelial=d["Endothelial"]["pseudobulk_CPM"],
        scrna_CPM_Myeloid=d["Myeloid"]["pseudobulk_CPM"],
        scrna_CPM_Tcells=d["T-cells"]["pseudobulk_CPM"],
        scrna_CPM_PVL=d["PVL"]["pseudobulk_CPM"],
        scrna_CPM_NormalEpi=d["Normal Epithelial"]["pseudobulk_CPM"],
        scrna_CPM_Bcells=d["B-cells"]["pseudobulk_CPM"],
        scrna_CPM_Plasmablasts=d["Plasmablasts"]["pseudobulk_CPM"],
        hpa_breast_top_celltype=htop,
        hpa_breast_top_nTPM=h.get(htop) if htop else None,
        hpa_breast_fibroblasts_nTPM=hfib,
        hpa_breast_glandular_nTPM=hgl,
        hpa_fibroblast_over_glandular=round(hfib / hgl, 2) if (hfib is not None and hgl) else None,
    ))
with open(RES + "celltype_source.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)
print("celltype_source.csv rows:", len(out))

focus = ["COL1A1","COL3A1","NFKB1","RELA","SP1","ETS1","VEGFA","EZH2",
         "EPCAM","PTPRC","DCN","FAP","MYC","TP53","E2F1","CCND2","STAT3","TGFBR2","HIF1A"]
print("\n%-8s %-18s %9s %9s %8s %8s %8s  %-16s %8s" % (
    "gene","scRNA top","CAF_CPM","Epi_CPM","CAF/Epi","%CAF+","%Epi+","HPA breast top","fib/gland"))
for g in focus:
    r = next((x for x in out if x["gene"] == g), None)
    if r is None: print(g, "NOT FOUND"); continue
    print("%-8s %-18s %9.1f %9.1f %8s %8.1f %8.1f  %-16s %8s" % (
        g, r["scrna_top_celltype"], r["scrna_CAF_CPM"], r["scrna_CancerEpi_CPM"],
        r["scrna_CAF_over_CancerEpi"], r["scrna_pct_CAFs_expressing"],
        r["scrna_pct_CancerEpi_expressing"], str(r["hpa_breast_top_celltype"]),
        r["hpa_fibroblast_over_glandular"]))
