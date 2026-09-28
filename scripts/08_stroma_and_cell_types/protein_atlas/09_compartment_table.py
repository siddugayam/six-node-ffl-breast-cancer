#!/usr/bin/env python3
"""Part C: the compartment test at protein level. One row per gene, combining
HPA IHC annotator summaries, normal-breast/soft-tissue cell-type staining,
breast-cancer tumour-cell staining and its recorded subcellular location, and
the HPA single-cell fibroblast-vs-epithelial ratios (pan-tissue and breast)."""
import csv, os, collections
R="/path/to/revision/results/v6"
def rd(f): return list(csv.DictReader(open(os.path.join(R,f))))
summ  = {r["gene"]:r for r in rd("hpa_gene_summaries_v25.csv")}
norm  = {r["gene"]:r for r in rd("hpa_normal_breast_by_celltype.csv")}
bc    = {r["gene"]:r for r in rd("hpa_breast_cancer_ihc_nodes.csv")}
scpan = {r["gene"]:r for r in rd("hpa_compartment_singlecell.csv")}
scbr  = {r["gene"]:r for r in rd("hpa_breast_singlecell_summary_v23.csv")}
pat   = [r for r in rd("hpa_cancer_ihc_patients_v25.csv") if r["cancer"]=="Breast cancer"]
loc=collections.defaultdict(collections.Counter)
for r in pat: loc[r["gene"]][r["location"] or "-"]+=1

STROMAL=["COL1A1","COL3A1","FN1","PDGFRB","POSTN"]
TFS=["NFKB1","RELA","SP1","ETS1"]
rows=[]
for g in STROMAL+TFS:
    rows.append(dict(
      gene=g, arm=("collagen / ECM arm (predicted stromal)" if g in STROMAL else "TF arm (predicted not stromal)"),
      hpa_ihc_reliability=summ[g]["ihc_reliability"],
      hpa_normal_summary=summ[g]["normal_tissue_summary"],
      hpa_cancer_summary=summ[g]["cancer_summary"],
      normal_breast_glandular=norm[g]["breast_glandular"],
      normal_breast_myoepithelial=norm[g]["breast_myoepithelial"],
      normal_breast_adipocytes=norm[g]["breast_adipocytes"],
      soft_tissue_fibroblasts=f'{norm[g]["soft_tissue1_fibroblasts"]}/{norm[g]["soft_tissue2_fibroblasts"]}',
      breast_cancer_tumour_cells_H_M_L_ND=f'{bc[g]["high"]}/{bc[g]["medium"]}/{bc[g]["low"]}/{bc[g]["not_detected"]}',
      breast_cancer_n_patients=bc[g]["n_patients"],
      breast_cancer_recorded_location="; ".join(f"{k}={v}" for k,v in loc[g].most_common()),
      sc_pan_fibroblast_nCPM=scpan[g]["fibroblasts_nCPM"],
      sc_pan_max_breast_epithelial_nCPM=scpan[g]["max_breast_epithelial_nCPM"],
      sc_pan_fibroblast_over_epithelial=scpan[g]["fibroblast_over_max_breast_epithelial"],
      sc_breast_fibroblast_nTPM=scbr[g]["max_fibroblast_nTPM"],
      sc_breast_epithelial_nTPM=scbr[g]["max_breast_epithelial_nTPM"],
      sc_breast_fibroblast_over_epithelial=scbr[g]["fibroblast_over_epithelial"],
      sc_specificity=scpan[g]["sc_specificity_category"],
      sc_enhanced_cell_types=scpan[g]["sc_enhanced_cell_types"]))
with open(os.path.join(R,"hpa_compartment_test.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("hpa_compartment_test.csv",len(rows))
for r in rows:
    print(f"{r['gene']:7s} {r['arm'][:28]:30s} fibr/epi pan={r['sc_pan_fibroblast_over_epithelial']:>6s} breast={r['sc_breast_fibroblast_over_epithelial']:>6s}  tumourcells H/M/L/ND={r['breast_cancer_tumour_cells_H_M_L_ND']:>12s}")
