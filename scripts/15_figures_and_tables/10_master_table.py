#!/usr/bin/env python3
"""One master table: the 20 protein-coding prioritised nodes (10 TFs + 10 genes),
HPA v25.1 - breast-cancer IHC, normal-breast IHC by cell type, pathology-atlas
prognostics, subcellular location and antibody reliability. The 6 comparator
genes used for the compartment test are appended and labelled."""
import csv, os
R="/path/to/revision/results/v6"
def rd(f): return {r["gene"]:r for r in csv.DictReader(open(os.path.join(R,f)))}
b,n,p,fl,sc = (rd("hpa_breast_cancer_ihc_nodes.csv"), rd("hpa_normal_breast_by_celltype.csv"),
               rd("hpa_prognostic_breast.csv"), rd("hpa_antibody_reliability_flags.csv"),
               rd("hpa_breast_singlecell_summary_v23.csv"))
order=["E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3",
       "CCND2","COL1A1","STAT5A","MYBL2","FN1","PDGFRB","MET","CXCL12","MMP14","PLAU",
       "COL3A1","POSTN","NFKB1","RELA","SP1","ETS1"]
rows=[]
for g in order:
    rows.append(dict(gene=g, role=b[g]["role"], hpa_version="25.1",
      antibody=b[g]["selected_antibody"], antibody_matches_HPA_table=b[g]["selected_matches_HPA_pathology_table"],
      bc_high=b[g]["high"], bc_medium=b[g]["medium"], bc_low=b[g]["low"], bc_not_detected=b[g]["not_detected"],
      bc_n_patients=b[g]["n_patients"], bc_pct_high_or_medium=b[g]["pct_high_or_medium"],
      bc_n_antibodies=b[g]["n_antibodies_with_breast_IHC"],
      normal_breast_glandular=n[g]["breast_glandular"], normal_breast_myoepithelial=n[g]["breast_myoepithelial"],
      normal_breast_adipocytes=n[g]["breast_adipocytes"],
      normal_soft_tissue_fibroblasts=f'{n[g]["soft_tissue1_fibroblasts"]}/{n[g]["soft_tissue2_fibroblasts"]}',
      hpa_prognostic_call=p[g]["hpa_TCGA_call"], hpa_prognostic_direction=p[g]["hpa_TCGA_direction"] or "-",
      hpa_prognostic_p=p[g]["hpa_TCGA_p"],
      hpa_validation_call=p[g]["hpa_validation_call"], hpa_validation_p=p[g]["hpa_validation_p"],
      ihc_reliability=fl[g]["ihc_reliability"], if_reliability=fl[g]["if_reliability"],
      subcellular_main=fl[g]["subcellular_main"], subcellular_additional=fl[g]["subcellular_additional"],
      FLAG_gene_level_uncertain=fl[g]["FLAG_gene_level_uncertain"], FLAG_no_ihc=fl[g]["FLAG_no_ihc"],
      sc_breast_fibroblast_nTPM=sc[g]["max_fibroblast_nTPM"], sc_breast_epithelial_nTPM=sc[g]["max_breast_epithelial_nTPM"],
      sc_breast_fibroblast_over_epithelial=sc[g]["fibroblast_over_epithelial"],
      hpa_normal_summary=n[g]["hpa_normal_tissue_summary"], hpa_cancer_summary=b[g]["hpa_cancer_summary"]))
with open(os.path.join(R,"hpa_master_nodes.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("hpa_master_nodes.csv",len(rows))
