#!/usr/bin/env python3
"""HPA breast-tissue-specific single-cell clusters (rna_single_cell_type_tissue.tsv, v23):
fibroblast vs breast-epithelial nTPM for the 26 genes."""
import csv, os, collections
D="/path/to/revision/data/hpa"
R="/path/to/revision/results/v6"
TFS=["E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3"]
GEN=["CCND2","COL1A1","STAT5A","MYBL2","FN1","PDGFRB","MET","CXCL12","MMP14","PLAU"]
EX =["COL3A1","POSTN","NFKB1","RELA","SP1","ETS1"]
ROLE={g:"TF (prioritised)" for g in TFS}; ROLE.update({g:"gene (prioritised)" for g in GEN})
ROLE.update({g:"comparator (compartment test)" for g in EX})
ALL=TFS+GEN+EX
keep=collections.defaultdict(list)
with open(os.path.join(D,"rna_single_cell_type_tissue.tsv"),newline="",encoding="utf-8") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["Tissue"]=="breast" and r["Gene name"] in ALL:
            keep[r["Gene name"]].append(r)
long_rows=[]
for g in ALL:
    for r in keep.get(g,[]):
        long_rows.append(dict(gene=g,role=ROLE[g],tissue=r["Tissue"],cluster=r["Cluster"],
                              cell_type=r["Cell type"],read_count=r["Read count"],nTPM=r["nTPM"]))
with open(os.path.join(R,"hpa_breast_singlecell_clusters_v23.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(long_rows[0].keys())); w.writeheader(); w.writerows(long_rows)
print("hpa_breast_singlecell_clusters_v23.csv",len(long_rows))

summ=[]
for g in ALL:
    rs=keep.get(g,[])
    fib=[float(r["nTPM"]) for r in rs if r["Cell type"]=="fibroblasts"]
    gla=[float(r["nTPM"]) for r in rs if r["Cell type"] in
         ("breast glandular cells","breast myoepithelial cells")]
    adi=[float(r["nTPM"]) for r in rs if r["Cell type"]=="adipocytes"]
    mf=max(fib) if fib else None; mg=max(gla) if gla else None; ma=max(adi) if adi else None
    summ.append(dict(gene=g,role=ROLE[g],
        n_breast_clusters=len(rs),
        max_fibroblast_nTPM=mf, max_breast_epithelial_nTPM=mg, max_adipocyte_nTPM=ma,
        fibroblast_over_epithelial=(round(mf/mg,1) if mf and mg and mg>0 else "")))
with open(os.path.join(R,"hpa_breast_singlecell_summary_v23.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
print("hpa_breast_singlecell_summary_v23.csv",len(summ))
for s in summ:
    print(f"{s['gene']:8s} fib={str(s['max_fibroblast_nTPM']):>8s} epi={str(s['max_breast_epithelial_nTPM']):>8s} adip={str(s['max_adipocyte_nTPM']):>8s} ratio={str(s['fibroblast_over_epithelial']):>7s}")
