#!/usr/bin/env python3
import os,csv,json,datetime,hashlib
BASE="/path/to/revision"
OUT=BASE+"/results/v6"
rows=[]
def add(item, source, detail):
    rows.append(dict(item=item, source=source, detail=detail))
add("retrieval date","2026-09-09","all GEO / Bioconductor retrievals in this run")
add("GEO search engine","NCBI E-utilities esearch/esummary on db=gds",
    "5 systematic queries (breast x miRNA; breast x ncRNA DataSet Type; +outcome terms; "
    "+serum/plasma terms; +neoadjuvant terms), GSE entry type, Homo sapiens; "
    "882 unique GSE returned; 150 breast+miRNA series with n>=40 triaged; "
    "155 series probed by downloading the full SOFT-brief record and reading every "
    "deposited !Sample_characteristics_ch1 field")
add("ArrayExpress / BioStudies","not separately searched",
    "ArrayExpress has mirrored its array-based records into GEO/BioStudies; every "
    "candidate named in the task and every hit above came back through the GEO route. "
    "This is a gap in coverage and is reported as such rather than claimed.")
add("miRNA Cox cohorts (6)","TCGA-BRCA; GSE19783; GSE22216; GSE37405; GSE59829; GSE78870",
    "1,732 patients, 483 events in the primary-endpoint pool")
add("cohorts new in this run","GSE59829 (n=123, DMFS); GSE78870 (n=106, TTP)",
    "the other four were already assembled in results/v4 (miR-130a only)")
add("miRNA binary-outcome cohorts (6)","GSE57897; GSE103161; GSE97811; GSE28321; GSE40267; GSE26659",
    "794 patients, 306 events; no follow-up time deposited, so logistic not Cox")
add("liquid-biopsy cohorts (5)","GSE73002; GSE44281; GSE110651; GSE118782; GSE68373",
    "serum n=4,113 + 419 + 147; plasma n=40 + 64")
add("mRNA cohorts new in this run (7)",
    "MAINZ/GSE11121; TRANSBIG/GSE7390; VDX/GSE2034+GSE5327; UPP/GSE3494; UNT/GSE2990; "
    "NKI/van de Vijver; STK/GSE1456",
    "1,579 patients, 478 events; six from the Bioconductor breastCancer* experiment "
    "packages, GSE1456 parsed from its GEO series matrix (GPL96 arm)")
add("pan-cancer matrices","UCSC Xena PANCAN on disk in /path/to/home/Desktop/DD/R_GPR/ML/",
    "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz (20,531 genes) and "
    "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz (743 miRNAs); "
    "10,509 samples carry both assays; 9,812 primary tumours in 32 cohorts; 30 cohorts with n>=50")
add("CAF-A signature","results/v2/cafA_gene_list_v3.txt (125 symbols)",
    "121 measurable in the pan-cancer matrix; z-scored across all pan-cancer primary "
    "tumours so cohort means are comparable")
add("meta-analysis method","DerSimonian-Laird random effects, metafor::rma(method='DL')",
    "per-SD log hazard ratios; I2, tau2, Q and its p reported; fixed-effect and Egger "
    "regression reported alongside in mirna_meta_FE_vs_RE_and_egger.csv")
add("probe -> mature miRNA rule","dominant non-star arm",
    "all probes whose annotated name matches the target after stripping hsa-, arm "
    "(-3p/-5p), precursor index and platform suffixes; the passenger ('*') arm is used "
    "only if no guide-arm probe exists; among the remainder the highest-mean probe is "
    "taken; every probe considered is recorded in cohorts_mirna_probe_mapping.csv")
add("cross-check against results/v4","miR-130a per-cohort HRs reproduce exactly",
    "TCGA 0.8679, GSE19783 0.9122, GSE22216 1.2256, GSE37405 0.8688 - identical to "
    "results/v4/multicohort_B_mir130a_cox_all.csv")
with open(OUT+"/cohorts_provenance.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["item","source","detail"]); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote cohorts_provenance.csv", len(rows), "rows")
