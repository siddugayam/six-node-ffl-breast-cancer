#!/usr/bin/env python3
import pandas as pd, json, os, hashlib
BASE="/path/to/revision"; RES=f"{BASE}/results/v6"
rows=[
 dict(arm="PDX", resource="NCI Patient-Derived Models Repository (PDMR)",
      accession="cBioPortal study pancan_pdmr_2025; profile pancan_pdmr_2025_rna_seq_v2_mrna",
      url="https://www.cbioportal.org/study/summary?id=pancan_pdmr_2025",
      access="public REST API, retrieved 2026-09-09",
      content="RSEM RNA-seq V2 mRNA, human genome; 6,272 samples pan-cancer, 300 breast-cancer samples used",
      detail="14 models with a patient ORIGINATOR specimen and >=1 PDX passage; 265 PDX samples; 21 PDOrg/PDC in-vitro samples"),
 dict(arm="PDX", resource="TCGA-BRCA bulk primary tumours (local Xena matrix)",
      accession="data/brca_gene_expr.rds (symbol-keyed 20,530 x 1,211)",
      url="local", access="on disk",
      content="log2 RSEM normalised counts", detail="1,097 Primary Tumor samples used as the external bulk reference"),
 dict(arm="PDX", resource="Wu et al. 2021 breast atlas pseudobulk (for CAF/carcinoma ratio)",
      accession="GSE176078 / local cache/celltype/scrna_pseudobulk_all.tsv",
      url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
      access="on disk", content="per-cell-type summed UMI, 9 cell types",
      detail="CAF vs Cancer Epithelial CPM ratio computed for all genes"),
 dict(arm="spatial", resource="10x Genomics public Visium human breast cancer",
      accession="V1_Breast_Cancer_Block_A_Section_1; V1_Breast_Cancer_Block_A_Section_2 (Space Ranger 1.1.0); Visium_FFPE_Human_Breast_Cancer (Space Ranger 1.3.0)",
      url="https://cf.10xgenomics.com/samples/spatial-exp/",
      access="public download, retrieved 2026-09-09",
      content="filtered_feature_bc_matrix.h5 + spatial/ (tissue image, positions, scalefactors)",
      detail="3,798 / 3,985 / 2,516 spots passing >=500 UMI"),
 dict(arm="spatial", resource="Wu et al. 2021 Nat Genet spatial transcriptomics (Visium)",
      accession="Zenodo record 4739739 (filtered_count_matrices, spatial, metadata)",
      url="https://zenodo.org/records/4739739",
      access="public download, retrieved 2026-09-09",
      content="6 Visium sections: 1142243F, 1160920F, CID4290, CID4465, CID44971, CID4535, with the authors' per-spot pathology annotation",
      detail="4,755 / 4,890 / 2,432 / 1,211 / 1,162 / 1,127 spots passing >=500 UMI"),
]
pd.DataFrame(rows).to_csv(f"{RES}/v6_data_provenance.csv", index=False)
print(pd.DataFrame(rows)[["arm","resource","accession"]].to_string(index=False))
