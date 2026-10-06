# Xenografts

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_fetch_pdmr.py` | Fetch NCI PDMR (cBioPortal study pancan_pdmr_2025) RNA-seq V2 RSEM expression for breast-cancer samples, for a reference gene panel + 40 target genes. |
| `03_pdx_analysis.py` | PDX natural experiment. NCI Patient-Derived Models Repository (PDMR), cBioPortal study pancan_pdmr_2025, profile pancan_pdmr_2025_rna_seq_v2_mrna (RSEM, human genome). |
| `04_pdx_foldchange_and_caf.py` | (a) Matched PDX-vs-originator log2 fold change, corrected for the compositional renormalisation that loss of the stromal transcriptome causes (each pair's log2 ratios are centred on the median log2 ratio of the reference panel). |
| `08_pdx_passage_trend.py` | Dose-response: does the human collagen signal fall progressively with PDX passage, as residual human stroma is replaced by mouse? |
