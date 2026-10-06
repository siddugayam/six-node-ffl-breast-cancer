# Encode accessibility

## Scripts in this folder

| Script | What it does |
|---|---|
| `atac_01_prepare.py` | Collapse TCGA-BRCA ATAC-seq (Corces 2018) technical replicates to sample level. |
| `atac_02_tf_target_compartment.R` | Which TF->target correlations in TCGA-BRCA survive adjustment for stromal content? |
| `atac_03_motif_scan.py` | PWM scan of all TCGA-BRCA ATAC peaks (Corces 2018 BRCA peak set) with JASPAR 2024 CORE vertebrate PFMs. |
| `atac_04_encode_fetch.py` | Fetch ENCODE hg38 peak calls for breast carcinoma lines, breast epithelium and fibroblasts, to ask whether the COL1A1/COL3A1 promoters are accessible in carcinoma cells vs fibroblasts. |
| `atac_05_accessibility.py` | (a) Are COL1A1/COL3A1 regulatory elements accessible in TCGA-BRCA tumours, and does that accessibility track stromal content? |
| `atac_06_encode_compartment.py` | (a, supporting) ENCODE chromatin accessibility in breast carcinoma cells vs fibroblasts at the COL1A1/COL3A1 promoters, and the extent to which the TCGA-BRCA tumour ATAC peak set captures fibroblast-specific elements. |
| `atac_06b_capture_matched.py` | Reproducibility-matched version of the compartment-capture comparison: elements called in ALL 5 samples of one compartment and in NONE of the other. |
| `atac_07_remap.py` | (b, supporting) Measured TF occupancy (ReMap 2022, hg38, all cell types) at the COL1A1 and COL3A1 promoters, placed against the local background of all RefSeq promoters inside the same cached ReMap intervals. |
| `atac_08_remap_binding.py` | (b)/(c) Measured genome-wide TF occupancy: ReMap 2022 non-redundant hg38 peaks for ETS1, NFKB1, RELA, SP1 at promoters (TSS +/-1 kb). |
| `atac_09_motif_enrichment.py` | (b)/(c) Do ETS1/NFKB1/RELA/SP1 motifs occur in accessible chromatin at COL1A1/COL3A1, and at targets whose correlation survives compartment adjustment? |
| `atac_10_fibroblast_motifs.py` | (b) In the compartment where the collagens are actually transcribed: do ETS1/NFKB1/RELA/SP1 motifs occur in the fibroblast-accessible COL1A1 and COL3A1 promoters? |
| `atac_11_locus_stroma.py` | (a) Locus-level test: does aggregate ATAC accessibility at a gene locus track RNA-derived stromal content in the same 74 TCGA-BRCA patients? |
| `atac_12_remap_hot.py` | (b, supporting) Genome-wide HOT-promoter index from ReMap 2022 non-redundant hg38 (all TFs, all cell types): how many distinct TFs have a ChIP-seq peak in each RefSeq-Select promoter (TSS +/-1 kb)? |
| `atac_compartment_specificity_tests_rebuild.py` | Fisher exact tests of fibroblast against epithelial accessibility for each gene and window, recomputed from the accessibility matrix of atac_06 and compared with results/v6/atac_encode_compartment_specificity_tests.csv. |
