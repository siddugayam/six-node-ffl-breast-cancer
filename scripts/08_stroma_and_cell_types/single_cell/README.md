# Single cell

## Scripts in this folder

| Script | What it does |
|---|---|
| `05_breast_singlecell.py` | HPA breast-tissue-specific single-cell clusters (rna_single_cell_type_tissue.tsv, v23): fibroblast vs breast-epithelial nTPM for the 26 genes. |
| `19_mir29_target_set.R` | 19) Build the miR-29 family target set and its TCGA-BRCA correlations, so that the miR-29 axis can be tested in the functional-genomic screens and single-cell analyses. |
| `21_scrna_reference_build.py` | Build a breast-tumour-specific deconvolution reference from Wu et al. 2021 (GSE176078). |
| `30_sc_atlas2_pseudobulk.py` | Second (and much larger) breast atlas: the Human Breast Cancer Single Cell Atlas global object (621,200 cells, 138 donors, 8 constituent studies including Pal et al. 2021 and Bassez et al. 2021). |
| `30a_wu_pseudobulk.sh` | Stream the Wu et al. 2021 (GSE176078) 29,733 x 100,064 count matrix once and accumulate summed UMI per gene for three groupings: celltype_minor (CAF subsets), patient\|celltype_major and patient\|celltype_minor. |
| `31_sc_caf_subtypes.py` | CAF subtypes in the Human Breast Cancer Single Cell Atlas stromal compartment. |
| `31b_wu_caf_subtypes.py` | Cross-check in the atlas the project already uses (Wu et al. 2021, GSE176078), which carries its own myCAF-like / MSC iCAF-like CAF labels. |
| `31c_caf_summary_table.py` | Summary: assign each CAF subtype to myCAF / iCAF / apCAF by marker score, and state quantitatively which subtype carries the collagen programme and which carries the TFs. |
| `32_sc_celltype_enrichment.py` | (confirmation) and 2H (in which cell types are the miR-29 target sets expressed?). |
| `32b_sc_setcontrols.py` | Controls. (1) Is the CAF-weighting of the miR-29 target set just the collagen genes themselves? |
| `32c_wu_compartment_check.py` | Independent replication in the Wu et al. 2021 atlas (the project's original atlas), using its own celltype_minor pseudobulk, to check that the compartment bias of the miR-29 target sets is not specific to the Human Breast Cancer … |
| `33_sc_donor_pseudobulk.py` | Donor x cell-type pseudobulk from the global atlas, so that the TF -> COL1A1 association can be tested WITHIN the CAF compartment alone. |
| `34_sc_tf_collagen_within_caf.R` | Does the TF -> COL1A1 association exist WITHIN the CAF compartment alone? |
| `34b_caf_composition_control.py` | Control: is the surviving within-CAF TF -> COL1A1 association just CAF-subtype composition (donors with more myCAF-like cells have more collagen AND more RUNX2)? |
| `35_sc_cellchat_liana.py` | Cell-cell communication on the Human Breast Cancer Single Cell Atlas, focused on the CAF -> cancer-cell axis, collagen/ECM and TGF-beta. |
| `36_sc_cellchat.R` | CellChat v2 on the Human Breast Cancer Single Cell Atlas (subsampled): CAF subtypes to cancer cells, with the ECM/collagen and TGF-β pathways reported. |
| `36a_cellchat_input.py` | Export a subsampled, log-normalised expression matrix + cell labels for CellChat. |
| `36b_sc_cellchat_small.R` | As 36_sc_cellchat.R, on the smaller subsample in cache/v7/cc_small. |
| `ct_02_scrna_extract.sh` | Single streaming pass over the 178M-entry GSE176078 sparse matrix. |
| `ct_02b_pseudobulk_all.sh` | Second streaming pass: per-(gene x celltype_major) UMI sums and n-cells-expressing, for ALL 29,733 genes. |
