# Deconvolution and mediation

## Scripts in this folder

| Script | What it does |
|---|---|
| `05_mediation.R` | Formal causal-mediation analysis -- how much of TF -> COL1A1 is carried by CAF (stromal) content? |
| `20_deconv_install.R` | Install every deconvolution package we intend to use. |
| `22_deconv_prep.R` | Build the matrices every deconvolution method needs. |
| `23_cibersort_impl.R` | Faithful open re-implementation of the CIBERSORT nu-support-vector-regression deconvolution of Newman et al. 2015 (Nat Methods 12:453-457). |
| `24_wu_signature_matrix.R` | Build a BREAST-TUMOUR-SPECIFIC signature matrix for 9 major cell types from the Wu et al. 2021 (GSE176078) atlas, following the CIBERSORT signature-matrix construction recipe: (i) per-cell-type differential expression across … |
| `25_run_deconvolution.R` | Run every deconvolution method that installed, on the 1,097 TCGA-BRCA primary tumours. |
| `26_instaprism_bayesprism.R` | BayesPrism-family deconvolution (InstaPrism: the exact BayesPrism model solved with a deterministic fixed-point iteration instead of Gibbs sampling) using the Wu et al. 2021 breast-tumour atlas as reference. |
| `27_assemble_scores.R` | Assemble every per-sample cell-type estimate into one table, add the signature-score methods the project already used (so the new methods are compared with the incumbent on identical samples), and write: … |
| `28_mediation_by_method.R` | THE DECISIVE RE-TEST. Repeat the causal-mediation analysis of the TF -> collagen and miR-29 -\| collagen edges once per deconvolution method, using THAT METHOD'S fibroblast/stromal estimate as the single mediator. |
| `29_correlations.R` | (B) pairwise Spearman correlation matrix between every fibroblast/stromal estimate (D) correlation of every cell-type fraction from every method with COL1A1, COL3A1 and the miR-29 family (and, for context, ETS1/NFKB1) |
| `30_simulation_benchmark.R` | Which fibroblast estimate is actually accurate? |
| `31_gene_lengths.py` | Union-exon length per gene symbol from GENCODE v36 basic, used to test whether the log2(norm_count) -> TPM approximation in 22_deconv_prep.R changes any fibroblast estimate. |
| `32_tpm_length_sensitivity.R` | Xena TCGA-BRCA HiSeqV2 is log2(RSEM normalised count + 1): counts, not TPM, so the linear matrix used for EPIC / quanTIseq / CIBERSORT is not length-corrected. |
| `33_cibersort_lm22_variants.R` | CIBERSORT/LM22 goodness-of-fit under both quantile-normalisation settings, with the published doPerm null. |
| `34_collagenfree_variants.R` | CIRCULARITY GUARD. Several published fibroblast signatures literally contain the outcome variable: MCP-counter's Fibroblasts score is the mean of 8 markers of which COL1A1, COL3A1, COL6A1 and COL6A2 are four; the two … |
| `35_consensus_mediator.R` | Build two method-agnostic consensus fibroblast axes and append them to the mediator set, so the mediation can also be reported against an estimate that does not belong to any single method. |
| `42_deconv2_external_refs.R` | Two EXTERNAL, independently generated estimates of tumour composition for the same TCGA-BRCA tumours, neither produced by me: (1) Thorsson et al. 2018 (Immunity 48:812) PanImmune official CIBERSORT LM22 relative fractions, run by … |
| `43_deconv2_newmethods.R` | Further deconvolution methods, so that the mediation analysis spans several algorithm families and not only several signature sets. |
| `44_deconv2_pathology_stroma.py` | Retrieve the PATHOLOGIST-SCORED slide annotations for every TCGA-BRCA case from the GDC API (percent_stromal_cells, percent_tumor_cells, percent_lymphocyte_infiltration, percent_necrosis, percent_normal_cells). |
| `45_deconv2_assemble.R` | (1) run the collagen-free variants of the new least-squares / dtangle methods, so that every new fibroblast estimate also exists in a form whose signature cannot contain COL1A1 or COL3A1; (2) assemble EVERY fibroblast/stromal … |
| `46_deconv2_mediation.R` | THE DECISIVE RE-TEST, extended. Repeat the causal-mediation analysis of the TF -> collagen and miR-29 -\| collagen edges once per fibroblast/stromal estimate, now including * five new algorithm families run on the same Wu breast … |
| `47_deconv2_bayesprism.R` | BayesPrism (v2.2.3) with the Wu et al. 2021 breast single-cell atlas as reference, compared with the InstaPrism estimates. |
| `48_deconv2_celltype_correlations.R` | Correlate EVERY cell-type fraction from EVERY method with COL1A1, COL3A1, the miR-29 family and the two TFs, so the paper can state which compartment each component of the circuit tracks. |
| `49_deconv2_benchmark_newmethods.R` | An INDEPENDENT accuracy benchmark that includes the new algorithm families. |
| `50_deconv2_matched_collagen_control.R` | Matched collagen controls for deconvolution: the same 641-gene Wu signature without its eight COL* genes, or without COL1A1, run through the same estimators. |
| `51_deconv2_contradictions.R` | Which methods disagree, and why. For every mediation axis, classify each mediator's verdict, find the minority verdicts, and test whether the disagreement is explained by (i) how strongly the mediator itself tracks the outcome … |
| `52_deconv2_methylation_epidish.R` | A FIBROBLAST ESTIMATE FROM A DIFFERENT ASSAY. |
| `53_deconv2_consolidate.R` | Summary table of the multi-method deconvolution analysis. |
