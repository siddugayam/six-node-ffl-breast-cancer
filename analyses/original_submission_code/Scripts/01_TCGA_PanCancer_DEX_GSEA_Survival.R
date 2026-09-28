# ==============================================================================
# MASTER SCRIPT 01: TCGA Pan-Cancer Data Audit, Differential Expression,
#                   Hallmark GSEA, and Survival Analysis
# ==============================================================================
#
# CONSOLIDATES (7 original scripts):
#   - stat_and_header_chek.R             → Section 1: File audit
#   - First_Stage_verification.R         → Section 1: Sample overlap check
#   - ID_Genes_TFs_miRNAs.R              → Section 2: Data filtering & TF ID
#   - TCGA_PAN_CANCER_BRCA_HALLMARKS_XENA_DEX_DATA.R → Section 3: BRCA DEX
#   - dessese_wise_results.R             → Section 4: Disease-wise DEX loop
#   - Untitled.R (complete pipeline)     → Sections 5–7: GSEA + survival
#   - TGCA_PAN_CAN_BRCA_METAVSNORM_HALLMARK_SURV.R → Section 3 (mode toggle)
#   NOTE: DESEQ_BRCA_TCGA_PANCAN-1.R and TCGA_PAN_CANCER_BRCA_HALLMARKS_XENA.R
#         were superseded by the XENA_DEX_DATA version (includes dup fix) and merged.
#
# PURPOSE:
#   Full pipeline from raw Xena files → differential expression → Hallmark GSEA
#   → survival analysis across all TCGA cancer types.
#   Can be run in three modes (see SECTION 0):
#     "BRCA_TvN"   — Breast cancer: Primary Tumor vs. Solid Tissue Normal
#     "BRCA_MvN"   — Breast cancer: Metastatic vs. Solid Tissue Normal
#     "PAN_CANCER" — All TCGA cancer types with sufficient normal samples
#
# WORKFLOW:
#   [1] File audit  →  [2] Load & clean data  →  [3] BRCA or Pan-Cancer DEX
#   → [4] Disease-wise limma loop  →  [5] Hallmark GSEA  →  [6] Survival
#   → [7] Pan-cancer summary heatmap
#
# REQUIRED INPUT FILES (place in working directory):
#   EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
#     → Gene expression matrix (log2-normalized), from UCSC Xena Pan-Cancer
#   pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena
#     → miRNA expression matrix (same platform), from UCSC Xena
#   TCGA_phenotype_denseDataOnlyDownload.tsv.gz
#     → Sample phenotype/type annotations, from UCSC Xena
#   Survival_SupplementalTable_S1_20171025_xena_sp
#     → Clinical survival data (OS time, OS event), from UCSC Xena
#   DatabaseExtract_v_1.01.csv
#     → Human TF list (humantfs.ccbr.utoronto.ca), column "Is TF?" = "Yes"
#
# OUTPUTS (written to ./TCGA_Analysis_Results/):
#   00_audit/          — File headers, sample counts, overlap summary
#   01_filtered_data/  — TF/gene/miRNA subsets of pan-cancer matrix
#   02_DEX/            — Per-disease limma DEX tables
#   03_GSEA/           — Per-cancer GSEA NES plots + heatmaps
#   04_survival/       — Kaplan-Meier survival plots per cancer
#   Pan-Cancer_Hallmark_Activity_Heatmap.png  — Summary heatmap
#
# KEY PARAMETERS (edit in SECTION 0):
#   RUN_MODE              — "BRCA_TvN", "BRCA_MvN", or "PAN_CANCER"
#   NORMAL_SAMPLE_THRESHOLD — Minimum normal samples to include a cancer type
#   AUDIT_ONLY            — Set TRUE to only run the file audit (Section 1)
#   SAVE_CLEANED_DATA     — Set TRUE to write cleaned matrices to disk
#
# DEPENDENCIES:
#   CRAN:         dplyr, ggplot2, msigdbr, tidyr, pheatmap, data.table,
#                 tibble, survival, survminer
#   Bioconductor: limma, fgsea, GSVA
#
# AUTHORS:  Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad
# VERSION:  1.0 (consolidated)
# ==============================================================================


# ==============================================================================
# SECTION 0: Configuration  ← EDIT THESE
# ==============================================================================

RUN_MODE               <- "PAN_CANCER"  # "BRCA_TvN" | "BRCA_MvN" | "PAN_CANCER"
NORMAL_SAMPLE_THRESHOLD <- 25           # min normal samples to include a cancer type
AUDIT_ONLY             <- FALSE         # TRUE = only run Section 1 (file audit)
SAVE_CLEANED_DATA      <- TRUE          # TRUE = write cleaned matrices to disk
TOP_PATHWAYS_TO_PLOT   <- 15           # top N pathways in individual GSEA plots

# Input file paths
GENE_FILE    <- "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
MIRNA_FILE   <- "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
PHENO_FILE   <- "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
SURVIVAL_FILE <- "Survival_SupplementalTable_S1_20171025_xena_sp"
TF_LIST_FILE <- "DatabaseExtract_v_1.01.csv"

# Column names in phenotype file
CANCER_TYPE_COL <- "_primary_disease"
SAMPLE_TYPE_COL <- "sample_type"

# Output
OUT_BASE <- "TCGA_Analysis_Results"


# ==============================================================================
# SECTION 1: Package Installation and Loading
# ==============================================================================

cat("=== TCGA PAN-CANCER PIPELINE ===\n")
cat("Mode:", RUN_MODE, "\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("dplyr", "ggplot2", "msigdbr", "tidyr", "pheatmap",
               "data.table", "tibble", "survival", "survminer")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

bioc_pkgs <- c("limma", "fgsea", "GSVA")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

dir.create(OUT_BASE, showWarnings = FALSE)


# ==============================================================================
# SECTION 2: File Audit and Sample Overlap Check
#   (consolidates stat_and_header_chek.R + First_Stage_verification.R)
# ==============================================================================

cat("--- [1] File Audit ---\n")

audit_dir <- file.path(OUT_BASE, "00_audit")
dir.create(audit_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# get_sample_type()
#   Decodes TCGA barcode characters 14-15 to sample type classification.
#   Codes 01-09 = Tumor; 10-19 = Normal; otherwise Unknown.
# Args: barcodes — character vector of TCGA sample barcodes
# Returns: character vector ("Tumor" | "Normal" | "Unknown")
# ------------------------------------------------------------------------------
get_sample_type <- function(barcodes) {
  type_bits <- as.numeric(substr(barcodes, 14, 15))
  dplyr::case_when(
    type_bits >= 1  & type_bits <= 9  ~ "Tumor",
    type_bits >= 10 & type_bits <= 19 ~ "Normal",
    TRUE                              ~ "Unknown"
  )
}

audit_results <- list()

if (file.exists(GENE_FILE)) {
  gene_h <- names(data.table::fread(GENE_FILE, nrows = 0))
  audit_results$gene_samples <- length(gene_h) - 1
  cat("Gene expression — Samples:", audit_results$gene_samples, "| ID col:", gene_h[1], "\n")
} else {
  cat("WARNING: Gene expression file not found:", GENE_FILE, "\n")
}

if (file.exists(MIRNA_FILE)) {
  mir_h <- names(data.table::fread(MIRNA_FILE, nrows = 0))
  audit_results$mirna_samples <- length(mir_h) - 1
  cat("miRNA expression — Samples:", audit_results$mirna_samples, "| ID col:", mir_h[1], "\n")

  if (exists("gene_h")) {
    common_gm <- intersect(gene_h[-1], mir_h[-1])
    audit_results$gene_mirna_overlap <- length(common_gm)
    cat("Gene-miRNA sample overlap:", length(common_gm), "\n")

    type_tbl <- table(get_sample_type(common_gm))
    cat("Sample types in overlap — Tumor:", type_tbl["Tumor"],
        "| Normal:", type_tbl["Normal"], "\n")
  }
} else {
  cat("WARNING: miRNA file not found:", MIRNA_FILE, "\n")
}

if (file.exists(PHENO_FILE)) {
  pheno_preview <- data.table::fread(PHENO_FILE, nrows = 5)
  audit_results$pheno_cols <- paste(names(pheno_preview), collapse = ", ")
  cat("Phenotype file columns:", audit_results$pheno_cols, "\n")
}

write.csv(as.data.frame(audit_results), file.path(audit_dir, "file_audit_summary.csv"),
          row.names = TRUE)
cat("Audit summary saved to 00_audit/\n\n")

if (AUDIT_ONLY) {
  cat("AUDIT_ONLY = TRUE. Stopping after Section 1.\n")
  stop("(Intentional stop — set AUDIT_ONLY = FALSE to run full pipeline)")
}


# ==============================================================================
# SECTION 3: Load, Clean, and Filter Pan-Cancer Data
#   (consolidates ID_Genes_TFs_miRNAs.R + data cleaning from Untitled.R)
# ==============================================================================

cat("--- [2] Loading and cleaning data ---\n")

# --- 3A. Load TF reference list ---
if (file.exists(TF_LIST_FILE)) {
  tf_meta <- data.table::fread(TF_LIST_FILE)
  human_tf_symbols <- tf_meta[`Is TF?` == "Yes", `HGNC symbol`]
  cat("TF list loaded:", length(human_tf_symbols), "human TFs\n")
} else {
  human_tf_symbols <- character(0)
  cat("WARNING: TF list file not found. TF separation will be skipped.\n")
}

# --- 3B. Load phenotype data ---
cat("Loading phenotype data...\n")
pheno <- data.table::fread(PHENO_FILE)

# --- 3C. Load expression data (handle duplicate column names) ---
cat("Loading gene expression matrix (may take a few minutes)...\n")
expression_data <- data.table::fread(GENE_FILE)

orig_cols   <- colnames(expression_data)
unique_cols <- make.unique(orig_cols, sep = ".")
if (any(orig_cols != unique_cols)) {
  cat("Duplicate column names found and repaired.\n")
  colnames(expression_data) <- unique_cols
}

# --- 3D. Load survival data ---
if (file.exists(SURVIVAL_FILE)) {
  cat("Loading survival data...\n")
  clinical_data <- data.table::fread(SURVIVAL_FILE)
} else {
  clinical_data <- NULL
  cat("WARNING: Survival file not found. Survival analysis will be skipped.\n")
}

# --- 3E. Identify common samples across all files ---
cat("Identifying common samples across all data files...\n")
expression_samples  <- colnames(expression_data)[-1]
sample_info_samples <- pheno$sample
survival_samples    <- if (!is.null(clinical_data)) clinical_data$sample else expression_samples

core_samples <- Reduce(intersect, list(expression_samples, sample_info_samples, survival_samples))
cat("Core samples (present in all files):", length(core_samples), "\n")

# --- 3F. De-duplicate: one sample per patient per group ---
final_pheno <- pheno %>%
  dplyr::filter(sample %in% core_samples) %>%
  dplyr::mutate(Patient_ID = substr(sample, 1, 12)) %>%
  dplyr::group_by(Patient_ID, .data[[SAMPLE_TYPE_COL]]) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

cat("Samples after de-duplication:", nrow(final_pheno), "\n")
final_clean_samples <- final_pheno$sample
final_survival_data <- if (!is.null(clinical_data)) dplyr::filter(clinical_data, sample %in% final_clean_samples) else NULL

# --- 3G. Build final expression matrix (remove zero-expression genes) ---
cat("Building cleaned expression matrix...\n")
expr_df <- expression_data %>%
  dplyr::select(sample, dplyr::all_of(intersect(final_clean_samples, colnames(expression_data)[-1]))) %>%
  dplyr::rename(Gene = sample)

if (any(duplicated(expr_df$Gene))) {
  cat("Averaging duplicated gene symbols...\n")
  expr_df <- expr_df %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop")
}

expr_mat <- expr_df %>% tibble::column_to_rownames("Gene") %>% as.matrix()
zero_genes  <- apply(expr_mat, 1, function(r) any(r == 0))
expr_mat    <- expr_mat[!zero_genes, ]
cat("Genes retained after zero-expression filter:", nrow(expr_mat), "\n")

# --- 3H. Separate into TF / Gene / miRNA subsets ---
if (length(human_tf_symbols) > 0) {
  cat("Separating expression data into TFs and target genes...\n")
  id_col <- names(expression_data)[1]
  gene_ids  <- names(data.table::fread(GENE_FILE, nrows = 0))[-1]  # reuse cached headers
  tf_gene_names    <- intersect(rownames(expr_mat), human_tf_symbols)
  nontf_gene_names <- setdiff(rownames(expr_mat), human_tf_symbols)
  cat("TFs identified:", length(tf_gene_names), "| Non-TF genes:", length(nontf_gene_names), "\n")
}

# --- 3I. Save cleaned data ---
if (SAVE_CLEANED_DATA) {
  filt_dir <- file.path(OUT_BASE, "01_filtered_data")
  dir.create(filt_dir, showWarnings = FALSE)
  write.csv(as.data.frame(expr_mat), file.path(filt_dir, "final_expression_matrix.csv"))
  write.csv(final_pheno,             file.path(filt_dir, "final_phenotype_data.csv"),    row.names = FALSE)
  if (!is.null(final_survival_data))
    write.csv(final_survival_data,   file.path(filt_dir, "final_survival_data.csv"),     row.names = FALSE)
  cat("Cleaned data saved to 01_filtered_data/\n")
}

# --- 3J. Generate sample count summary table ---
counts_before <- pheno %>%
  dplyr::filter(.data[[SAMPLE_TYPE_COL]] %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  dplyr::group_by(.data[[CANCER_TYPE_COL]], .data[[SAMPLE_TYPE_COL]]) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = SAMPLE_TYPE_COL, values_from = count, values_fill = 0) %>%
  dplyr::rename(Disease = all_of(CANCER_TYPE_COL),
                Tumor_Before = `Primary Tumor`, Normal_Before = `Solid Tissue Normal`)

counts_after <- final_pheno %>%
  dplyr::filter(.data[[SAMPLE_TYPE_COL]] %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  dplyr::group_by(.data[[CANCER_TYPE_COL]], .data[[SAMPLE_TYPE_COL]]) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = SAMPLE_TYPE_COL, values_from = count, values_fill = 0) %>%
  dplyr::rename(Disease = all_of(CANCER_TYPE_COL),
                Tumor_After = `Primary Tumor`, Normal_After = `Solid Tissue Normal`)

summary_table <- dplyr::full_join(counts_before, counts_after, by = "Disease") %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~tidyr::replace_na(., 0)),
                Meets_Threshold = ifelse(Normal_After >= NORMAL_SAMPLE_THRESHOLD, "Yes", "No")) %>%
  dplyr::arrange(desc(Tumor_After))

write.csv(summary_table, file.path(OUT_BASE, "sample_count_summary.csv"), row.names = FALSE)
cat("Sample count summary saved.\n\n")


# ==============================================================================
# SECTION 4: Disease-wise Limma Differential Expression
#   (consolidates dessese_wise_results.R + BRCA DEX from XENA_DEX_DATA version)
# ==============================================================================

cat("--- [3] Differential Expression Analysis ---\n")

dex_dir <- file.path(OUT_BASE, "02_DEX")
dir.create(dex_dir, showWarnings = FALSE)

# Determine which cancers and comparison to run based on RUN_MODE
if (RUN_MODE == "BRCA_TvN") {
  cancers_to_run  <- "breast invasive carcinoma"
  group_levels    <- c("Solid Tissue Normal", "Primary Tumor")
  comparison_name <- "PrimaryTumor_vs_Normal"
} else if (RUN_MODE == "BRCA_MvN") {
  cancers_to_run  <- "breast invasive carcinoma"
  group_levels    <- c("Solid Tissue Normal", "Metastatic")
  comparison_name <- "Metastatic_vs_Normal"
} else {
  # PAN_CANCER: run all cancers with enough samples
  cancers_to_run  <- summary_table %>% dplyr::filter(Meets_Threshold == "Yes") %>% dplyr::pull(Disease)
  group_levels    <- c("Solid Tissue Normal", "Primary Tumor")
  comparison_name <- "PrimaryTumor_vs_Normal"
}

cat("Cancers to analyse:", length(cancers_to_run), "\n")

# ------------------------------------------------------------------------------
# run_limma_de()
#   Runs a two-group limma + eBayes differential expression analysis.
#   Contrast: group 2 - group 1 (as defined by group_levels).
#
# Args:
#   mat        — Numeric expression matrix (genes × samples)
#   groups     — Factor vector with levels = group_levels
#   out_file   — Path to write the DEX results CSV
#
# Returns: topTable data frame with logFC, t, P.Value, adj.P.Val
# ------------------------------------------------------------------------------
run_limma_de <- function(mat, groups, out_file) {
  design         <- model.matrix(~ 0 + groups)
  colnames(design) <- make.names(levels(groups))
  contrast_str   <- paste0(make.names(levels(groups))[2], "-", make.names(levels(groups))[1])
  cont_matrix    <- limma::makeContrasts(contrasts = contrast_str, levels = design)

  fit  <- limma::lmFit(mat, design)
  fit  <- limma::contrasts.fit(fit, cont_matrix)
  fit  <- limma::eBayes(fit)
  res  <- limma::topTable(fit, coef = 1, number = Inf, adjust.method = "BH") %>%
    tibble::rownames_to_column("Gene")
  data.table::fwrite(res, out_file)
  return(res)
}

all_dex_results <- list()

for (cancer in cancers_to_run) {
  cat("  DEX:", cancer, "\n")

  cancer_pheno <- final_pheno %>%
    dplyr::filter(.data[[CANCER_TYPE_COL]] == cancer,
                  .data[[SAMPLE_TYPE_COL]] %in% group_levels)

  if (nrow(cancer_pheno) < 10) {
    cat("    Skipping (fewer than 10 matched samples)\n")
    next
  }

  avail_samples <- intersect(cancer_pheno$sample, colnames(expr_mat))
  cancer_mat    <- expr_mat[, avail_samples]
  cancer_groups <- factor(cancer_pheno[[SAMPLE_TYPE_COL]][cancer_pheno$sample %in% avail_samples],
                          levels = group_levels)

  cancer_clean <- gsub("[^[:alnum:]]", "_", cancer)
  out_file     <- file.path(dex_dir, paste0(cancer_clean, "_", comparison_name, "_DEX.csv"))

  tryCatch({
    res <- run_limma_de(cancer_mat, cancer_groups, out_file)
    all_dex_results[[cancer]] <- res
    cat("    Saved:", nrow(res), "genes\n")
  }, error = function(e) cat("    ERROR:", conditionMessage(e), "\n"))
}


# ==============================================================================
# SECTION 5: Hallmark GSEA per Cancer Type
#   (consolidates Comprehensive Pan_Cancer Hallmark GSEA Workflow.R + Untitled.R)
# ==============================================================================

cat("\n--- [4] Hallmark GSEA ---\n")

gsea_dir <- file.path(OUT_BASE, "03_GSEA")
dir.create(gsea_dir, showWarnings = FALSE)

# Load Hallmark gene sets once (MSigDB H collection)
hallmark_df   <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

all_gsea_results <- list()

for (cancer in names(all_dex_results)) {
  cat("  GSEA:", cancer, "\n")
  de_res <- all_dex_results[[cancer]]

  # Build ranked gene list (ranked by moderated t-statistic)
  ranked <- de_res$t
  names(ranked) <- de_res$Gene
  ranked <- sort(ranked, decreasing = TRUE)
  ranked <- ranked[!is.na(ranked)]

  tryCatch({
    gsea_res <- fgsea::fgsea(
      pathways   = hallmark_list,
      stats      = ranked,
      minSize    = 15,
      maxSize    = 500,
      nPermSimple = 1000
    )
    gsea_res$cancer_type    <- cancer
    all_gsea_results[[cancer]] <- gsea_res

    # Individual GSEA bar plot (top activated + suppressed)
    cancer_gsea_dir <- file.path(gsea_dir, gsub("[ /]", "_", cancer))
    dir.create(cancer_gsea_dir, showWarnings = FALSE)

    top_paths <- dplyr::bind_rows(
      gsea_res %>% dplyr::filter(NES > 0) %>% dplyr::arrange(padj) %>% head(TOP_PATHWAYS_TO_PLOT),
      gsea_res %>% dplyr::filter(NES < 0) %>% dplyr::arrange(padj) %>% head(TOP_PATHWAYS_TO_PLOT)
    )

    p_gsea <- ggplot2::ggplot(top_paths, ggplot2::aes(reorder(pathway, NES), NES)) +
      ggplot2::geom_col(ggplot2::aes(fill = NES > 0)) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "darkblue"),
                                 name = "", labels = c("Enriched in Normal", "Enriched in Tumor")) +
      ggplot2::labs(x = "Pathway", y = "Normalized Enrichment Score (NES)",
                    title = paste("Hallmark GSEA —", cancer)) +
      ggplot2::theme_minimal(base_size = 10)

    ggplot2::ggsave(
      file.path(cancer_gsea_dir, paste0(gsub("[ /]", "_", cancer), "_Hallmark_GSEA.png")),
      p_gsea, width = 12, height = 10, bg = "white"
    )

    # Leading-edge heatmap for top pathway
    top_pw   <- gsea_res %>% dplyr::arrange(padj) %>% head(1)
    le_genes <- unlist(top_pw$leadingEdge)

    if (length(le_genes) > 2 && cancer %in% unique(final_pheno[[CANCER_TYPE_COL]])) {
      cancer_pheno <- final_pheno %>%
        dplyr::filter(.data[[CANCER_TYPE_COL]] == cancer,
                      .data[[SAMPLE_TYPE_COL]] %in% group_levels)
      cancer_samps <- intersect(cancer_pheno$sample, colnames(expr_mat))
      heat_mat     <- expr_mat[intersect(le_genes, rownames(expr_mat)), cancer_samps, drop = FALSE]

      if (nrow(heat_mat) > 2) {
        ann_col <- data.frame(SampleType = cancer_pheno[[SAMPLE_TYPE_COL]][cancer_pheno$sample %in% cancer_samps])
        rownames(ann_col) <- cancer_samps
        tryCatch(
          pheatmap::pheatmap(
            heat_mat,
            main           = paste("Leading Edge:", gsub("HALLMARK_", "", top_pw$pathway)),
            annotation_col = ann_col,
            scale          = "row",
            show_colnames  = FALSE,
            fontsize_row   = 6,
            filename       = file.path(cancer_gsea_dir, paste0(gsub("[ /]", "_", cancer), "_LeadingEdge_Heatmap.png"))
          ), error = function(e) NULL
        )
      }
    }
    cat("    GSEA saved\n")
  }, error = function(e) cat("    ERROR in GSEA:", conditionMessage(e), "\n"))
}


# ==============================================================================
# SECTION 6: Survival Analysis (GSVA-based pathway stratification)
#   (consolidates TCGA_PAN_CANCR_BRCA_HALLMARK_SURVIVAL.R +
#    TGCA_PAN_CAN_BRCA_METAVSNORM_HALLMARK_SURV.R)
# ==============================================================================

if (!is.null(final_survival_data) && nrow(final_survival_data) > 0) {
  cat("\n--- [5] Survival Analysis ---\n")

  surv_dir <- file.path(OUT_BASE, "04_survival")
  dir.create(surv_dir, showWarnings = FALSE)

  for (cancer in names(all_gsea_results)) {
    cat("  Survival:", cancer, "\n")

    cancer_pheno <- final_pheno %>%
      dplyr::filter(.data[[CANCER_TYPE_COL]] == cancer,
                    .data[[SAMPLE_TYPE_COL]] == "Primary Tumor") %>%
      dplyr::inner_join(final_survival_data, by = "sample")

    if (nrow(cancer_pheno) < 10) { cat("    Skipping (insufficient tumor samples)\n"); next }

    tryCatch({
      # GSVA score for each sample against all Hallmark pathways
      tumor_mat    <- expr_mat[, intersect(cancer_pheno$sample, colnames(expr_mat)), drop = FALSE]
      gsva_param   <- GSVA::gsvaParam(tumor_mat, hallmark_list, kcdf = "Gaussian")
      gsva_scores  <- GSVA::gsva(param = gsva_param, verbose = FALSE)

      # Stratify by activity of the top GSEA pathway
      top_pw_name  <- all_gsea_results[[cancer]] %>%
        dplyr::arrange(padj) %>% dplyr::pull(pathway) %>% dplyr::first()

      pw_scores <- gsva_scores[top_pw_name, ]

      surv_df <- cancer_pheno %>%
        dplyr::mutate(
          gsva_score     = pw_scores[match(sample, names(pw_scores))],
          strata         = ifelse(gsva_score > median(gsva_score, na.rm = TRUE), "High", "Low"),
          OS.time_years  = OS.time / 365.25
        ) %>%
        dplyr::filter(!is.na(strata), !is.na(OS.time_years), !is.na(OS))

      if (nrow(surv_df) < 5) { cat("    Skipping survival (too few valid samples)\n"); next }

      surv_fit <- survival::survfit(
        survival::Surv(OS.time_years, OS) ~ strata, data = surv_df
      )
      km_plot <- survminer::ggsurvplot(
        surv_fit, data = surv_df, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
        legend.labs  = c("High Activity", "Low Activity"),
        legend.title = "Pathway Activity",
        title        = paste("KM Survival —", gsub("HALLMARK_", "", top_pw_name)),
        subtitle     = cancer, xlab = "Time (Years)"
      )
      ggplot2::ggsave(
        file.path(surv_dir, paste0(gsub("[ /]", "_", cancer), "_KM_Survival.png")),
        print(km_plot), width = 8, height = 8, bg = "white"
      )
      cat("    Survival plot saved\n")
    }, error = function(e) cat("    ERROR in survival:", conditionMessage(e), "\n"))
  }
}


# ==============================================================================
# SECTION 7: Pan-Cancer Summary Heatmap
#   (consolidates final section of Untitled.R)
# ==============================================================================

cat("\n--- [6] Pan-Cancer Summary Heatmap ---\n")

if (length(all_gsea_results) > 1) {
  pan_df <- dplyr::bind_rows(all_gsea_results)

  nes_mat <- pan_df %>%
    dplyr::select(pathway, cancer_type, NES) %>%
    tidyr::pivot_wider(names_from = cancer_type, values_from = NES, values_fill = 0) %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()

  rownames(nes_mat) <- gsub("HALLMARK_", "", rownames(nes_mat))

  pheatmap::pheatmap(
    nes_mat,
    main         = "Hallmark Pathway Activity Across TCGA Cancers",
    fontsize_row = 8, fontsize_col = 8,
    cluster_rows = TRUE, cluster_cols = TRUE,
    color        = colorRampPalette(c("blue", "white", "red"))(100),
    filename     = file.path(OUT_BASE, "Pan-Cancer_Hallmark_Activity_Heatmap.png"),
    width = 15, height = 10
  )
  cat("Pan-cancer heatmap saved.\n")
} else {
  cat("Not enough cancers for pan-cancer heatmap.\n")
}

cat("\n=== PIPELINE COMPLETE ===\n")
cat("All results saved to:", OUT_BASE, "\n")
