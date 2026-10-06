# ==============================================================================
# MASTER SCRIPT 02: miRNA Target Gene Set Enrichment Analysis (GSEA)
# ==============================================================================
#
# CONSOLIDATES (2 original scripts):
#   - miRNA Target Gene Set Enrichment Analysis.R  → Section 4 (single miRNA mode)
#   - Automated miRNA Target Gene Expression Analysis.R → Section 5 (multi-miRNA loop)
#
# PURPOSE:
#   For one or more miRNAs of interest, this script:
#     1. Loads TCGA BRCA expression and phenotype data from UCSC Xena
#     2. Performs differential expression (limma) between Tumor and Normal
#     3. Builds a ranked gene list from the t-statistics
#     4. Assembles a custom gene set from the miRNA's known MSigDB targets
#     5. Runs GSEA (fgsea) to test whether those targets are collectively
#        enriched in the direction expected for each miRNA
#     6. Generates per-miRNA plots: GSEA enrichment curve, volcano,
#        and expression barplot of validated targets
#
# TWO MODES (set in SECTION 0):
#   ANALYSIS_MODE = "SINGLE"  — Analyse one miRNA (MIRNA_OF_INTEREST)
#   ANALYSIS_MODE = "MULTI"   — Loop over all miRNAs in MIRNA_LIST
#
# REQUIRED INPUT FILES:
#   EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
#     → Gene expression matrix from UCSC Xena
#   TCGA_phenotype_denseDataOnlyDownload.tsv.gz
#     → Sample phenotype annotations from UCSC Xena
#
# OUTPUTS (written to ./miRNA_GSEA_Results/{miRNA_name}/):
#   {miRNA}_GSEA_Enrichment_Plot.png  — Classic GSEA enrichment curve
#   {miRNA}_Volcano_Plot.png          — Volcano plot of all genes
#   {miRNA}_Target_Expression.png     — Target gene expression barplot
#   {miRNA}_target_DE_results.csv     — DEX table for target genes only
#
# DEPENDENCIES:
#   CRAN:         dplyr, ggplot2, msigdbr, tidyr, data.table, tibble, ggrepel
#   Bioconductor: limma, fgsea
#
# AUTHORS:  Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad
# VERSION:  1.0 (consolidated)
# ==============================================================================


# ==============================================================================
# SECTION 0: Configuration  ← EDIT THESE
# ==============================================================================

ANALYSIS_MODE <- "MULTI"         # "SINGLE" | "MULTI"

# For SINGLE mode: the one miRNA to analyse
MIRNA_OF_INTEREST <- "hsa-miR-29a-3p"

# For MULTI mode: list of miRNAs to loop over
MIRNA_LIST <- c(
  "hsa-miR-21-5p",
  "hsa-miR-124-3p",
  "hsa-miR-34a-5p",
  "hsa-miR-101-3p",
  "hsa-let-7b-5p",
  "hsa-let-7e-5p",
  "hsa-miR-29a-3p",
  "hsa-miR-29b-3p",
  "hsa-miR-29c-3p"
)

# Input file paths (Xena Pan-Cancer)
GENE_FILE  <- "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
PHENO_FILE <- "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"

# Cancer type filter
CANCER_TYPE     <- "breast invasive carcinoma"
CANCER_TYPE_COL <- "_primary_disease"
SAMPLE_TYPE_COL <- "sample_type"
TUMOR_LABEL     <- "Primary Tumor"
NORMAL_LABEL    <- "Solid Tissue Normal"

# MSigDB target gene set prefix (C3 miRNA targets collection)
MIRNA_MSIG_CATEGORY    <- "C3"
MIRNA_MSIG_SUBCATEGORY <- "MIR:MIR_Legacy"   # legacy miRNA target gene sets

# Output
OUT_BASE <- "miRNA_GSEA_Results"


# ==============================================================================
# SECTION 1: Package Installation and Loading
# ==============================================================================

cat("=== miRNA TARGET GSEA PIPELINE ===\n")
cat("Mode:", ANALYSIS_MODE, "\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("dplyr", "ggplot2", "msigdbr", "tidyr", "data.table", "tibble", "ggrepel")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

bioc_pkgs <- c("limma", "fgsea")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

dir.create(OUT_BASE, showWarnings = FALSE)


# ==============================================================================
# SECTION 2: Load and Prepare BRCA Expression Data
#   (shared by both modes — data loaded once, reused for each miRNA)
# ==============================================================================

cat("--- [1] Loading BRCA expression and phenotype data ---\n")

# Load phenotype and filter for BRCA Tumor/Normal
pheno <- data.table::fread(PHENO_FILE)

brca_pheno <- pheno %>%
  dplyr::filter(.data[[CANCER_TYPE_COL]] == CANCER_TYPE,
                .data[[SAMPLE_TYPE_COL]] %in% c(TUMOR_LABEL, NORMAL_LABEL)) %>%
  dplyr::mutate(
    sample_type = factor(.data[[SAMPLE_TYPE_COL]], levels = c(NORMAL_LABEL, TUMOR_LABEL))
  ) %>%
  dplyr::distinct(sample, .keep_all = TRUE)

cat("BRCA samples — Tumor:", sum(brca_pheno[[SAMPLE_TYPE_COL]] == TUMOR_LABEL),
    "| Normal:", sum(brca_pheno[[SAMPLE_TYPE_COL]] == NORMAL_LABEL), "\n")

# Load only the BRCA columns from the large expression matrix
cat("Loading gene expression (BRCA columns only)...\n")
all_expr_headers <- names(data.table::fread(GENE_FILE, nrows = 0))
brca_cols        <- intersect(all_expr_headers[-1], brca_pheno$sample)
expr_raw         <- data.table::fread(GENE_FILE, select = c(all_expr_headers[1], brca_cols))

# Handle duplicate column names
orig_cols   <- colnames(expr_raw)
unique_cols <- make.unique(orig_cols, sep = ".")
if (any(orig_cols != unique_cols)) colnames(expr_raw) <- unique_cols

# Build expression matrix
expr_mat <- expr_raw %>%
  dplyr::rename(Gene = 1) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

# Align phenotype to matrix columns
brca_pheno <- brca_pheno %>% dplyr::filter(sample %in% colnames(expr_mat))
expr_mat   <- expr_mat[, brca_pheno$sample]
cat("Expression matrix built:", nrow(expr_mat), "genes ×", ncol(expr_mat), "samples\n\n")


# ==============================================================================
# SECTION 3: Differential Expression (Limma — run once for all miRNAs)
# ==============================================================================

cat("--- [2] Running limma differential expression (Tumor vs Normal) ---\n")

# ------------------------------------------------------------------------------
# run_brca_limma()
#   Performs two-group limma analysis (Tumor vs. Normal) on the BRCA matrix.
#   Returns full topTable result, sorted by t-statistic (ascending).
#
# Returns: data frame with Gene, logFC, t, P.Value, adj.P.Val
# ------------------------------------------------------------------------------
run_brca_limma <- function() {
  groups  <- factor(brca_pheno[[SAMPLE_TYPE_COL]], levels = c(NORMAL_LABEL, TUMOR_LABEL))
  design  <- model.matrix(~ groups)
  fit     <- limma::lmFit(expr_mat, design)
  fit     <- limma::eBayes(fit)
  res     <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "t") %>%
    tibble::rownames_to_column("Gene")
  return(res)
}

de_results  <- run_brca_limma()
ranked_full <- setNames(de_results$t, de_results$Gene)
ranked_full <- sort(ranked_full[!is.na(ranked_full)], decreasing = TRUE)
cat("Ranked gene list:", length(ranked_full), "genes\n\n")


# ==============================================================================
# SECTION 4: miRNA-Specific GSEA and Visualisation
#   (each miRNA gets its own sub-directory with all three plot types)
# ==============================================================================

# Determine which miRNAs to run
mirnas_to_run <- if (ANALYSIS_MODE == "SINGLE") MIRNA_OF_INTEREST else MIRNA_LIST

# Load ALL MSigDB miRNA target gene sets once
cat("--- [3] Loading MSigDB miRNA target gene sets ---\n")
tryCatch({
  msig_mirna <- msigdbr::msigdbr(
    species     = "Homo sapiens",
    category    = MIRNA_MSIG_CATEGORY,
    subcategory = MIRNA_MSIG_SUBCATEGORY
  )
  cat("Loaded", length(unique(msig_mirna$gs_name)), "miRNA gene sets from MSigDB\n\n")
}, error = function(e) {
  cat("WARNING: Could not load MSigDB miRNA sets:", conditionMessage(e), "\n")
  msig_mirna <<- NULL
})


# ------------------------------------------------------------------------------
# mirna_to_msig_name()
#   Attempts to match a standard miRNA name (hsa-miR-XXX-Yp) to the MSigDB
#   gene set naming convention (MIR-XXX or MIRLET-XXX).
#   Returns the closest matching gene set name, or NULL if not found.
#
# Args:
#   mirna_name   — character
#   msig_df      — MSigDB data frame with gs_name column
# Returns: character gene set name or NULL
# ------------------------------------------------------------------------------
mirna_to_msig_name <- function(mirna_name, msig_df) {
  if (is.null(msig_df)) return(NULL)

  # Normalise: strip "hsa-", upper-case, convert "-" to "_"
  clean <- toupper(mirna_name)
  clean <- gsub("^HSA-", "", clean)
  clean <- gsub("-", "", clean)

  # Try prefix match in the gene set names
  all_sets <- unique(msig_df$gs_name)
  hits     <- grep(gsub("[0-9]+[ABab]?.*$", "", clean), all_sets, value = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) return(NULL)
  return(hits[1])
}


cat("--- [4] Running per-miRNA GSEA ---\n")

for (mirna in mirnas_to_run) {
  cat("\n> Processing:", mirna, "\n")

  mirna_dir <- file.path(OUT_BASE, gsub("[^[:alnum:]_-]", "_", mirna))
  dir.create(mirna_dir, showWarnings = FALSE)

  # --- Find matching MSigDB gene set ---
  gsname <- mirna_to_msig_name(mirna, msig_mirna)

  if (is.null(gsname)) {
    cat("  WARNING: No MSigDB gene set found for", mirna, "— skipping GSEA\n")
    # Still write the target expression plot if DE data exists
  } else {
    cat("  Matched MSigDB gene set:", gsname, "\n")

    # Assemble target gene set
    mirna_targets <- msig_mirna %>%
      dplyr::filter(gs_name == gsname) %>%
      dplyr::pull(gene_symbol) %>%
      unique()
    cat("  Target genes in set:", length(mirna_targets), "\n")

    gs_list <- list(mirna_targets)
    names(gs_list) <- gsname

    # --- GSEA ---
    tryCatch({
      gsea_res <- fgsea::fgsea(
        pathways    = gs_list,
        stats       = ranked_full,
        minSize     = 10,
        maxSize     = 5000,
        nPermSimple = 10000
      )

      cat("  NES:", round(gsea_res$NES, 3),
          "| padj:", signif(gsea_res$padj, 3), "\n")

      # Classic GSEA enrichment plot
      p_gsea <- fgsea::plotEnrichment(gs_list[[1]], ranked_full) +
        ggplot2::labs(
          title    = paste("GSEA Enrichment —", mirna),
          subtitle = paste0("Targets of ", mirna, " | NES = ", round(gsea_res$NES, 3),
                            " | padj = ", signif(gsea_res$padj, 3)),
          x        = "Rank in Ordered Gene List",
          y        = "Enrichment Score"
        ) +
        ggplot2::theme_bw(base_size = 12)

      ggplot2::ggsave(
        file.path(mirna_dir, paste0(mirna, "_GSEA_Enrichment_Plot.png")),
        p_gsea, width = 10, height = 6, dpi = 300
      )

      write.csv(as.data.frame(gsea_res),
                file.path(mirna_dir, paste0(mirna, "_GSEA_result.csv")),
                row.names = FALSE)
    }, error = function(e) cat("  ERROR in GSEA:", conditionMessage(e), "\n"))
  }

  # --- Volcano plot ---
  # Highlights targets of this miRNA in red
  tryCatch({
    vol_df <- de_results %>%
      dplyr::mutate(
        is_target  = Gene %in% (if (!is.null(gsname)) {
          msig_mirna %>% dplyr::filter(gs_name == gsname) %>% dplyr::pull(gene_symbol)
        } else character(0)),
        sig_label  = dplyr::case_when(
          is_target & adj.P.Val < 0.05 ~ "Target (sig.)",
          is_target                    ~ "Target",
          TRUE                         ~ "Other"
        ),
        label_gene = ifelse(is_target & adj.P.Val < 0.05, Gene, NA)
      )

    p_vol <- ggplot2::ggplot(vol_df, ggplot2::aes(x = logFC, y = -log10(adj.P.Val + 1e-300),
                                                   colour = sig_label)) +
      ggplot2::geom_point(alpha = 0.4, size = 0.8) +
      ggplot2::scale_colour_manual(
        values = c("Target (sig.)" = "red", "Target" = "orange", "Other" = "grey70")
      ) +
      ggrepel::geom_text_repel(ggplot2::aes(label = label_gene), size = 2.5,
                               max.overlaps = 20, na.rm = TRUE) +
      ggplot2::labs(
        title  = paste("Volcano Plot —", mirna),
        x      = "Log2 Fold Change (Tumor vs Normal)",
        y      = "-log10(adj. p-value)",
        colour = "Gene Category"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
      ggplot2::geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40")

    ggplot2::ggsave(
      file.path(mirna_dir, paste0(mirna, "_Volcano_Plot.png")),
      p_vol, width = 10, height = 8, dpi = 300
    )
  }, error = function(e) cat("  ERROR in volcano plot:", conditionMessage(e), "\n"))

  # --- Target expression barplot (top 20 significant targets) ---
  tryCatch({
    if (!is.null(gsname)) {
      targ_genes <- msig_mirna %>% dplyr::filter(gs_name == gsname) %>% dplyr::pull(gene_symbol)
    } else {
      targ_genes <- character(0)
    }

    if (length(targ_genes) > 0) {
      targ_de <- de_results %>%
        dplyr::filter(Gene %in% targ_genes) %>%
        dplyr::arrange(adj.P.Val) %>%
        head(20)

      write.csv(targ_de,
                file.path(mirna_dir, paste0(mirna, "_target_DE_results.csv")),
                row.names = FALSE)

      if (nrow(targ_de) > 0) {
        p_bar <- ggplot2::ggplot(targ_de,
                                 ggplot2::aes(x = reorder(Gene, logFC), y = logFC,
                                              fill = logFC > 0)) +
          ggplot2::geom_col() +
          ggplot2::scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                                     name = "", labels = c("Down in Tumor", "Up in Tumor")) +
          ggplot2::coord_flip() +
          ggplot2::labs(
            title    = paste("Top Target Gene Expression —", mirna),
            subtitle = "Top 20 targets by adj.p-value",
            x        = NULL, y = "Log2 Fold Change (Tumor vs Normal)"
          ) +
          ggplot2::theme_minimal(base_size = 11)

        ggplot2::ggsave(
          file.path(mirna_dir, paste0(mirna, "_Target_Expression.png")),
          p_bar, width = 9, height = 7, dpi = 300
        )
      }
    }
  }, error = function(e) cat("  ERROR in target barplot:", conditionMessage(e), "\n"))

  cat("  Done:", mirna, "\n")
}

cat("\n=== miRNA GSEA PIPELINE COMPLETE ===\n")
cat("Results saved to:", OUT_BASE, "\n")
