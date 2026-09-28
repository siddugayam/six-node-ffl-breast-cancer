# ==============================================================================
# SCRIPT 01: TCGA-BRCA miRNA Differential Expression (DESeq2 Preprocessing)
# ==============================================================================
#
# PURPOSE:
#   Runs differential expression analysis on TCGA-BRCA miRNA count data
#   downloaded from the UCSC Xena browser. Produces a formatted CSV file
#   that serves as the validated expression input for the main network
#   analysis pipeline (Script 02).
#
# WORKFLOW POSITION:
#   [01_BRCA_miRNA_DESeq2_Preprocessing.R]  <-- YOU ARE HERE
#            |
#            v  (produces BRCA_miRNA_DEX_results.csv)
#   [02_Signed_FFL_Network_Analysis.R]
#
# INPUTS (must be in working directory):
#   1. TCGA.BRCA.sampleMap_miRNA_HiSeq_gene.tsv
#        - miRNA expression matrix (log2-transformed counts) from UCSC Xena
#        - Download: https://xenabrowser.net/ > TCGA BRCA > miRNA HiSeq
#   2. TCGA.BRCA.sampleMap_BRCA_clinicalMatrix.tsv
#        - Clinical metadata including sample type labels
#        - Download: UCSC Xena > TCGA BRCA > clinical matrix
#
# OUTPUTS:
#   BRCA_miRNA_DEX_results.csv
#     Columns: Gene (miRNA ID), logFC (Log2 Fold Change), adj.P.Val (BH-adjusted p-value)
#     Significance thresholds applied downstream: |logFC| > 1.0 AND adj.P.Val < 0.05
#
# ANALYSIS:
#   - Comparison: Primary Tumor (n ~1097) vs. Solid Tissue Normal (n ~113)
#   - Method: DESeq2 (Love et al., Genome Biology 2014)
#   - Input values are reverse-transformed from log2(counts+1) before DESeq2
#
# DEPENDENCIES:
#   Bioconductor: DESeq2, SummarizedExperiment
#   CRAN:         tidyverse, data.table
#
# AUTHORS:  Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad
# VERSION:  2.0 (cleaned)
# REFERENCE: TCGA-BRCA cohort via UCSC Xena Pan-Cancer Atlas
# ==============================================================================


# ==============================================================================
# SECTION 1: Load Required Packages
# ==============================================================================

cat("--- [1/4] Loading and installing packages ---\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages if missing
bioc_packages <- c("DESeq2", "SummarizedExperiment")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# Install CRAN packages if missing
cran_packages <- c("tidyverse", "data.table")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


# ==============================================================================
# SECTION 2: File Configuration
# ==============================================================================

# --- Input file names (adjust if your downloaded files have different names) ---
COUNTS_FILE <- "TCGA.BRCA.sampleMap_miRNA_HiSeq_gene.tsv"
META_FILE   <- "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix.tsv"

# --- Output file name (consumed by Script 02) ---
OUTPUT_FILE <- "BRCA_miRNA_DEX_results.csv"

# --- Minimum read-count filter (miRNAs with fewer total reads are excluded) ---
MIN_COUNT_FILTER <- 10

# --- Sample type labels used in the clinical matrix ---
TUMOR_LABEL  <- "Primary Tumor"
NORMAL_LABEL <- "Solid Tissue Normal"

cat("--- [2/4] Checking input files ---\n")

if (!file.exists(COUNTS_FILE)) {
  stop(paste(
    "File not found:", COUNTS_FILE,
    "\nPlease download the TCGA-BRCA miRNA HiSeq matrix from",
    "https://xenabrowser.net/ and place it in your working directory."
  ))
}

if (!file.exists(META_FILE)) {
  stop(paste(
    "File not found:", META_FILE,
    "\nPlease download the TCGA-BRCA clinical matrix from",
    "https://xenabrowser.net/ and place it in your working directory."
  ))
}


# ==============================================================================
# SECTION 3: Load and Format Data
# ==============================================================================

cat("--- [3/4] Loading and formatting data ---\n")

# --- Load the counts matrix ---
# The Xena file stores values as log2(counts + 1); DESeq2 requires raw integers.
# We reverse-transform with: round(2^x - 1)
counts_raw <- fread(COUNTS_FILE, data.table = FALSE)

counts_matrix <- counts_raw %>%
  tibble::column_to_rownames("sample") %>%          # miRNA IDs as row names
  as.matrix() %>%
  (function(x) round(2^x - 1))()                   # reverse log2 transform

# --- Load and filter clinical metadata ---
meta_all <- read.delim(META_FILE)

meta_filtered <- meta_all %>%
  dplyr::select(sample = sampleID, sample_type) %>%
  dplyr::filter(sample_type %in% c(TUMOR_LABEL, NORMAL_LABEL)) %>%
  dplyr::mutate(
    sample_type = factor(sample_type, levels = c(NORMAL_LABEL, TUMOR_LABEL))
    # Reference level = Normal (so LogFC > 0 means upregulated in Tumor)
  ) %>%
  dplyr::distinct(sample, .keep_all = TRUE) %>%
  tibble::column_to_rownames("sample")

# --- Align matrix columns to metadata rows ---
# Only retain samples present in both objects, and ensure matching order
common_samples       <- intersect(colnames(counts_matrix), rownames(meta_filtered))
counts_matrix        <- counts_matrix[, common_samples]
meta_filtered        <- meta_filtered[common_samples, , drop = FALSE]

cat(paste(
  "Matched samples:", length(common_samples),
  "(", sum(meta_filtered$sample_type == TUMOR_LABEL), "Tumor /",
  sum(meta_filtered$sample_type == NORMAL_LABEL), "Normal )\n"
))


# ==============================================================================
# SECTION 4: DESeq2 Differential Expression Analysis
# ==============================================================================

cat("--- [4/4] Running DESeq2 ---\n")

# Build the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = meta_filtered,
  design    = ~ sample_type           # single-factor: Tumor vs. Normal
)

# Pre-filter: remove miRNAs with very low total counts across all samples
keep <- rowSums(counts(dds)) >= MIN_COUNT_FILTER
dds  <- dds[keep, ]
cat(paste("miRNAs retained after count filter:", sum(keep), "\n"))

# Run the full DESeq2 pipeline (size factors → dispersions → Wald test)
dds <- DESeq(dds)

# Extract results for Primary Tumor vs. Solid Tissue Normal
res <- results(
  dds,
  name = "sample_type_Primary_Tumor_vs_Solid_Tissue_Normal",
  alpha = 0.05                        # FDR threshold for independent filtering
)

cat(paste(
  "Significant miRNAs (|logFC|>1, padj<0.05):",
  sum(!is.na(res$padj) & abs(res$log2FoldChange) > 1 & res$padj < 0.05),
  "\n"
))

# ==============================================================================
# SECTION 5: Format and Save Output
# ==============================================================================

# Format output to match the column names expected by Script 02
final_dex <- as.data.frame(res) %>%
  tibble::rownames_to_column("miRNA_ID") %>%
  dplyr::rename(
    Gene      = miRNA_ID,
    logFC     = log2FoldChange,
    adj.P.Val = padj
  ) %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::arrange(adj.P.Val)            # sort by significance for easy review

write.csv(final_dex, OUTPUT_FILE, row.names = FALSE)

cat(paste("\n--- SUCCESS ---\n"))
cat(paste("Output saved to:", OUTPUT_FILE, "\n"))
cat(paste("This file is ready to use as DEX_FILE_NAME in Script 02.\n"))
