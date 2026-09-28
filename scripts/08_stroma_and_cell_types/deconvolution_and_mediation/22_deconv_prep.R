#!/usr/bin/env Rscript
## 22_deconv_prep.R -- build the matrices every deconvolution method needs.
## Input : data/brca_gene_expr.rds  (Xena TCGA-BRCA HiSeqV2, log2(norm_count+1), 20530 x 1211)
## Output: cache/v3/deconv/bulk_log2.rds     log2 scale, primary tumours only
##         cache/v3/deconv/bulk_tpm.rds      linear, columns rescaled to 1e6
##         cache/v3/deconv/bulk_tpm.txt      tab file for command-line tools
## The linear matrix is an APPROXIMATION to TPM: Xena HiSeqV2 is RSEM upper-quartile
## normalised counts, not TPM. Un-logging and rescaling to 1e6 gives a per-sample
## compositional vector on the same scale as TPM but without effective-length
## correction. This is the standard practice for running EPIC/quanTIseq/CIBERSORT on
## TCGA Xena data and is recorded here so the manuscript can state it.
suppressPackageStartupMessages({library(matrixStats)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/22_deconv_prep.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

g  <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
ph <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
logf("raw matrix ", nrow(g), " x ", ncol(g), " range ", paste(range(g), collapse = " .. "))
stopifnot(!any(duplicated(rownames(g))))               # RULE 4 guard
logf("duplicated rownames: 0  (asserted)")
num_rn <- grepl("^[0-9]+$", rownames(g))
logf("purely numeric rownames (unmapped Entrez, kept as-is): ", sum(num_rn))

tum <- ph$sample[ph$sample_type == "Primary Tumor"]
ss  <- sort(intersect(colnames(g), tum))
logf("primary tumours with gene assay: n = ", length(ss))
G <- g[, ss, drop = FALSE]

## drop genes that are zero in every tumour -- they carry no information and some
## deconvolution packages choke on all-zero rows
nz <- rowSums(G) > 0
logf("genes with all-zero expression across tumours, dropped: ", sum(!nz))
G  <- G[nz, , drop = FALSE]
logf("log2 matrix kept: ", nrow(G), " x ", ncol(G))

TPM <- 2^G - 1
TPM[TPM < 0] <- 0
cs  <- colSums(TPM)
logf("pre-rescale column sums: median ", round(median(cs)), " min ", round(min(cs)), " max ", round(max(cs)))
TPM <- sweep(TPM, 2, cs, "/") * 1e6
logf("post-rescale column sums all == 1e6: ", all(abs(colSums(TPM) - 1e6) < 1e-6))
logf("TPM COL1A1 median = ", round(median(TPM["COL1A1", ]), 1),
     " ; ACTB median = ", round(median(TPM["ACTB", ]), 1),
     " ; PTPRC median = ", round(median(TPM["PTPRC", ]), 1))

saveRDS(G,   file.path(BASE, "cache/v3/deconv/bulk_log2.rds"))
saveRDS(TPM, file.path(BASE, "cache/v3/deconv/bulk_tpm.rds"))
write.table(data.frame(Gene = rownames(TPM), TPM, check.names = FALSE),
            file.path(BASE, "cache/v3/deconv/bulk_tpm.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = rownames(G), G, check.names = FALSE),
            file.path(BASE, "cache/v3/deconv/bulk_log2.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines(ss, file.path(BASE, "cache/v3/deconv/samples_primary_tumour.txt"))
logf("WROTE cache/v3/deconv/{bulk_log2.rds,bulk_tpm.rds,bulk_tpm.txt,bulk_log2.txt,samples_primary_tumour.txt}")
logf("DONE 22")
