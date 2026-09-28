#!/usr/bin/env Rscript
## 32_tpm_length_sensitivity.R
## Xena TCGA-BRCA HiSeqV2 is log2(RSEM normalised count + 1): counts, not TPM, so the
## linear matrix used for EPIC / quanTIseq / CIBERSORT is not length-corrected. Here we
## build a properly length-corrected TPM (count / union-exon length, renormalised to
## 1e6) and re-run the TPM-dependent methods on it, to test whether the approximation
## changes any fibroblast estimate.
suppressPackageStartupMessages({library(matrixStats); library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/32_tpm_sensitivity.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")
TPM <- readRDS(file.path(D, "bulk_tpm.rds"))
L   <- fread(file.path(D, "gene_lengths_gencode_v36.tsv"))
len <- setNames(L$union_exon_length, L$gene)
gi  <- intersect(rownames(TPM), names(len))
logf("genes with a GENCODE v36 union-exon length: ", length(gi), " of ", nrow(TPM))
A <- TPM[gi, ]
B <- A / len[gi]
B <- sweep(B, 2, colSums(B), "/") * 1e6
logf("length-corrected TPM built; COL1A1 median ", round(median(B["COL1A1", ]), 1),
     " (uncorrected ", round(median(A["COL1A1", ]), 1), ")")
saveRDS(B, file.path(D, "bulk_tpm_lengthcorrected.rds"))

res <- list()
try({ library(EPIC)
  e1 <- EPIC(bulk = A, reference = "TRef", withOtherCells = TRUE)$cellFractions
  e2 <- EPIC(bulk = B, reference = "TRef", withOtherCells = TRUE)$cellFractions
  res$EPIC_CAFs <- c(rho = cor(e1[, "CAFs"], e2[, "CAFs"], method = "spearman"),
                     med_uncorr = median(e1[, "CAFs"]), med_corr = median(e2[, "CAFs"]))
  res$EPIC_Endothelial <- c(rho = cor(e1[, "Endothelial"], e2[, "Endothelial"], method = "spearman"),
                            med_uncorr = median(e1[, "Endothelial"]), med_corr = median(e2[, "Endothelial"]))
  saveRDS(e2, file.path(D, "out_lengthcorrected_EPIC.rds")) }, silent = FALSE)
try({ library(quantiseqr)
  q1 <- run_quantiseq(A, signature_matrix = "TIL10", is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE)
  q2 <- run_quantiseq(B, signature_matrix = "TIL10", is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE)
  for (cc in setdiff(names(q1), "Sample"))
    res[[paste0("quanTIseq_", cc)]] <- c(rho = cor(q1[[cc]], q2[[cc]], method = "spearman"),
                                         med_uncorr = median(q1[[cc]]), med_corr = median(q2[[cc]])) }, silent = FALSE)
try({ source(file.path(BASE, "scripts/v3/23_cibersort_impl.R"))
  SIG <- as.matrix(read.delim(file.path(D, "wu2021_signature_matrix.txt"), row.names = 1, check.names = FALSE))
  r1 <- run_cibersort(SIG, A, QN = FALSE, cores = 4)$fractions
  r2 <- run_cibersort(SIG, B, QN = FALSE, cores = 4)$fractions
  for (cc in colnames(r1))
    res[[paste0("CBSX_Wu_", cc)]] <- c(rho = cor(r1[, cc], r2[, cc], method = "spearman"),
                                       med_uncorr = median(r1[, cc]), med_corr = median(r2[, cc]))
  saveRDS(r2, file.path(D, "out_lengthcorrected_CBSX_Wu.rds")) }, silent = FALSE)

R <- data.table(feature = names(res), do.call(rbind, res))
fwrite(R, file.path(BASE, "results/v3/deconv_tpm_length_sensitivity.csv"))
logf("\n=== effect of length-correcting the TPM approximation ===")
for (i in seq_len(nrow(R))) logf(sprintf("%-28s spearman(uncorrected, corrected) = %+.4f   median %.4f -> %.4f",
  R$feature[i], R$rho[i], R$med_uncorr[i], R$med_corr[i]))
logf("DONE 32")
