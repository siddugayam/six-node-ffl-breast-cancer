#!/usr/bin/env Rscript
# Agilent Feature Extraction files -> normalised matrices (SETTINGS.md: "normexp + quantile for Agilent").
# usage: Rscript _agilent.R single <out.tsv> <files...>       single-channel: normexp (offset 16), quantile; log2 E
#        Rscript _agilent.R two <out_prefix> <files...>       two-colour: normexp, loess within arrays, Aquantile
#                                                              -> <prefix>.M.tsv (log2 R/G) and <prefix>.A.tsv
# Control features (ControlType != 0) are dropped; replicate probes are averaged by ProbeName.
# Column names = the GSM prefix of each file name.
suppressMessages(library(limma))
a <- commandArgs(TRUE); mode <- a[1]; out <- a[2]; files <- a[-(1:2)]
gsm <- sub("^(GSM[0-9]+).*$", "\\1", basename(files))
if (mode == "single") {
  RG <- read.maimages(files, source = "agilent.median", green.only = TRUE, verbose = FALSE)
  bc <- backgroundCorrect(RG, method = "normexp", offset = 16, verbose = FALSE)
  n <- normalizeBetweenArrays(bc, method = "quantile")
  keep <- n$genes$ControlType == 0
  E <- avereps(n$E[keep, , drop = FALSE], ID = n$genes$ProbeName[keep]); colnames(E) <- gsm
  write.table(data.frame(id = rownames(E), E, check.names = FALSE), out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("probes", nrow(E), "arrays", ncol(E), "single-channel\n")
} else {
  RG <- read.maimages(files, source = "agilent.median", verbose = FALSE)
  bc <- backgroundCorrect(RG, method = "normexp", offset = 16, verbose = FALSE)
  MA <- normalizeWithinArrays(bc, method = "loess"); MA <- normalizeBetweenArrays(MA, method = "Aquantile")
  keep <- MA$genes$ControlType == 0
  M <- avereps(MA$M[keep, , drop = FALSE], ID = MA$genes$ProbeName[keep]); A <- avereps(MA$A[keep, , drop = FALSE], ID = MA$genes$ProbeName[keep])
  colnames(M) <- colnames(A) <- gsm
  write.table(data.frame(id = rownames(M), M, check.names = FALSE), paste0(out, ".M.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(id = rownames(A), A, check.names = FALSE), paste0(out, ".A.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("probes", nrow(M), "arrays", ncol(M), "two-colour\n")
}
