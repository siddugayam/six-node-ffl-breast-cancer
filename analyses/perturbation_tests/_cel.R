#!/usr/bin/env Rscript
# RMA of Affymetrix CEL files with oligo (SETTINGS.md: "RMA for Affymetrix CEL").  Gene/exon ST and HTA arrays are
# summarised to transcript clusters (target = "core"); 3' IVT arrays to probe sets.
# usage: Rscript _cel.R <out_matrix.tsv> <cel files...>   (column names = the GSM prefix of each file name)
suppressMessages(library(oligo))
a <- commandArgs(TRUE); out <- a[1]; cels <- a[-1]
raw <- read.celfiles(cels)
e <- if (is(raw, "GeneFeatureSet")) rma(raw, target = "core") else rma(raw)
m <- exprs(e); colnames(m) <- sub("^(GSM[0-9]+).*$", "\\1", basename(colnames(m)))
write.table(data.frame(id = rownames(m), m, check.names = FALSE), out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("probesets", nrow(m), "samples", ncol(m), class(raw)[1], annotation(raw), "\n")
