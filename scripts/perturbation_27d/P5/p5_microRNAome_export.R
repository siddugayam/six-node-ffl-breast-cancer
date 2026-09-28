#!/usr/bin/env Rscript
# P5: export from the Bioconductor microRNAome package (McCall et al. 2017) the counts of the four miRNAs and each
# sample's total miRNA count (column sum over all miRNAs), with CellType and Class, to the 27d cache.
suppressMessages({library(microRNAome); library(SummarizedExperiment)})
data(microRNAome)
X <- assay(microRNAome, "counts"); cd <- as.data.frame(colData(microRNAome))
mirs <- c("hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p", "hsa-miR-101-3p")
stopifnot(all(mirs %in% rownames(X)))
out <- data.frame(sample = colnames(X), CellType = cd$CellType, Class = cd$Class, total_miRNA_counts = colSums(X), t(X[mirs, , drop = FALSE]), check.names = FALSE)
a <- commandArgs(TRUE)
write.table(out, a[1], sep = "\t", quote = FALSE, row.names = FALSE)
cat(readLines(system.file("CITATION", package = "microRNAome")), sep = "\n")   # citation() fails on the package's CITATION file
cat("version", as.character(packageVersion("microRNAome")), "licence", packageDescription("microRNAome")$License, "\n")
