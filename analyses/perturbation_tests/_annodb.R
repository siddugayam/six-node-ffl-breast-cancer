#!/usr/bin/env Rscript
# id -> SYMBOL from a Bioconductor annotation package (used when the GPL table has no gene symbols, e.g. GPL23159).
# usage: Rscript _annodb.R <package.db> <ids.txt> <out.tsv>
a <- commandArgs(TRUE); suppressMessages(library(a[1], character.only = TRUE))
ids <- readLines(a[2]); db <- get(a[1])
m <- suppressMessages(AnnotationDbi::select(db, keys = ids, columns = "SYMBOL", keytype = "PROBEID"))
m <- m[!duplicated(m$PROBEID), ]
write.table(data.frame(id = m$PROBEID, symbol = ifelse(is.na(m$SYMBOL), "", m$SYMBOL)), a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat("mapped", sum(!is.na(m$SYMBOL)), "of", length(ids), "\n")
