#!/usr/bin/env Rscript
# One-sample limma on two-colour log ratios (perturbed / control), for paired two-colour designs.
# usage: Rscript _de2.R <spec.json> <out.tsv>
# spec: {"M": path, "A": path or null, "gene_map": path, "orient": {"<GSM>": 1 | -1}, "extra_pseudo": [...]}
# Probe -> gene: the probe with the highest mean A (or highest mean |M| if A is absent).
suppressMessages({library(limma); library(jsonlite)})
a <- commandArgs(TRUE); sp <- fromJSON(a[1])
rd <- function(p) { x <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE); m <- as.matrix(x[, -1, drop = FALSE]); rownames(m) <- as.character(x[[1]]); storage.mode(m) <- "double"; m }
M <- rd(sp$M); o <- unlist(sp$orient); M <- sweep(M[, names(o), drop = FALSE], 2, o, "*")
A <- if (!is.null(sp$A)) rd(sp$A)[rownames(M), names(o), drop = FALSE] else NULL
gm <- read.delim(sp$gene_map, stringsAsFactors = FALSE, colClasses = "character"); sym <- gm$symbol[match(rownames(M), gm$id)]
ok <- !is.na(sym) & sym != "" & rowSums(!is.na(M)) >= 2          # limma handles missing values within a row
M <- M[ok, , drop = FALSE]; if (!is.null(A)) A <- A[ok, , drop = FALSE]; sym <- sym[ok]
score <- if (!is.null(A)) rowMeans(A, na.rm = TRUE) else rowMeans(abs(M), na.rm = TRUE)
ord <- order(sym, -score); first <- ord[!duplicated(sym[ord])]
M <- M[first, , drop = FALSE]; rownames(M) <- sym[first]
Aexp <- if (!is.null(A)) { x <- A[first, , drop = FALSE]; rownames(x) <- sym[first]; rowMeans(x, na.rm = TRUE) } else rep(NA, nrow(M))
names(Aexp) <- rownames(M)
if (!is.null(sp$extra_pseudo) && all(sp$extra_pseudo %in% rownames(M))) {
  M <- rbind(M, colMeans(M[sp$extra_pseudo, , drop = FALSE])); rownames(M)[nrow(M)] <- paste(sp$extra_pseudo, collapse = "_")
  Aexp <- c(Aexp, NA); names(Aexp)[length(Aexp)] <- rownames(M)[nrow(M)]
}
fit <- eBayes(lmFit(M, matrix(1, ncol(M), 1, dimnames = list(NULL, "int"))))
tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
pct <- 100 * rank(Aexp, na.last = "keep") / sum(!is.na(Aexp))
res <- data.frame(gene = rownames(tt), logFC = tt$logFC, SE = tt$logFC / tt$t, t = tt$t, P = tt$P.Value, adjP = tt$adj.P.Val, AveExpr = Aexp[rownames(tt)],
                  ctrl_mean_log2 = Aexp[rownames(tt)], expressed = if (all(is.na(Aexp))) TRUE else Aexp[rownames(tt)] >= median(Aexp, na.rm = TRUE),
                  n_pert = ncol(M), n_ctrl = ncol(M), ctrl_pct = pct[rownames(tt)], stringsAsFactors = FALSE)
write.table(res, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat("genes", nrow(res), "arrays", ncol(M), "\n")
