#!/usr/bin/env Rscript
# Differential expression for one dataset of analyses/perturbation_tests (SETTINGS.md, "Differential expression").
# usage: Rscript _de.R <spec.json> <out.tsv>
# spec: {"kind": "counts" | "log_matrix" | "linear_matrix",
#        "matrix": path (TSV; first column gene or probe id; one column per sample),
#        "gene_map": optional TSV path (id -> symbol, columns id, symbol),
#        "samples": {"<column name>": "pert" | "ctrl", ...},
#        "tpm": optional TSV path (NCBI TPM table, same ids as matrix),
#        "extra_pseudo": optional list of two symbols to average into a pseudo-gene row (P1b)}
# counts: edgeR filterByExpr (defaults) + TMM + voom; matrices: log2 (if linear or max > 50) + limma.
# Probe/id -> gene: the id with the highest mean expression across all samples of the contrast.
suppressMessages({library(limma); library(edgeR); library(jsonlite)})
a <- commandArgs(TRUE); sp <- fromJSON(a[1]); out <- a[2]
M <- read.delim(sp$matrix, check.names = FALSE, stringsAsFactors = FALSE)
ids <- as.character(M[[1]]); M <- as.matrix(M[, -1, drop = FALSE]); rownames(M) <- ids
grp <- unlist(sp$samples); cols <- names(grp)
stopifnot(all(cols %in% colnames(M)))
M <- M[, cols, drop = FALSE]; storage.mode(M) <- "double"
g <- factor(ifelse(grp == "pert", "pert", "ctrl"), levels = c("ctrl", "pert"))
design <- model.matrix(~ g)
sym <- ids
if (!is.null(sp$gene_map)) {
  gm <- read.delim(sp$gene_map, stringsAsFactors = FALSE, colClasses = "character")
  sym <- gm$symbol[match(ids, gm$id)]
}
names(sym) <- ids
if (sp$kind == "counts") {
  keep_nonzero <- !is.na(sym) & sym != ""
  M <- M[keep_nonzero, , drop = FALSE]
  y <- DGEList(M); keep <- filterByExpr(y, group = g)
  if (!is.null(sp$force_keep)) {            # sensitivity option: keep named genes that filterByExpr would drop (counts > 0)
    fk <- names(sym)[sym %in% sp$force_keep]; fk <- intersect(fk, rownames(y))
    keep[fk] <- keep[fk] | rowSums(y$counts[fk, , drop = FALSE]) > 0
  }
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  v <- voom(y, design)
  E <- v$E; W <- v$weights
} else {
  X <- M[!is.na(sym[rownames(M)]) & sym[rownames(M)] != "", , drop = FALSE]
  if (sp$kind == "linear_matrix" || max(X, na.rm = TRUE) > 50) X <- log2(pmax(X, 0) + 1)
  X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
  E <- X; W <- NULL
}
# one row per gene symbol: the id with the highest mean expression
s <- sym[rownames(E)]; mu <- rowMeans(E)
o <- order(s, -mu); first <- o[!duplicated(s[o])]
E <- E[first, , drop = FALSE]; if (!is.null(W)) W <- W[first, , drop = FALSE]
rownames(E) <- s[first]; if (!is.null(W)) rownames(W) <- s[first]
ctrl_mean <- rowMeans(E[, g == "ctrl", drop = FALSE])
expressed <- if (sp$kind == "counts") setNames(rep(TRUE, nrow(E)), rownames(E)) else ctrl_mean >= median(ctrl_mean)
if (!is.null(sp$extra_pseudo)) {
  pg <- sp$extra_pseudo
  if (all(pg %in% rownames(E))) {
    E <- rbind(E, PSEUDO = colMeans(E[pg, , drop = FALSE]))
    if (!is.null(W)) W <- rbind(W, PSEUDO = colMeans(W[pg, , drop = FALSE]))
    ctrl_mean <- c(ctrl_mean, PSEUDO = NA); expressed <- c(expressed, PSEUDO = NA)
    rownames(E)[nrow(E)] <- paste(pg, collapse = "_")
    if (!is.null(W)) rownames(W)[nrow(W)] <- paste(pg, collapse = "_")
    names(ctrl_mean)[length(ctrl_mean)] <- names(expressed)[length(expressed)] <- paste(pg, collapse = "_")
  }
}
fit <- if (is.null(W)) lmFit(E, design) else lmFit(E, design, weights = W)
fit <- eBayes(fit, trend = is.null(W))
tt <- topTable(fit, coef = "gpert", number = Inf, sort.by = "none")
res <- data.frame(gene = rownames(tt), logFC = tt$logFC, SE = tt$logFC / tt$t, t = tt$t, P = tt$P.Value, adjP = tt$adj.P.Val,
                  AveExpr = tt$AveExpr, ctrl_mean_log2 = ctrl_mean[rownames(tt)], expressed = expressed[rownames(tt)],
                  n_pert = sum(g == "pert"), n_ctrl = sum(g == "ctrl"), stringsAsFactors = FALSE)
if (!is.null(sp$tpm)) {
  T <- read.delim(sp$tpm, check.names = FALSE, stringsAsFactors = FALSE)
  tid <- as.character(T[[1]]); T <- as.matrix(T[, cols[g == "ctrl"], drop = FALSE]); rownames(T) <- tid
  tsym <- sym[tid]; tm <- rowMeans(T)
  agg <- tapply(tm, tsym, max)
  res$ctrl_TPM <- as.numeric(agg[res$gene])
}
res$ctrl_pct <- 100 * rank(res$ctrl_mean_log2, na.last = "keep") / sum(!is.na(res$ctrl_mean_log2))
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("genes", nrow(res), "pert", sum(g == "pert"), "ctrl", sum(g == "ctrl"), "\n")
