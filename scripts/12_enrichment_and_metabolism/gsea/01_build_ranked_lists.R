#!/usr/bin/env Rscript
# v4/01 -- Build the three ranked lists for GSEA
# Outputs: results/v4/rank_*.csv  and  results/v4/ranklist_provenance.csv
suppressPackageStartupMessages({
  library(limma); library(data.table); library(matrixStats)
})
options(stringsAsFactors = FALSE)
ROOT <- "/path/to/revision"
RES  <- file.path(ROOT, "results", "v4")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
set.seed(20260908)

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---------------------------------------------------------------- load
expr  <- readRDS(file.path(ROOT, "data", "brca_gene_expr.rds"))
mir   <- readRDS(file.path(ROOT, "data", "brca_mirna_expr.rds"))
pheno <- readRDS(file.path(ROOT, "data", "brca_pheno.rds"))

## ---- MANDATORY ASSERTIONS (limma silent-failure guard, rule #4) -------
stopifnot(is.matrix(expr))
if (any(duplicated(rownames(expr))))
  stop("FATAL: duplicated rownames in expr -> limma topTable would return integer indices")
if (is.null(rownames(expr))) stop("FATAL: expr has no rownames")
frac_symbol <- mean(grepl("^[A-Za-z][A-Za-z0-9._-]*$", rownames(expr)))
msg("expr ", nrow(expr), " x ", ncol(expr),
    " | duplicated rownames = 0 | fraction of rownames that look like symbols = ",
    signif(frac_symbol, 4))
if (frac_symbol < 0.99) stop("FATAL: rownames do not look like gene symbols")

tum <- pheno$sample[pheno$sample_type == "Primary Tumor"]
nor <- pheno$sample[pheno$sample_type == "Solid Tissue Normal"]
tum <- intersect(tum, colnames(expr)); nor <- intersect(nor, colnames(expr))
msg("tumour n=", length(tum), " normal n=", length(nor))

## expression filter: non-zero variance across the 1,211 DE samples.
## This exactly reproduces the manuscript's tested universe (n = 20,250).
keep <- matrixStats::rowVars(expr[, c(tum, nor), drop = FALSE]) > 0
msg("genes passing non-zero-variance filter: ", sum(keep), " of ", nrow(expr))
E <- expr[keep, c(tum, nor), drop = FALSE]
grp <- factor(c(rep("Tumor", length(tum)), rep("Normal", length(nor))),
              levels = c("Normal", "Tumor"))

## ================================================================ (1) T vs N
des <- model.matrix(~ grp)
fit <- eBayes(lmFit(E, des))
tt  <- topTable(fit, coef = "grpTumor", number = Inf, sort.by = "none")
tt$feature <- rownames(E)
## guard: topTable must not have handed back integers
if (all(grepl("^[0-9]+$", rownames(tt)))) stop("FATAL: topTable returned integer indices")
stopifnot(identical(rownames(tt), rownames(E)))
r1 <- data.table(feature = tt$feature, stat = tt$t, logFC = tt$logFC,
                 AveExpr = tt$AveExpr, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
setorder(r1, -stat)
fwrite(r1, file.path(RES, "rank_tumour_vs_normal.csv"))
msg("(1) tumour-vs-normal: ", nrow(r1), " genes; t range ",
    paste(signif(range(r1$stat), 4), collapse = " .. "))

## sanity cross-check against the v1 DE table
old <- tryCatch(fread(file.path(ROOT, "results", "BRCA_DEX_genes.csv")), error = function(e) NULL)
if (!is.null(old)) {
  mrg <- merge(r1[, .(feature, stat)], old[, .(feature, t)], by = "feature")
  msg("    cross-check vs results/BRCA_DEX_genes.csv: n=", nrow(mrg),
      " Pearson r(t,t) = ", signif(cor(mrg$stat, mrg$t), 6))
}

## ================================================================ corr lists
com  <- intersect(colnames(expr), colnames(mir))
ctum <- intersect(com, pheno$sample[pheno$sample_type == "Primary Tumor"])
msg("primary tumours with paired gene+miRNA: ", length(ctum))

Ec <- expr[keep, ctum, drop = FALSE]
## drop zero-variance genes for correlation
vv <- rowVars(Ec); Ec <- Ec[vv > 0, , drop = FALSE]
msg("genes with non-zero variance in the ", length(ctum), " paired tumours: ", nrow(Ec))

spearman_rank <- function(mat, y) {
  ## rank-transform then Pearson == Spearman; ties handled by average ranks
  Rm <- t(apply(mat, 1, rank))
  ry <- rank(y)
  Rm <- Rm - rowMeans(Rm)
  sm <- sqrt(rowSums(Rm^2))
  ry <- ry - mean(ry); sy <- sqrt(sum(ry^2))
  as.numeric((Rm %*% ry) / (sm * sy))
}

make_corr_rank <- function(mirna_id, tag) {
  stopifnot(mirna_id %in% rownames(mir))
  y   <- as.numeric(mir[mirna_id, ctum])
  rho <- spearman_rank(Ec, y)
  n   <- length(ctum)
  tstat <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
  p   <- 2 * pt(-abs(tstat), df = n - 2)
  dt  <- data.table(feature = rownames(Ec), rho = rho, stat = rho,
                    t_stat = tstat, P.Value = p, fdr = p.adjust(p, "BH"),
                    n_samples = n, mirna = mirna_id)
  setorder(dt, -rho)
  fwrite(dt, file.path(RES, paste0("rank_", tag, ".csv")))
  msg("    ", tag, " (", mirna_id, "): rho range ",
      paste(signif(range(rho), 4), collapse = " .. "),
      " | FDR<0.05: ", sum(dt$fdr < 0.05),
      " (neg ", sum(dt$fdr < 0.05 & dt$rho < 0), " / pos ",
      sum(dt$fdr < 0.05 & dt$rho > 0), ")")
  invisible(dt)
}
msg("(2) Spearman correlation ranking")
c29  <- make_corr_rank("hsa-miR-29a-3p",  "mir29a_corr")

## ================================================================ (3) signed influence
sif <- file.path(ROOT, "results", "signed_influence_scores.csv")
if (file.exists(sif)) {
  s <- fread(sif)
  msg("(3) signed_influence_scores.csv found: ", nrow(s), " nodes; cols: ",
      paste(names(s), collapse = ", "))
  msg("    node types: ", paste(names(table(s$type)), table(s$type), sep = "=", collapse = " "))
  ## keep protein-coding nodes (Gene + TF) -- gene sets are gene-symbol based
  sp <- s[type %in% c("Gene", "TF")]
  out <- data.table(feature = sp$name, stat = sp$influence_curated,
                    influence_as_published = sp$influence_as_published,
                    influence_curated = sp$influence_curated,
                    node_type = sp$type, is_seed = sp$is_seed)
  out <- out[!is.na(stat)]
  setorder(out, -stat)
  fwrite(out, file.path(RES, "rank_ffl_signed_influence.csv"))
  msg("    kept ", nrow(out), " protein-coding network nodes; stat range ",
      paste(signif(range(out$stat), 4), collapse = " .. "))
  ## also keep the as-published variant for a sensitivity check
  out2 <- copy(out); out2[, stat := influence_as_published]
  out2 <- out2[!is.na(stat)]; setorder(out2, -stat)
  fwrite(out2, file.path(RES, "rank_ffl_signed_influence_aspublished.csv"))
} else {
  stop("signed_influence_scores.csv missing -- fallback not implemented")
}

## ================================================================ provenance
prov <- data.table(
  ranked_list = c("tumour_vs_normal", "mir29a_corr",
                  "ffl_signed_influence"),
  statistic = c("limma moderated t (Tumor vs Normal)",
                "Spearman rho vs hsa-miR-29a-3p",
                "signed influence (curated network)"),
  n_features = c(nrow(r1), nrow(c29),
                 nrow(fread(file.path(RES, "rank_ffl_signed_influence.csv")))),
  n_samples = c(length(tum) + length(nor), length(ctum),
                NA_integer_),
  detail = c(paste0(length(tum), " tumour vs ", length(nor), " normal"),
             paste0(length(ctum), " paired primary tumours"),
             "587 network nodes -> protein-coding subset")
)
fwrite(prov, file.path(RES, "ranklist_provenance.csv"))
print(prov)
msg("DONE 01")
