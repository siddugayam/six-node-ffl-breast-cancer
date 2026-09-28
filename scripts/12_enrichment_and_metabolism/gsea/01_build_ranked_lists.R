#!/usr/bin/env Rscript
# v4/01 -- Build the five ranked lists for GSEA
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
msg("(2)(3) Spearman correlation rankings")
c130 <- make_corr_rank("hsa-miR-130a-3p", "mir130a_corr")
c29  <- make_corr_rank("hsa-miR-29a-3p",  "mir29a_corr")

## sanity vs v3 file
v3 <- tryCatch(fread(file.path(ROOT, "results", "v3", "mir130a_tcga_target_correlations.csv")),
               error = function(e) NULL)
if (!is.null(v3) && "rho_3p" %in% names(v3)) {
  mm <- merge(c130[, .(gene = feature, rho_new = rho)], v3[, .(gene, rho_3p)], by = "gene")
  msg("    cross-check vs v3 mir130a_tcga_target_correlations: n=", nrow(mm),
      " Pearson r = ", signif(cor(mm$rho_new, mm$rho_3p), 6))
}

## ================================================================ (4) signed influence
sif <- file.path(ROOT, "results", "signed_influence_scores.csv")
if (file.exists(sif)) {
  s <- fread(sif)
  msg("(4) signed_influence_scores.csv found: ", nrow(s), " nodes; cols: ",
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

## ================================================================ (5) 130a high vs low tertile
y130 <- as.numeric(mir["hsa-miR-130a-3p", ctum])
qs   <- quantile(y130, c(1/3, 2/3))
grp2 <- ifelse(y130 <= qs[1], "low", ifelse(y130 >= qs[2], "high", "mid"))
msg("(5) miR-130a-3p tertiles: ", paste(names(table(grp2)), table(grp2), sep = "=", collapse = " "),
    " | cutpoints ", paste(signif(qs, 5), collapse = ", "))
sel  <- grp2 %in% c("low", "high")
Eh   <- expr[keep, ctum[sel], drop = FALSE]
Eh   <- Eh[rowVars(Eh) > 0, , drop = FALSE]
gh   <- factor(grp2[sel], levels = c("low", "high"))
d2   <- model.matrix(~ gh)
f2   <- eBayes(lmFit(Eh, d2))
t2   <- topTable(f2, coef = "ghhigh", number = Inf, sort.by = "none")
if (all(grepl("^[0-9]+$", rownames(t2)))) stop("FATAL: topTable returned integer indices (contrast 5)")
stopifnot(identical(rownames(t2), rownames(Eh)))
r5 <- data.table(feature = rownames(Eh), stat = t2$t, logFC = t2$logFC,
                 AveExpr = t2$AveExpr, P.Value = t2$P.Value, adj.P.Val = t2$adj.P.Val)
setorder(r5, -stat)
fwrite(r5, file.path(RES, "rank_mir130a_high_vs_low.csv"))
msg("    ", nrow(r5), " genes; FDR<0.05: ", sum(r5$adj.P.Val < 0.05),
    " (up in high ", sum(r5$adj.P.Val < 0.05 & r5$logFC > 0),
    " / down in high ", sum(r5$adj.P.Val < 0.05 & r5$logFC < 0), ")")
## also record the sample-level group assignment for the survival arm
fwrite(data.table(sample = ctum, mir130a_3p = y130, tertile = grp2),
       file.path(RES, "mir130a_tertile_assignment.csv"))

## ================================================================ provenance
prov <- data.table(
  ranked_list = c("tumour_vs_normal", "mir130a_corr", "mir29a_corr",
                  "ffl_signed_influence", "mir130a_high_vs_low"),
  statistic = c("limma moderated t (Tumor vs Normal)",
                "Spearman rho vs hsa-miR-130a-3p",
                "Spearman rho vs hsa-miR-29a-3p",
                "signed influence (curated network)",
                "limma moderated t (miR-130a-3p high vs low tertile)"),
  n_features = c(nrow(r1), nrow(c130), nrow(c29),
                 nrow(fread(file.path(RES, "rank_ffl_signed_influence.csv"))), nrow(r5)),
  n_samples = c(length(tum) + length(nor), length(ctum), length(ctum),
                NA_integer_, sum(sel)),
  detail = c(paste0(length(tum), " tumour vs ", length(nor), " normal"),
             paste0(length(ctum), " paired primary tumours"),
             paste0(length(ctum), " paired primary tumours"),
             "587 network nodes -> protein-coding subset",
             paste0(sum(grp2 == "high"), " high vs ", sum(grp2 == "low"), " low"))
)
fwrite(prov, file.path(RES, "ranklist_provenance.csv"))
print(prov)
msg("DONE 01")
