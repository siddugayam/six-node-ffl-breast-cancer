#!/usr/bin/env Rscript
## 42_farmer_scores_tcga.R
## Per-sample Farmer stromal score (and every comparison signature) in TCGA-BRCA primary
## tumours by ssGSEA (Barbie), GSVA (Gaussian kcdf) and the mean-of-z score that Farmer
## et al. actually used ("the average of the 50 genes"), then correlate with
##   (a) our CAF / stromal estimates,
##   (b) COL1A1, COL3A1 and the collagen module,
##   (c) the 3-node and higher-order FFL module scores,
##   (d) the 30 prioritised nodes.
suppressPackageStartupMessages({
  library(data.table); library(GSVA); library(matrixStats)
})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v5/42_farmer_scores_tcga.log")
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); set.seed(42)

sets  <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
gex   <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mex   <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
pheno <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
tum   <- pheno$sample[pheno$sample_type == "Primary Tumor"]
S     <- sort(intersect(colnames(gex), tum))
logf("TCGA-BRCA primary tumours with gene expression: n = ", length(S))
X <- gex[, S, drop = FALSE]
keep <- rowSums(is.na(X)) == 0 & rowSds(X) > 0
X <- X[keep, , drop = FALSE]
logf("genes used for scoring (no NA, non-zero variance): ", nrow(X), " of ", nrow(gex))
logf("expression matrix range: ", paste(round(range(X), 2), collapse = " .. "), " (log2 scale)")

## keep only gene-level (protein-coding) members
gsets <- lapply(sets, function(s) intersect(s, rownames(X)))
gsets <- gsets[sapply(gsets, length) >= 2]
logf("signatures scored: ", length(gsets))

## ------------------------------------------------------------------ three scorings ----
logf("running ssGSEA (Barbie) ...")
ss  <- gsva(ssgseaParam(X, gsets, minSize = 2, normalize = TRUE), verbose = FALSE)
logf("running GSVA (Gaussian kcdf) ...")
gv  <- gsva(gsvaParam(X, gsets, kcdf = "Gaussian", minSize = 2), verbose = FALSE)
Z   <- (X - rowMeans(X)) / rowSds(X)
mz  <- t(sapply(gsets, function(s) colMeans(Z[s, , drop = FALSE])))
logf("score matrices: ssGSEA ", paste(dim(ss), collapse = "x"),
     " | GSVA ", paste(dim(gv), collapse = "x"), " | meanZ ", paste(dim(mz), collapse = "x"))

## concordance of the three scorings for the Farmer signature
for (nm in c("FARMER_STROMAL", "FFL_higher_order_only", "FFL_3node")) {
  logf(sprintf("%-24s rho(ssGSEA,meanZ)=%.3f  rho(GSVA,meanZ)=%.3f  rho(ssGSEA,GSVA)=%.3f", nm,
       cor(ss[nm, ], mz[nm, ], method = "spearman"),
       cor(gv[nm, ], mz[nm, ], method = "spearman"),
       cor(ss[nm, ], gv[nm, ], method = "spearman")))
}
saveRDS(list(ssgsea = ss, gsva = gv, meanz = mz, samples = S, sets = gsets),
        file.path(BASE, "cache/v5/farmer_tcga_scores.rds"))
scdt <- data.table(sample = S)
for (nm in rownames(ss)) scdt[[paste0("ssGSEA_", nm)]] <- as.numeric(ss[nm, ])
for (nm in rownames(mz)) scdt[[paste0("meanZ_",  nm)]] <- as.numeric(mz[nm, ])
fwrite(scdt, file.path(BASE, "results/v5/farmer_tcga_persample_scores.csv"))
logf("WROTE results/v5/farmer_tcga_persample_scores.csv  (", nrow(scdt), " samples)")

## ------------------------------------------------------- (a) our CAF/stromal estimates --
own <- fread(file.path(BASE, "results/multiomics/tcga_stromal_scores.csv"))
setkey(own, sample)
own <- own[S]
stopifnot(identical(own$sample, S))
ownv <- c("ssgsea_stromal", "ssgsea_stromal_noCOL", "ESTIMATE_StromalScore", "ESTIMATE_ImmuneScore",
          "score_CAF_marker", "score_CAF_scRNA", "score_epithelial", "score_immune", "score_endothelial")
farmer_ss <- as.numeric(ss["FARMER_STROMAL", ]); farmer_mz <- as.numeric(mz["FARMER_STROMAL", ])
rows <- list(); k <- 0L
ct <- function(a, b) { h <- cor.test(a, b, method = "spearman", exact = FALSE)
                       list(rho = unname(h$estimate), p = h$p.value) }
for (v in ownv) {
  y <- own[[v]]
  for (sc in c("ssGSEA", "meanZ")) {
    x <- if (sc == "ssGSEA") farmer_ss else farmer_mz
    h <- ct(x, y); k <- k + 1L
    rows[[k]] <- data.table(block = "our_stromal_estimate", farmer_score = sc, partner = v,
                            n = sum(complete.cases(x, y)), rho = h$rho, p = h$p)
  }
}
## ------------------------------------------------ (b) collagens and the collagen module --
for (v in c("COL1A1", "COL3A1", "COL1A2", "COL5A2", "FN1", "DCN", "POSTN")) {
  if (!v %in% rownames(X)) next
  y <- as.numeric(X[v, ])
  for (sc in c("ssGSEA", "meanZ")) {
    x <- if (sc == "ssGSEA") farmer_ss else farmer_mz
    h <- ct(x, y); k <- k + 1L
    rows[[k]] <- data.table(block = "gene", farmer_score = sc, partner = v,
                            n = length(S), rho = h$rho, p = h$p)
  }
}
## ------------------------------------------- (c) FFL / other signature module scores ----
for (v in setdiff(rownames(ss), "FARMER_STROMAL")) {
  for (sc in c("ssGSEA", "meanZ")) {
    M <- if (sc == "ssGSEA") ss else mz
    x <- as.numeric(M["FARMER_STROMAL", ]); y <- as.numeric(M[v, ])
    h <- ct(x, y); k <- k + 1L
    rows[[k]] <- data.table(block = "signature_score", farmer_score = sc, partner = v,
                            n = length(S), rho = h$rho, p = h$p)
  }
}
res <- rbindlist(rows)
res[, q_BH := p.adjust(p, "BH"), by = .(block, farmer_score)]
setorder(res, block, farmer_score, -rho)
fwrite(res, file.path(BASE, "results/v5/farmer_correlations_tcga.csv"))
logf("WROTE results/v5/farmer_correlations_tcga.csv  rows=", nrow(res))
logf("\n=== Farmer stromal score (ssGSEA) vs our estimates and the FFL modules ===")
pr <- res[farmer_score == "ssGSEA"][order(-abs(rho))]
for (i in seq_len(nrow(pr)))
  logf(sprintf("  %-22s %-28s rho=%+.3f  p=%.3g", pr$block[i], pr$partner[i], pr$rho[i], pr$p[i]))

## --------------------------------------------------- (d) the 30 prioritised nodes -------
top30 <- fread(file.path(BASE, "results/v5/top10_per_class.csv"))[, .(type, name)]
Sb <- sort(intersect(S, colnames(mex)))
logf("\nsamples with gene + miRNA assay: ", length(Sb))
fs_g <- setNames(farmer_ss, S); fm_g <- setNames(farmer_mz, S)
nrows <- list(); k2 <- 0L
for (i in seq_len(nrow(top30))) {
  nm <- top30$name[i]; ty <- top30$type[i]
  if (ty == "miRNA") { if (!nm %in% rownames(mex)) next
    ss_ <- Sb; y <- as.numeric(mex[nm, ss_]) } else { if (!nm %in% rownames(X)) next
    ss_ <- S;  y <- as.numeric(X[nm, ss_]) }
  for (sc in c("ssGSEA", "meanZ")) {
    x <- if (sc == "ssGSEA") fs_g[ss_] else fm_g[ss_]
    h <- ct(as.numeric(x), y); k2 <- k2 + 1L
    nrows[[k2]] <- data.table(node = nm, type = ty, farmer_score = sc, n = length(ss_),
                              rho = h$rho, p = h$p)
  }
}
nres <- rbindlist(nrows)
nres[, q_BH := p.adjust(p, "BH"), by = farmer_score]
setorder(nres, farmer_score, -rho)
fwrite(nres, file.path(BASE, "results/v5/farmer_vs_prioritised_nodes.csv"))
logf("WROTE results/v5/farmer_vs_prioritised_nodes.csv rows=", nrow(nres))
logf("\n=== Farmer stromal score (ssGSEA) vs the 30 prioritised nodes ===")
pn <- nres[farmer_score == "ssGSEA"]
for (i in seq_len(nrow(pn)))
  logf(sprintf("  %-6s %-16s rho=%+.3f  p=%.3g  q=%.3g", pn$type[i], pn$node[i], pn$rho[i], pn$p[i], pn$q_BH[i]))
logf("DONE 42")
