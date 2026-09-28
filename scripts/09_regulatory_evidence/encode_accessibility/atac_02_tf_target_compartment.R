#!/usr/bin/env Rscript
## Which TF->target correlations in TCGA-BRCA survive adjustment for stromal content?
## Used to (i) place COL1A1/COL3A1 on that spectrum and (ii) pick positive-control
## targets for the ATAC motif analysis.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
setwd(BASE)

X   <- readRDS("data/brca_gene_expr_symbol.rds")          # genes x samples, log2
sc  <- fread("results/v5/farmer_tcga_persample_scores.csv")
tft <- fread("data/layer_TF_target.tsv")

TFS <- c("ETS1","NFKB1","RELA","SP1")
common <- intersect(colnames(X), sc$sample)
cat("RNA samples:", ncol(X), " scored samples:", nrow(sc), " common:", length(common), "\n")
X  <- X[, common]; sc <- sc[match(common, sc$sample)]

# stromal covariates: collagen-free by construction (avoids circularity with COL1A1/3A1)
S1 <- sc$meanZ_ESTIMATE_STROMAL_noCOL
S2 <- sc$meanZ_CAF_scRNA_50
stopifnot(!any(is.na(S1)), !any(is.na(S2)))

pcor_spearman <- function(x, y, z) {           # partial Spearman via rank residuals
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  ex <- residuals(lm(rx ~ rz)); ey <- residuals(lm(ry ~ rz))
  n  <- length(rx)
  r  <- cor(ex, ey)
  df <- n - 3
  tt <- r * sqrt(df / (1 - r^2))
  list(rho = r, p = 2 * pt(-abs(tt), df))
}

res <- list()
for (tf in TFS) {
  if (!tf %in% rownames(X)) { cat("TF missing from RNA:", tf, "\n"); next }
  tg <- unique(tft[source == tf & target %in% rownames(X), target])
  tg <- setdiff(tg, tf)
  cat(tf, ": ", length(tg), " curated targets present in RNA\n", sep = "")
  xv <- as.numeric(X[tf, ])
  for (g in tg) {
    yv <- as.numeric(X[g, ])
    if (sd(yv) == 0) next
    ct  <- cor.test(xv, yv, method = "spearman", exact = FALSE)
    p1  <- pcor_spearman(xv, yv, S1)
    p2  <- pcor_spearman(xv, yv, S2)
    res[[length(res) + 1]] <- data.table(
      TF = tf, target = g, n = length(xv),
      rho_raw = unname(ct$estimate), p_raw = ct$p.value,
      rho_adj_estimateNoCol = p1$rho, p_adj_estimateNoCol = p1$p,
      rho_adj_cafScrna = p2$rho, p_adj_cafScrna = p2$p,
      rho_gene_stroma = cor(rank(yv), rank(S1)),
      rho_tf_stroma   = cor(rank(xv), rank(S1)))
  }
}
R <- rbindlist(res)
R[, attenuation := 1 - abs(rho_adj_estimateNoCol) / abs(rho_raw)]
R[, fdr_raw := p.adjust(p_raw, "BH")]
R[, fdr_adj := p.adjust(p_adj_estimateNoCol, "BH")]
R[, sign_consistent := sign(rho_raw) == sign(rho_adj_estimateNoCol)]
R[, class := fifelse(fdr_raw < 0.05 & abs(rho_raw) >= 0.30 & sign_consistent &
                     fdr_adj < 0.05 & abs(rho_adj_estimateNoCol) >= 0.30 &
                     attenuation <= 0.25, "compartment_robust",
             fifelse(fdr_raw < 0.05 & abs(rho_raw) >= 0.30 &
                     (abs(rho_adj_estimateNoCol) < 0.15 | attenuation >= 0.50 | !sign_consistent),
                     "compartment_explained", "intermediate"))]
## positive-control set for the ATAC motif analysis: sign-consistent, minimally attenuated
R[, poscontrol := fdr_adj < 0.05 & sign_consistent & attenuation <= 0.25 &
                  abs(rho_adj_estimateNoCol) >= 0.20 & abs(rho_raw) >= 0.20]
setorder(R, TF, -rho_raw)
fwrite(R, "results/v6/atac_tf_target_compartment_adjusted.csv")

cat("\n== class counts ==\n"); print(R[, .N, by = .(TF, class)][order(TF, class)])
cat("\n== collagen edges ==\n")
print(R[target %in% c("COL1A1","COL3A1"),
        .(TF, target, rho_raw, rho_adj_estimateNoCol, rho_adj_cafScrna, attenuation, class)])
cat("\n== positive controls per TF ==\n"); print(R[poscontrol == TRUE, .N, by = TF])
print(R[poscontrol == TRUE][order(TF, -abs(rho_adj_estimateNoCol)), paste(target, collapse=", "), by = TF])
cat("\n== compartment-robust, top per TF ==\n")
print(R[class == "compartment_robust"][order(TF, -abs(rho_adj_estimateNoCol))][,
        head(.SD, 8), by = TF, .SDcols = c("target","rho_raw","rho_adj_estimateNoCol","attenuation")])
cat("\nDONE\n")
