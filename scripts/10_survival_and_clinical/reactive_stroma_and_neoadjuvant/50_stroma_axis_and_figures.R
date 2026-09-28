#!/usr/bin/env Rscript
## 50_stroma_axis_and_figures.R
##  (a) how much of each FFL module belongs to the UNION of published stromal/ECM programmes
##  (b) a single "reactive-stroma axis" (PC1 of eight independent stromal scores) and the
##      fraction of each module score's variance that this axis explains
##  (c) figures: forest plot of pCR odds ratios; module score vs Farmer score
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v5/50_stroma_axis.log")
logf <- function(...) { m <- paste0(...); cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE) }
cat("", file = LOG); set.seed(8)
sets <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
sc   <- readRDS(file.path(BASE, "cache/v5/farmer_tcga_scores.rds"))
gex  <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
X <- gex[, sc$samples, drop = FALSE]; X <- X[rowSums(is.na(X)) == 0 & rowSds(X) > 0, ]
UNIV <- rownames(X)
STROMAL <- c("FARMER_STROMAL","WEST_DTF_FIBROMATOSIS","WINSLOW_STROMAL_SIG1","FINAK_SDPP",
             "FINAK_REF_STROMA_I","FARMER2005_STROMAL_CLUSTER4","ESTIMATE_STROMAL_noCOL",
             "CAF_scRNA_50","CAF_MARKER_PANEL","NABA_CORE_MATRISOME","TRIULZI_ECM","HALLMARK_EMT")
U <- intersect(unique(unlist(sets[STROMAL])), UNIV)
logf("union of the 12 published stromal / ECM / EMT programmes: ", length(U), " genes")
MODS <- c("FFL_3node","FFL_4node","FFL_5node","FFL_6node_all","FFL_higher_order_only","MS_exemplar_module")
rows <- list(); k <- 0L
for (m in MODS) {
  A <- intersect(sets[[m]], UNIV); int <- intersect(A, U)
  p <- phyper(length(int)-1L, length(U), length(UNIV)-length(U), length(A), lower.tail = FALSE)
  k <- k + 1L
  rows[[k]] <- data.table(module = m, n_module = length(A), n_union = length(U),
    n_in_union = length(int), frac_of_module_stromal = length(int)/length(A),
    expected = length(A)*length(U)/length(UNIV),
    fold = length(int)/(length(A)*length(U)/length(UNIV)), p_hypergeometric = p,
    genes = paste(sort(int), collapse=";"))
}
UN <- rbindlist(rows)
fwrite(UN, file.path(BASE, "results/v5/farmer_module_union_stromal.csv"))
logf("\n=== fraction of each FFL module belonging to a published stromal/ECM/EMT programme ===")
for (i in seq_len(nrow(UN))) { r <- UN[i]
  logf(sprintf("  %-24s %3d/%3d = %.1f%% (expected %.1f%%, fold %.1f, hyperg p=%.3g)",
    r$module, r$n_in_union, r$n_module, 100*r$frac_of_module_stromal,
    100*r$expected/r$n_module, r$fold, r$p_hypergeometric)) }

## ---------------------------------------------- (b) the reactive-stroma axis (PC1) ------
AX <- c("FARMER_STROMAL","WEST_DTF_FIBROMATOSIS","WINSLOW_STROMAL_SIG1","ESTIMATE_STROMAL_noCOL",
        "CAF_scRNA_50","CAF_MARKER_PANEL","NABA_CORE_MATRISOME","TRIULZI_ECM")
S <- t(scale(t(sc$ssgsea[AX, , drop = FALSE])))
pc <- prcomp(t(S), center = TRUE, scale. = TRUE)
axis1 <- pc$x[, 1]
if (cor(axis1, sc$ssgsea["FARMER_STROMAL", ]) < 0) axis1 <- -axis1
ve <- summary(pc)$importance[2, 1]
logf("\n=== reactive-stroma axis = PC1 of ", length(AX), " independent stromal scores; ",
     "variance explained = ", round(100*ve, 1), "% ; rho with Farmer = ",
     round(cor(axis1, sc$ssgsea["FARMER_STROMAL", ], method="spearman"), 3))
ax_rows <- list(); k2 <- 0L
for (nm in rownames(sc$ssgsea)) {
  y <- as.numeric(sc$ssgsea[nm, ])
  h <- cor.test(y, axis1, method = "spearman", exact = FALSE)
  k2 <- k2 + 1L
  ax_rows[[k2]] <- data.table(score = nm, rho_vs_stroma_axis = unname(h$estimate),
    p = h$p.value, r2_pearson = cor(y, axis1)^2)
}
AXT <- rbindlist(ax_rows); setorder(AXT, -r2_pearson)
AXT[, variance_explained_pct := round(100*r2_pearson, 1)]
fwrite(AXT, file.path(BASE, "results/v5/farmer_reactive_stroma_axis.csv"))
logf("\n=== variance of each score explained by the reactive-stroma axis ===")
for (i in seq_len(nrow(AXT))) logf(sprintf("  %-28s rho=%+.3f  R2=%.3f  (%.1f%% of variance)",
  AXT$score[i], AXT$rho_vs_stroma_axis[i], AXT$r2_pearson[i], AXT$variance_explained_pct[i]))
saveRDS(axis1, file.path(BASE, "cache/v5/farmer_stroma_axis_tcga.rds"))

## ------------------------------------------------------------------- (c) figures --------
dir.create(file.path(BASE, "figures/v5"), showWarnings = FALSE, recursive = TRUE)
res  <- fread(file.path(BASE, "results/v5/farmer_neoadjuvant_pcr.csv"))
meta <- fread(file.path(BASE, "results/v5/farmer_neoadjuvant_meta.csv"))
pick <- c("FARMER_STROMAL","MS_exemplar_module","COLLAGEN_PAIR","WINSLOW_STROMAL_SIG1",
          "WEST_DTF_FIBROMATOSIS","NABA_CORE_MATRISOME","FFL_3node","FFL_higher_order_only",
          "FFL_6node_all","CHANG_WOUND_UP_VANTVEER")
pdf(file.path(BASE, "figures/v5/farmer_neoadjuvant_forest.pdf"), width = 9.5, height = 13.5)
DISP <- c(FARMER_STROMAL="Farmer stroma-related (50)",
          MS_exemplar_module="exemplar 6-node module",
          COLLAGEN_PAIR="COL1A1 + COL3A1",
          WINSLOW_STROMAL_SIG1="Winslow stromal sig. 1",
          WEST_DTF_FIBROMATOSIS="West/Beck DTF fibromatosis",
          NABA_CORE_MATRISOME="NABA core matrisome",
          FFL_3node="FFL 3-node module",
          FFL_higher_order_only="FFL higher-order-only",
          FFL_6node_all="FFL 6-node module (all)",
          CHANG_WOUND_UP_VANTVEER="Chang wound-activated (up)")
op <- par(mar = c(4.5, 14.5, 3, 2))
lab <- c(); est <- c(); lo <- c(); hi <- c(); col <- c()
for (s in pick) {
  sub <- res[score == s]
  for (i in seq_len(nrow(sub))) { lab <- c(lab, paste0("   ", sub$cohort[i]))
    est <- c(est, sub$OR_per_SD[i]); lo <- c(lo, sub$OR_low[i]); hi <- c(hi, sub$OR_high[i])
    col <- c(col, "grey35") }
  m <- meta[score == s]
  lab <- c(lab, paste0(DISP[[s]], "  \u2014 summary")); est <- c(est, m$OR_fixed)
  lo <- c(lo, m$fixed_lo); hi <- c(hi, m$fixed_hi); col <- c(col, "#b2182b")
  lab <- c(lab, ""); est <- c(est, NA); lo <- c(lo, NA); hi <- c(hi, NA); col <- c(col, NA)
}
n <- length(lab); yy <- n:1
plot(NA, xlim = c(0.3, 3.2), ylim = c(0.5, n + 0.5), log = "x", yaxt = "n",
     xlab = "odds ratio of pathological complete response per +1 SD of score",
     ylab = "", main = "Stromal and FFL module scores vs neoadjuvant chemotherapy response\n(7 GEO cohorts, 1,383 evaluable patients; OR < 1 = resistance)")
abline(v = 1, col = "grey60", lty = 2)
segments(lo, yy, hi, yy, col = col, lwd = 2)
points(est, yy, pch = ifelse(col == "#b2182b", 18, 16),
       cex = ifelse(col == "#b2182b", 1.6, 0.9), col = col)
axis(2, at = yy, labels = lab, las = 1, cex.axis = 0.6, tick = FALSE, gap.axis = -1)
par(op); dev.off()
logf("WROTE figures/v5/farmer_neoadjuvant_forest.pdf")

png(file.path(BASE, "figures/v5/farmer_vs_module_scatter.png"), width = 1500, height = 520, res = 130)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
f <- as.numeric(sc$ssgsea["FARMER_STROMAL", ])
for (m in c("FFL_3node","FFL_higher_order_only","MS_exemplar_module")) {
  y <- as.numeric(sc$ssgsea[m, ])
  plot(f, y, pch = 16, cex = 0.4, col = "#2c7fb877",
       xlab = "Farmer stroma-related score (ssGSEA)", ylab = paste(m, "score"),
       main = sprintf("%s\nSpearman rho = %.2f", m, cor(f, y, method = "spearman")))
  abline(lm(y ~ f), col = "#b2182b", lwd = 2)
}
par(op); dev.off()
logf("WROTE figures/v5/farmer_vs_module_scatter.png")
logf("DONE 50")
