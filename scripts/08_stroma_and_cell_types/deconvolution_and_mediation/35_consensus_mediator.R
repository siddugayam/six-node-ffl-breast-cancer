#!/usr/bin/env Rscript
## 35_consensus_mediator.R -- build two method-agnostic consensus fibroblast axes and
## append them to the mediator set, so the mediation can also be reported against an
## estimate that does not belong to any single method.
##   CONSENSUS_PC1      first principal component of the 23 z-scored estimates
##   CONSENSUS_rankmean mean within-cohort rank across the 23 estimates
## A collagen-safe consensus is also built from only those estimates whose signature
## contains neither COL1A1 nor COL3A1.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/35_consensus.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")
F0 <- readRDS(file.path(D, "fibroblast_estimates_all.rds"))
## drop the exact duplicate (my MCP-counter reproduction) so it is not double-weighted
F0 <- F0[, setdiff(colnames(F0), "MCPcounter_Fibroblasts_recomputed"), drop = FALSE]
logf("estimates entering the consensus: ", ncol(F0))
Z  <- scale(F0)
pc <- prcomp(Z, center = FALSE, scale. = FALSE)
v  <- pc$sdev^2 / sum(pc$sdev^2)
logf("variance explained: PC1 ", round(100 * v[1], 1), "%  PC2 ", round(100 * v[2], 1),
     "%  PC3 ", round(100 * v[3], 1), "%")
p1 <- pc$x[, 1]
if (cor(p1, F0[, "ESTIMATE_StromalScore"]) < 0) p1 <- -p1     # orient positively with stroma
logf("PC1 loadings (oriented so that higher = more stroma):")
ld <- pc$rotation[, 1] * ifelse(cor(pc$x[, 1], F0[, "ESTIMATE_StromalScore"]) < 0, -1, 1)
for (j in names(sort(ld, decreasing = TRUE))) logf(sprintf("  %-34s %+.3f", j, ld[j]))
rk <- apply(F0, 2, rank); rm1 <- rowMeans(rk)
CLEAN <- c("xCell_Fibroblasts","xCell_StromaScore","xCell_MicroenvironmentScore",
           "MCPcounter_Fibroblasts_nocol","EPIC_CAFs_nocol","ConsensusTME_Fibroblasts_nocol",
           "CBSX_Wu_nocol_CAFs","InstaPrism_Wu_nocollagen_CAFs","ssGSEA_ESTIMATE_stromal_clean",
           "meanZ_CAF_A_estimate","meanZ_CAF_scRNA","ssGSEA_CAF_scRNA","meanZ_FARMER_clean",
           "meanZ_CAF_B_markers","NNLS_Wu_CAFs")
CLEAN <- intersect(CLEAN, colnames(F0))
logf("\ncollagen-safe consensus built from ", length(CLEAN), " estimates: ", paste(CLEAN, collapse = ", "))
logf("  (a signature is 'collagen-safe' if it contains neither COL1A1 nor COL3A1; the ",
     "meanZ/ssGSEA sets already exclude COL* and network nodes)")
rm2 <- rowMeans(apply(F0[, CLEAN, drop = FALSE], 2, rank))
NEW <- cbind(CONSENSUS_PC1 = p1, CONSENSUS_rankmean = rm1, CONSENSUS_rankmean_collagensafe = rm2)
logf("\nSpearman between consensus axes: PC1 vs rankmean ", round(cor(p1, rm1, method = "spearman"), 4),
     " ; rankmean vs collagen-safe rankmean ", round(cor(rm1, rm2, method = "spearman"), 4))
gcol <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
for (j in colnames(NEW)) logf(sprintf("  %-34s rho COL1A1 %+.3f  rho COL3A1 %+.3f", j,
  cor(NEW[, j], gcol["COL1A1", rownames(F0)], method = "spearman"),
  cor(NEW[, j], gcol["COL3A1", rownames(F0)], method = "spearman")))
ALL <- cbind(readRDS(file.path(D, "fibroblast_estimates_all.rds")), NEW[rownames(F0), ])
saveRDS(ALL, file.path(D, "fibroblast_estimates_all.rds"))
fwrite(data.table(sample = rownames(ALL), as.data.table(ALL)),
       file.path(BASE, "results/v3/deconv_fibroblast_estimates.csv"))
logf("fibroblast estimate table now has ", ncol(ALL), " columns")
logf("DONE 35")
