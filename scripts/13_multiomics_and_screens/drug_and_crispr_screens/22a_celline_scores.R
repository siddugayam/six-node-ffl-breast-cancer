#!/usr/bin/env Rscript
## 22a) Cell-line module scores used for the drug-sensitivity analysis (Part 1C).
## Scores are built from DepMap 24Q4 expression (log2 TPM+1), z-scored across breast
## lines, then averaged over the gene set. Three scores:
##   COLLAGEN        : the fibrillar-collagen programme
##   MIR29_TARGET    : miR-29 targets anti-correlated with miR-29 in TCGA
##   MIR29_STRONG    : miR-29 targets in the STRONG_lowthroughput evidence tier
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

E <- fread(file.path(REV,"data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1_24Q4.csv"))
setnames(E, 1, "ModelID")
gs <- sub(" \\(\\d+\\)$","", setdiff(names(E),"ModelID"))
stopifnot(!any(duplicated(gs)))
msg("DepMap 24Q4 expression:", nrow(E), "models x", length(gs), "genes")
mod <- fread(file.path(REV,"data/depmap/Model_24Q4.csv"))
bre <- mod[OncotreeLineage=="Breast", .(ModelID, CellLineName, StrippedCellLineName, OncotreeSubtype, CCLEName, SangerModelID)]
Eb <- E[ModelID %in% bre$ModelID]
msg("breast models with expression:", nrow(Eb))
X <- as.matrix(Eb[, -1]); colnames(X) <- gs; rownames(X) <- Eb$ModelID
## keep genes expressed in breast lines
keep <- colMeans(X > 0.5) >= 0.25
X <- X[, keep, drop=FALSE]
msg("genes retained (log2TPM+1 > 0.5 in >=25% of breast lines):", ncol(X))
Z <- scale(X)                                            # z across breast lines, per gene
Z[!is.finite(Z)] <- NA

sets <- list(
  COLLAGEN       = c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL11A1"),
  MIR29_TARGET   = fread(file.path(OUT,"screens_mir29_anticorrelated_genes.csv"))$gene,
  MIR29_STRONG   = fread(file.path(OUT,"screens_mir29_target_set.csv"))[tier=="STRONG_lowthroughput", gene]
)
S <- data.table(ModelID=rownames(X))
for (nm in names(sets)) {
  g <- intersect(sets[[nm]], colnames(Z))
  msg("set", nm, ":", length(sets[[nm]]), "genes ->", length(g), "measured in breast lines")
  S[[nm]] <- rowMeans(Z[, g, drop=FALSE], na.rm=TRUE)
}
S <- merge(S, bre, by="ModelID")
setorder(S, -COLLAGEN)
fwrite(S, file.path(OUT,"screens_celline_module_scores.csv"))
msg("wrote screens_celline_module_scores.csv:", nrow(S), "breast lines")
msg("score correlations across breast lines (Spearman):")
print(round(cor(S[, .(COLLAGEN, MIR29_TARGET, MIR29_STRONG)], method="spearman"),3))
msg("top 5 COLLAGEN-high lines:"); print(head(S[, .(StrippedCellLineName, OncotreeSubtype, COLLAGEN, MIR29_TARGET)],5))
msg("bottom 5:"); print(tail(S[, .(StrippedCellLineName, OncotreeSubtype, COLLAGEN, MIR29_TARGET)],5))
msg("DONE 22a")
