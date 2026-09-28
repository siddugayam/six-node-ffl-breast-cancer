## Independent re-derivation of the TCGA-BRCA differential expression tables (task A)
## Author: verification run.  Writes to results/v2/.
suppressPackageStartupMessages({library(limma)})
setwd("/path/to/revision")
OUT <- "results/v2"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

pheno <- readRDS("data/brca_pheno.rds")

run_de <- function(mat, label){
  st <- pheno$sample_type[match(colnames(mat), pheno$sample)]
  stopifnot(!any(is.na(st)))
  keep_s <- st %in% c("Primary Tumor","Solid Tissue Normal")
  mat <- mat[, keep_s, drop=FALSE]; st <- st[keep_s]
  grp <- factor(ifelse(st=="Primary Tumor","Tumor","Normal"), levels=c("Normal","Tumor"))
  cat(sprintf("[%s] samples: Tumor=%d Normal=%d\n", label, sum(grp=="Tumor"), sum(grp=="Normal")))
  ## drop constant / all-zero features only
  rv <- apply(mat, 1, function(x) var(x, na.rm=TRUE))
  keep_f <- is.finite(rv) & rv > 0
  cat(sprintf("[%s] features in: %d  tested: %d  dropped(zero variance): %d\n",
              label, nrow(mat), sum(keep_f), sum(!keep_f)))
  mat <- mat[keep_f, , drop=FALSE]
  design <- model.matrix(~grp)
  fit <- eBayes(lmFit(mat, design))
  tt <- topTable(fit, coef="grpTumor", number=Inf, sort.by="P")
  ## ---- HARD ASSERTION: rownames must be feature names, not integer row indices ----
  rn <- rownames(tt)
  if (all(grepl("^[0-9]+$", rn))) stop("topTable returned INTEGER ROW INDICES for ", label)
  ov <- mean(rn %in% rownames(mat))
  cat(sprintf("[%s] ASSERT: fraction of topTable rownames matching expression rownames = %.4f\n", label, ov))
  stopifnot(ov == 1)
  tt$feature <- rn
  tt <- tt[, c("feature","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
  tt
}

gexp <- readRDS("data/brca_gene_expr.rds")
cat("gene matrix:", dim(gexp), " duplicated rownames:", sum(duplicated(rownames(gexp))), "\n")
mexp <- readRDS("data/brca_mirna_expr_canonical.rds")
cat("miRNA matrix:", dim(mexp), " duplicated rownames:", sum(duplicated(rownames(mexp))), "\n")

de_g <- run_de(gexp, "gene")
de_m <- run_de(mexp, "miRNA")

thr <- function(tt, fdr=0.05, lfc=1) sum(tt$adj.P.Val < fdr & abs(tt$logFC) > lfc)
for (L in c(0,0.5,1,1.5,2)) cat(sprintf("gene  sig @FDR<0.05 |logFC|>%.1f : %d\n", L, thr(de_g,0.05,L)))
for (L in c(0,0.5,1,1.5,2)) cat(sprintf("miRNA sig @FDR<0.05 |logFC|>%.1f : %d\n", L, thr(de_m,0.05,L)))
cat("gene  sig @FDR<0.01 |logFC|>1 :", sum(de_g$adj.P.Val<0.01 & abs(de_g$logFC)>1), "\n")
cat("miRNA sig @FDR<0.01 |logFC|>1 :", sum(de_m$adj.P.Val<0.01 & abs(de_m$logFC)>1), "\n")
cat("gene  sig @FDR<0.05 |logFC|>=1:", sum(de_g$adj.P.Val<0.05 & abs(de_g$logFC)>=1), "\n")
cat("miRNA sig @FDR<0.05 |logFC|>=1:", sum(de_m$adj.P.Val<0.05 & abs(de_m$logFC)>=1), "\n")

write.csv(de_g, file.path(OUT,"v2_DE_genes.csv"), row.names=FALSE)
write.csv(de_m, file.path(OUT,"v2_DE_mirnas.csv"), row.names=FALSE)
saveRDS(list(gene=de_g, mirna=de_m), file.path(OUT,"v2_DE.rds"))
cat("DONE\n")
