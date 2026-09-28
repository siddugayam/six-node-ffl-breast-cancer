## 28_metabolic_ccle.R -----------------------------------------------------
## Stroma-free replication of part C: does miR-130a track the same metabolic
## programme across CCLE breast cancer cell lines (no stroma, no immune cells)?
## This is the version that actually predicts what an in vitro knockdown does.
suppressPackageStartupMessages({library(data.table); library(GSVA); library(BiocParallel)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
mir <- fread(file.path(ROOT,"results/v3/mir130a_ccle_breast_lines.csv"))
cat("CCLE breast lines with miR-130a:", nrow(mir), " range:", range(mir$val), "\n")

hdr <- fread(cmd=paste0("zcat ", shQuote(file.path(ROOT,"data/CCLE_rnaseq_tpm.txt.gz")), " | head -1"),
             header=FALSE, sep="\t")
cn <- as.character(hdr[1,])
br <- grep("_BREAST$", cn, value=TRUE)
cat("breast columns in CCLE RNA-seq:", length(br), "; overlap with miRNA file:",
    length(intersect(br, mir$line)), "\n")
sel <- c("gene_id", intersect(br, mir$line))
X <- fread(file.path(ROOT,"data/CCLE_rnaseq_tpm.txt.gz"), select=sel, showProgress=FALSE)
cat("expr loaded:", dim(X), "\n")
suppressPackageStartupMessages(library(org.Hs.eg.db))
ens <- sub("\\..*$","", X$gene_id)
m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=unique(ens), keytype="ENSEMBL", columns="SYMBOL"))
m <- as.data.table(m)[!is.na(SYMBOL)][, .SD[1], by=ENSEMBL]
sym <- setNames(m$SYMBOL, m$ENSEMBL)[ens]
M <- as.matrix(X[, -1]); rownames(M) <- sym
M <- M[!is.na(rownames(M)), , drop=FALSE]
M <- log2(M + 1)
## collapse duplicate symbols by highest mean expression
o <- order(rowMeans(M), decreasing=TRUE)
M <- M[o, ]; M <- M[!duplicated(rownames(M)), ]
stopifnot(!any(duplicated(rownames(M))))
sdv <- apply(M,1,sd); M <- M[sdv>0, ]
cat("symbol-keyed CCLE breast matrix:", dim(M), "\n")

sets <- readRDS(file.path(RES,"metabolic_genesets.rds"))
sets <- lapply(sets, intersect, rownames(M)); sets <- sets[lengths(sets)>=5]
ss <- gsva(ssgseaParam(M, sets, alpha=0.25, normalize=TRUE, minSize=5),
           BPPARAM=MulticoreParam(workers=8, RNGseed=1), verbose=FALSE)
cat("CCLE metabolic score matrix:", dim(ss), "\n")
saveRDS(ss, file.path(RES,"metabolic_ccle_ssgsea_scores.rds"))
fwrite(data.table(set_id=rownames(ss), as.data.table(ss)), file.path(RES,"metabolic_ccle_ssgsea_scores.csv"))

v <- setNames(log2(mir$val + 1), mir$line)[colnames(ss)]
stopifnot(!any(is.na(v)))
res <- rbindlist(lapply(rownames(ss), function(p) {
  ct <- suppressWarnings(cor.test(ss[p,], v, method="spearman", exact=FALSE))
  data.table(set_id=p, n_lines=length(v), rho=unname(ct$estimate), p=ct$p.value)
}))
res[, FDR := p.adjust(p,"BH")]
## TCGA comparison
tc <- fread(file.path(RES,"metabolic_feature_correlations.csv"))[
       score=="ssGSEA" & feature=="hsa-miR-130a-3p", .(set_id, rho_TCGA=rho, FDR_TCGA=FDR)]
res <- merge(res, tc, by="set_id", all.x=TRUE)
setnames(res, c("rho","p","FDR"), c("rho_CCLE","p_CCLE","FDR_CCLE"))
## per-line score percentiles for the two lines the group will use
for (ln in c("MDAMB231_BREAST","MCF7_BREAST")) {
  if (ln %in% colnames(ss))
    res[[paste0("pctile_", sub("_BREAST","",ln))]] <-
      sapply(res$set_id, function(p) round(100*mean(ss[p,] <= ss[p,ln]),1))
}
res <- res[order(-abs(rho_CCLE))]
fwrite(res, file.path(RES,"metabolic_ccle_mir130a_correlations.csv"))
cat("\nmiR-130a in the two planned lines: MDAMB231 =", mir[line=="MDAMB231_BREAST", val],
    " MCF7 =", mir[line=="MCF7_BREAST", val], "\n")
cat("\n-- CCLE: metabolic programmes tracking miR-130a across", length(v), "breast lines --\n")
print(res[order(-rho_CCLE)][1:15, .(set_id=substr(set_id,1,58), rho_CCLE=round(rho_CCLE,3),
      FDR_CCLE=signif(FDR_CCLE,3), rho_TCGA=round(rho_TCGA,3), pMB231=pctile_MDAMB231, pMCF7=pctile_MCF7)])
print(res[order(rho_CCLE)][1:15, .(set_id=substr(set_id,1,58), rho_CCLE=round(rho_CCLE,3),
      FDR_CCLE=signif(FDR_CCLE,3), rho_TCGA=round(rho_TCGA,3), pMB231=pctile_MDAMB231, pMCF7=pctile_MCF7)])
cat("\nsignificant at FDR<0.05 in CCLE:", res[FDR_CCLE<0.05, .N], "of", nrow(res), "\n")
ok <- res[!is.na(rho_TCGA)]
cat("TCGA vs CCLE rho agreement (all sets): Pearson", round(cor(ok$rho_CCLE, ok$rho_TCGA),3),
    "; Spearman", round(cor(ok$rho_CCLE, ok$rho_TCGA, method="spearman"),3),
    "; sign agreement", round(mean(sign(ok$rho_CCLE)==sign(ok$rho_TCGA)),3), "\n")
sig <- ok[FDR_TCGA<0.05]
cat("restricted to sets significant in TCGA (n=", nrow(sig), "): sign agreement",
    round(mean(sign(sig$rho_CCLE)==sign(sig$rho_TCGA)),3),
    "; binomial p =", signif(binom.test(sum(sign(sig$rho_CCLE)==sign(sig$rho_TCGA)), nrow(sig), 0.5)$p.value,3), "\n")
cat("DONE 28\n")
