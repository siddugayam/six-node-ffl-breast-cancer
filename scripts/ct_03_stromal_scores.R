#!/usr/bin/env Rscript
# Stromal / fibroblast / epithelial content scores for TCGA-BRCA primary tumours.
# Deliberately provides BOTH the canonical ESTIMATE stromal score (which contains COL3A1)
# and collagen-free variants, so the stromal adjustment of the collagen axis is not circular.
suppressMessages({library(estimate)})
set.seed(42)
RES <- "/path/to/revision/results/multiomics/"
CACHE <- "/path/to/revision/cache/celltype/"
DATA <- "/path/to/revision/data/"

expr  <- readRDS(paste0(DATA,"brca_gene_expr.rds"))
pheno <- readRDS(paste0(DATA,"brca_pheno.rds"))
tum <- pheno$sample[pheno$sample_type == "Primary Tumor"]
expr <- expr[, colnames(expr) %in% tum, drop=FALSE]
cat("tumour samples in expression matrix:", ncol(expr), " genes:", nrow(expr), "\n")

## ---------- 1. canonical ESTIMATE (package) ----------
gct_in  <- paste0(CACHE,"tcga_brca_expr_for_estimate.txt")
write.table(data.frame(GeneSymbol=rownames(expr), expr, check.names=FALSE),
            gct_in, sep="\t", quote=FALSE, row.names=FALSE)
filterCommonGenes(input.f=gct_in, output.f=paste0(CACHE,"est_filtered.gct"), id="GeneSymbol")
estimateScore(input.ds=paste0(CACHE,"est_filtered.gct"),
              output.ds=paste0(CACHE,"est_score.gct"), platform="illumina")
es <- read.delim(paste0(CACHE,"est_score.gct"), skip=2, row.names=1, check.names=FALSE)
es <- es[,-1, drop=FALSE]
est <- as.data.frame(t(es)); colnames(est) <- rownames(es)
rownames(est) <- gsub("\\.", "-", rownames(est))
cat("ESTIMATE output rows:", nrow(est), " cols:", paste(colnames(est), collapse=","), "\n")

## ---------- 2. own ssGSEA implementation (so gene sets can be edited) ----------
gmt <- readLines(system.file("extdata","SI_geneset.gmt", package="estimate"))
sets <- lapply(strsplit(gmt,"\t"), function(x) x[-c(1,2)])
names(sets) <- sapply(strsplit(gmt,"\t"), `[`, 1)
cat("gene sets:", names(sets), lengths(sets), "\n")

ssgsea <- function(X, gs, alpha=0.25) {
  # Barbie et al. ssGSEA, as used inside estimate::estimateScore
  gs <- intersect(gs, rownames(X)); if (length(gs) < 5) return(rep(NA_real_, ncol(X)))
  N <- nrow(X)
  apply(X, 2, function(v) {
    ord <- order(v, decreasing=TRUE)
    rnk <- rep(0, N); rnk[ord] <- N:1                 # rank, largest expression = N
    inset <- rownames(X) %in% gs
    ord2 <- order(rnk, decreasing=TRUE)
    hit <- inset[ord2]
    w <- (rnk[ord2] * hit)^alpha
    Pw  <- cumsum(w)/sum(w)
    Pnw <- cumsum(!hit)/sum(!hit)
    sum(Pw - Pnw)
  })
}
X <- expr[rowSums(is.na(expr))==0 & apply(expr,1,sd)>0, , drop=FALSE]
cat("genes used for ssGSEA:", nrow(X), "\n")

col_in_stromal <- grep("^COL", sets$StromalSignature, value=TRUE)
cat("collagen genes inside ESTIMATE StromalSignature:", paste(col_in_stromal, collapse=","), "\n")
sets$StromalSignature_noCOL <- setdiff(sets$StromalSignature, col_in_stromal)

sc <- list()
for (nm in c("StromalSignature","ImmuneSignature","StromalSignature_noCOL")) {
  sc[[nm]] <- ssgsea(X, sets[[nm]]); cat("ssGSEA done:", nm, "n genes matched:",
      length(intersect(sets[[nm]], rownames(X))), "\n")
}
scores <- data.frame(sample=colnames(X), ssgsea_stromal=sc$StromalSignature,
                     ssgsea_immune=sc$ImmuneSignature,
                     ssgsea_stromal_noCOL=sc$StromalSignature_noCOL,
                     stringsAsFactors=FALSE)
m <- match(scores$sample, rownames(est))
scores$ESTIMATE_StromalScore <- est$StromalScore[m]
scores$ESTIMATE_ImmuneScore  <- est$ImmuneScore[m]
scores$ESTIMATE_Score        <- est$ESTIMATEScore[m]
cat("VALIDATION own ssGSEA vs package StromalScore: spearman rho =",
    cor(scores$ssgsea_stromal, scores$ESTIMATE_StromalScore, method="spearman", use="complete.obs"), "\n")
cat("stromal vs stromal_noCOL: spearman rho =",
    cor(scores$ssgsea_stromal, scores$ssgsea_stromal_noCOL, method="spearman"), "\n")

## ---------- 3. marker-based scores (all collagen-free) ----------
zmean <- function(genes, X) {
  g <- intersect(genes, rownames(X))
  if (!length(g)) return(list(score=rep(NA_real_, ncol(X)), used=character(0)))
  Z <- t(scale(t(X[g,,drop=FALSE])))
  list(score=colMeans(Z, na.rm=TRUE), used=g)
}
panels <- list(
  CAF_marker = c("FAP","PDGFRA","PDGFRB","THY1","DCN","LUM","FBLN1","VCAN","SULF1","ISLR","CDH11","MMP2"),
  epithelial = c("EPCAM","KRT8","KRT18","KRT19","CDH1","KRT7"),
  immune     = c("PTPRC","CD2","CD3E","CD53","LCK","CD48"),
  endothelial= c("PECAM1","VWF","CDH5","CLDN5","ENG")
)
for (nm in names(panels)) {
  r <- zmean(panels[[nm]], X); scores[[paste0("score_",nm)]] <- r$score
  cat("panel", nm, "genes used:", paste(r$used, collapse=","), "\n")
}

## ---------- 4. scRNA-derived CAF signature (collagen-free), if available ----------
sig_f <- paste0(CACHE,"caf_signature_scrna.txt")
if (file.exists(sig_f)) {
  sig <- readLines(sig_f)
  r <- zmean(sig, X); scores$score_CAF_scRNA <- r$score
  cat("scRNA CAF signature genes supplied:", length(sig), " used:", length(r$used), "\n")
  writeLines(r$used, paste0(CACHE,"caf_signature_scrna_used.txt"))
} else cat("NOTE: scRNA CAF signature not yet available; column omitted\n")

write.csv(scores, paste0(RES,"tcga_stromal_scores.csv"), row.names=FALSE)
saveRDS(scores, paste0(CACHE,"tcga_stromal_scores.rds"))
cat("rows written:", nrow(scores), "cols:", ncol(scores), "\n")
print(round(cor(scores[,-1], method="spearman", use="pairwise.complete.obs"), 3))
