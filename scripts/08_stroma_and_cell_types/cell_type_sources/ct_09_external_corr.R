#!/usr/bin/env Rscript
# Independent-cohort replication of (i) the collagen/TF correlations and
# (ii) the miR-29 -> collagen axis, each with and without adjustment for CAF content.
suppressMessages({library(ppcor)})
set.seed(42)
RES  <- "/path/to/revision/results/multiomics/"
DATA <- "/path/to/revision/data/"
CACHE<- "/path/to/revision/cache/celltype/"
EXT  <- "/path/to/revision/cache/external2/"
CAFSIG <- readLines(paste0(CACHE,"caf_signature_scrna.txt"))

gene_pairs <- data.frame(
  source=c("COL1A1","NFKB1","RELA","SP1","ETS1","NFKB1","RELA","SP1","ETS1","EZH2","MYC"),
  target=c("COL3A1","COL1A1","COL1A1","COL1A1","COL1A1","COL3A1","COL3A1","COL3A1","COL3A1","COL1A1","COL1A1"),
  stringsAsFactors=FALSE)
mir_pairs <- data.frame(
  source=c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c",
           "hsa-let-7b","hsa-let-7e","hsa-miR-101"),
  target=c("COL1A1","COL1A1","COL1A1","COL3A1","COL3A1","COL3A1","COL3A1","COL3A1","EZH2"),
  stringsAsFactors=FALSE)

cafscore <- function(X) {
  g <- intersect(CAFSIG, rownames(X))
  Z <- t(scale(t(X[g,,drop=FALSE])))
  list(score=colMeans(Z, na.rm=TRUE), n_used=length(g), used=g)
}
docor <- function(cohort, platform, X, Y=NULL, pairs, caf) {
  out <- list()
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$source[i]; b <- pairs$target[i]
    xm <- if (!is.null(Y) && a %in% rownames(Y)) Y else X
    if (!(a %in% rownames(xm)) || !(b %in% rownames(X))) {
      out[[length(out)+1]] <- data.frame(cohort=cohort, platform=platform, source=a, target=b,
        n=NA, rho=NA, p=NA, rho_CAFadj=NA, p_CAFadj=NA, note="feature absent", stringsAsFactors=FALSE); next }
    x <- as.numeric(xm[a,]); y <- as.numeric(X[b,]); z <- caf
    ok <- complete.cases(x,y,z)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    pc <- suppressWarnings(pcor.test(x[ok], y[ok], z[ok], method="spearman"))
    out[[length(out)+1]] <- data.frame(cohort=cohort, platform=platform, source=a, target=b,
      n=sum(ok), rho=unname(ct$estimate), p=ct$p.value,
      rho_CAFadj=pc$estimate, p_CAFadj=pc$p.value, note="", stringsAsFactors=FALSE)
  }
  do.call(rbind, out)
}
ALL <- list(); INFO <- list()

## ---------------- TCGA-BRCA reference ----------------
expr <- readRDS(paste0(DATA,"brca_gene_expr.rds")); mirT <- readRDS(paste0(DATA,"brca_mirna_expr_canonical.rds"))
ph <- readRDS(paste0(DATA,"brca_pheno.rds")); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
sg <- intersect(colnames(expr), tum); sb <- intersect(sg, colnames(mirT))
cf <- cafscore(expr[,sg]); cat("TCGA CAF signature genes used:", cf$n_used, "\n")
ALL[[1]] <- docor("TCGA-BRCA (reference)","Illumina HiSeq RNA-seq", expr[,sg], NULL, gene_pairs, cf$score)
cfb <- cafscore(expr[,sb])
ALL[[2]] <- docor("TCGA-BRCA (reference)","Illumina HiSeq RNA-seq", expr[,sb], mirT[,sb], mir_pairs, cfb$score)
INFO[[1]] <- data.frame(cohort="TCGA-BRCA (reference)", accession="TCGA-BRCA",
  platform="Illumina HiSeq RNA-seq + miRNA-seq", n_tumours=length(sg), n_with_miRNA=length(sb),
  caf_sig_genes_used=cf$n_used, stringsAsFactors=FALSE)

## ---------------- METABRIC ----------------
mb <- readRDS(paste0(DATA,"metabric.rds")); M <- mb$M
cat("METABRIC:", dim(M), "\n")
cfm <- cafscore(M); cat("METABRIC CAF genes used:", cfm$n_used, "\n")
ALL[[3]] <- docor("METABRIC","Illumina HT-12 v3 microarray", M, NULL, gene_pairs, cfm$score)
INFO[[2]] <- data.frame(cohort="METABRIC", accession="cBioPortal brca_metabric",
  platform="Illumina HT-12 v3 microarray", n_tumours=ncol(M), n_with_miRNA=0,
  caf_sig_genes_used=cfm$n_used, stringsAsFactors=FALSE)

## ---------------- MET500 breast metastases ----------------
me <- readRDS(paste0(DATA,"met500_breast.rds")); ME <- me$M
cfe <- cafscore(ME); cat("MET500 CAF genes used:", cfe$n_used, " n:", ncol(ME), "\n")
ALL[[4]] <- docor("MET500 (breast metastases)","RNA-seq log2", ME, NULL, gene_pairs, cfe$score)
INFO[[3]] <- data.frame(cohort="MET500 (breast metastases)", accession="MET500",
  platform="RNA-seq (poly-A / exome capture), log2", n_tumours=ncol(ME), n_with_miRNA=0,
  caf_sig_genes_used=cfe$n_used, stringsAsFactors=FALSE)

## ---------------- GSE19783 : matched mRNA + miRNA ----------------
readsm <- function(f) {
  L <- readLines(gzfile(f)); i <- grep('^"ID_REF"', L); j <- grep('^!series_matrix_table_end', L)
  tb <- read.delim(text=L[i:(j-1)], row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
  ti <- grep('^!Sample_title', L)[1]; ga <- grep('^!Sample_geo_accession', L)[1]
  titles <- gsub('"','', strsplit(L[ti],"\t")[[1]][-1]); accs <- gsub('"','', strsplit(L[ga],"\t")[[1]][-1])
  list(tab=as.matrix(tb), title=setNames(titles, accs))
}
mrna <- readsm(paste0(EXT,"GSE19783-GPL6480_series_matrix.txt.gz"))
mirn <- readsm(paste0(EXT,"GSE19783-GPL8227_series_matrix.txt.gz"))
cat("GSE19783 mRNA:", dim(mrna$tab), " miRNA:", dim(mirn$tab), "\n")
an <- read.delim(gzfile(paste0(EXT,"GPL6480.annot.gz")), skip=grep("^ID\t", readLines(gzfile(paste0(EXT,"GPL6480.annot.gz")), n=200))-1,
                 header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)
sym <- setNames(an$Gene.symbol, an$ID); sym <- sym[sym!="" & !is.na(sym)]
sym <- sapply(strsplit(sym, "///"), `[`, 1)
pr <- intersect(rownames(mrna$tab), names(sym))
Xg <- mrna$tab[pr,]; g <- sym[pr]
Xg <- rowsum(Xg, g, na.rm=TRUE) / as.numeric(table(g)[rownames(rowsum(Xg,g))])
cat("GSE19783 genes after probe collapse:", nrow(Xg), "\n")
key <- function(s) sub("^(BreastCancer-[0-9]+).*", "\\1", s)
km <- key(mrna$title); ki <- key(mirn$title)
common <- intersect(km, ki); cat("GSE19783 matched mRNA/miRNA samples:", length(common), "\n")
Xg2 <- Xg[, names(km)[match(common, km)]]; Xi <- mirn$tab[, names(ki)[match(common, ki)]]
colnames(Xg2) <- common; colnames(Xi) <- common
Xi <- Xi[rowSums(!is.na(Xi))>=0.8*ncol(Xi), , drop=FALSE]
cat("GSE19783 miRNA features with >=80% non-missing:", nrow(Xi), "\n")
cfg <- cafscore(Xg2); cat("GSE19783 CAF genes used:", cfg$n_used, "\n")
ALL[[5]] <- docor("GSE19783","Agilent 44k mRNA (GPL6480)", Xg2, NULL, gene_pairs, cfg$score)
ALL[[6]] <- docor("GSE19783","Agilent miRNA v2 (GPL8227) + GPL6480", Xg2, Xi, mir_pairs, cfg$score)
INFO[[4]] <- data.frame(cohort="GSE19783", accession="GSE19783",
  platform="Agilent GPL6480 mRNA + Agilent GPL8227 miRNA", n_tumours=ncol(Xg2), n_with_miRNA=ncol(Xi),
  caf_sig_genes_used=cfg$n_used, stringsAsFactors=FALSE)

res <- do.call(rbind, ALL)
res$rho <- signif(res$rho,4); res$rho_CAFadj <- signif(res$rho_CAFadj,4)
write.csv(res, paste0(RES,"external_cohort_correlations.csv"), row.names=FALSE)
info <- do.call(rbind, INFO); write.csv(info, paste0(RES,"external_cohort_info.csv"), row.names=FALSE)
cat("\nexternal_cohort_correlations.csv rows:", nrow(res), "\n")
print(res[,c("cohort","source","target","n","rho","p","rho_CAFadj","p_CAFadj")], row.names=FALSE, digits=3)
saveRDS(list(GSE19783_mRNA=Xg2, GSE19783_miRNA=Xi), paste0(EXT,"gse19783.rds"))
