## ===========================================================================
## PART A. PAM50 subtype stratification of the miR-29 / TF -> collagen axis
## Independent implementation, written fresh 2026-09-09.
## TCGA-BRCA (Xena) + METABRIC.
## ===========================================================================
suppressPackageStartupMessages({library(matrixStats)})
setwd("/path/to/revision")
OUT <- "results/v2"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
options(width=200)
set.seed(20260909)

## ------------------------------------------------------------------ utils --
sp <- function(x,y){                     # Spearman rho + p + n on complete pairs
  k <- is.finite(x) & is.finite(y); n <- sum(k)
  if (n < 8) return(c(rho=NA, p=NA, n=n))
  ct <- suppressWarnings(cor.test(x[k], y[k], method="spearman", exact=FALSE))
  c(rho=unname(ct$estimate), p=unname(ct$p.value), n=n)
}
pcor <- function(x,y,z){                 # partial Spearman (rank -> Pearson partial)
  k <- is.finite(x)&is.finite(y)&is.finite(z); n <- sum(k)
  if (n < 12) return(c(rho=NA,p=NA,n=n))
  rx <- rank(x[k]); ry <- rank(y[k]); rz <- rank(z[k])
  rxy <- cor(rx,ry); rxz <- cor(rx,rz); ryz <- cor(ry,rz)
  r <- (rxy - rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
  tt <- r*sqrt((n-3)/(1-r^2)); p <- 2*pt(-abs(tt), df=n-3)
  c(rho=r, p=p, n=n)
}
fz  <- function(r) 0.5*log((1+r)/(1-r))
## Cochran Q heterogeneity test on Fisher-z with variance 1/(n-3)
hetQ <- function(r, n){
  k <- is.finite(r) & is.finite(n) & n > 4
  r <- r[k]; n <- n[k]
  if (length(r) < 2) return(list(Q=NA,df=NA,p=NA,I2=NA,k=length(r),pooled=NA,pooled_se=NA))
  z <- fz(r); w <- n-3
  zbar <- sum(w*z)/sum(w)
  Q <- sum(w*(z-zbar)^2); df <- length(r)-1
  p <- pchisq(Q, df=df, lower.tail=FALSE)
  I2 <- max(0, (Q-df)/Q)*100
  list(Q=Q, df=df, p=p, I2=I2, k=length(r), pooled=tanh(zbar), pooled_se=1/sqrt(sum(w)))
}

## ------------------------------------------------------------------- data --
gexp  <- readRDS("data/brca_gene_expr.rds")
mexp  <- readRDS("data/brca_mirna_expr_canonical.rds")
pheno <- readRDS("data/brca_pheno.rds")
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors=FALSE)
cat("gene matrix:", dim(gexp), " dup rownames:", sum(duplicated(rownames(gexp))), "\n")
cat("miRNA matrix:", dim(mexp), " dup rownames:", sum(duplicated(rownames(mexp))), "\n")
stopifnot(sum(duplicated(rownames(gexp)))==0, sum(duplicated(rownames(mexp)))==0)

cl <- read.delim("data/brca_clinicalMatrix.tsv", stringsAsFactors=FALSE,
                 check.names=FALSE, quote="", na.strings=c("NA",""))
stopifnot("PAM50Call_RNAseq" %in% colnames(cl), "sampleID" %in% colnames(cl))
pam <- setNames(cl$PAM50Call_RNAseq, cl$sampleID)
cat("PAM50 calls in clinical matrix:\n"); print(table(pam, useNA="ifany"))

tum  <- pheno$sample[pheno$sample_type=="Primary Tumor"]
samp <- sort(intersect(intersect(colnames(gexp), colnames(mexp)), tum))
cat("\nprimary tumours with BOTH assays:", length(samp), "\n")
sub  <- pam[samp]
cat("of which with a PAM50 call:", sum(!is.na(sub)), "\n")
use  <- samp[!is.na(sub)]
sub  <- sub[!is.na(sub)]
sub[sub=="Normal"] <- "Normal-like"
cat("ANALYSIS SET n =", length(use), "\n"); print(table(sub))

G <- gexp[, use, drop=FALSE]; M <- mexp[, use, drop=FALSE]

## ------------------------------------------------------- CAF / stromal ------
gmt <- readLines(system.file("extdata","SI_geneset.gmt", package="estimate"))
strom <- strsplit(gmt[grep("^StromalSignature", gmt)], "\t")[[1]][-c(1,2)]
strom <- unique(strom[strom!=""])
cat("\nESTIMATE StromalSignature genes:", length(strom), "\n")
net <- unique(nodes$name)
drop_col  <- grep("^COL", strom, value=TRUE)
drop_net  <- intersect(strom, net)
cafA_genes <- setdiff(strom, union(drop_col, drop_net))
cafA_genes <- intersect(cafA_genes, rownames(G))
cat("dropped collagens:", length(drop_col), paste(drop_col, collapse=","), "\n")
cat("dropped network nodes:", length(drop_net), paste(drop_net, collapse=","), "\n")
cat("CAF-A genes used:", length(cafA_genes), "\n")
stopifnot(length(grep("^COL", cafA_genes))==0, length(intersect(cafA_genes,net))==0)

cafB_genes <- intersect(setdiff(c("DCN","LUM","FAP","PDGFRB","THY1","POSTN"), net), rownames(G))
cat("CAF-B genes used:", paste(cafB_genes, collapse=","), "\n")
stopifnot(length(intersect(cafB_genes,net))==0)

zscore <- function(X){ mu <- rowMeans(X); s <- rowSds(X); (X-mu)/s }
CAFA <- colMeans(zscore(G[cafA_genes,,drop=FALSE]))
CAFB <- colMeans(zscore(G[cafB_genes,,drop=FALSE]))
cat("cor(CAF-A,CAF-B) spearman =", round(cor(CAFA,CAFB,method="spearman"),4), "\n")

## ------------------------------------------------------------- axis list ----
mirs  <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-101","hsa-let-7b")
tfs   <- c("ETS1","NFKB1","SP1","RELA","MYC","TP53")
cols  <- c("COL1A1","COL3A1")
stopifnot(all(mirs %in% rownames(M)), all(tfs %in% rownames(G)), all(cols %in% rownames(G)))

axes <- rbind(
  expand.grid(src=mirs[1:3], tgt=cols, stringsAsFactors=FALSE),
  data.frame(src="hsa-miR-101", tgt="EZH2", stringsAsFactors=FALSE),
  data.frame(src="hsa-let-7b", tgt="COL3A1", stringsAsFactors=FALSE),
  expand.grid(src=tfs, tgt="COL1A1", stringsAsFactors=FALSE),
  expand.grid(src=tfs, tgt="COL3A1", stringsAsFactors=FALSE),
  data.frame(src="COL1A1", tgt="COL3A1", stringsAsFactors=FALSE),
  data.frame(src="CAF_A",  tgt=c("COL1A1","COL3A1"), stringsAsFactors=FALSE)
)
axes$src_type <- ifelse(grepl("^hsa-", axes$src), "miRNA",
                 ifelse(axes$src=="CAF_A","score","gene"))
vecOf <- function(nm, ty, idx){
  if (ty=="miRNA") M[nm, idx] else if (ty=="score") CAFA[idx] else G[nm, idx]
}

## ------------------------------------------------- per-subtype computation --
subtypes <- c("LumA","LumB","Her2","Basal","Normal-like")
res <- list()
for (r in seq_len(nrow(axes))){
  s <- axes$src[r]; t <- axes$tgt[r]; ty <- axes$src_type[r]
  if (!(t %in% rownames(G))) next
  ## ALL
  x <- vecOf(s,ty,seq_along(use)); y <- G[t, ]
  a <- sp(x,y); pA <- pcor(x,y,CAFA); pB <- pcor(x,y,CAFB)
  res[[length(res)+1]] <- data.frame(cohort="TCGA-BRCA", stratum="ALL", axis=paste(s,"->",t),
      source=s, target=t, source_type=ty, n=a["n"], rho=a["rho"], p=a["p"],
      rho_partial_CAFA=pA["rho"], p_partial_CAFA=pA["p"],
      rho_partial_CAFB=pB["rho"], p_partial_CAFB=pB["p"], stringsAsFactors=FALSE)
  for (st in subtypes){
    k <- which(sub==st)
    x <- vecOf(s,ty,k); y <- G[t, k]
    a <- sp(x,y); pA <- pcor(x,y,CAFA[k]); pB <- pcor(x,y,CAFB[k])
    res[[length(res)+1]] <- data.frame(cohort="TCGA-BRCA", stratum=st, axis=paste(s,"->",t),
        source=s, target=t, source_type=ty, n=a["n"], rho=a["rho"], p=a["p"],
        rho_partial_CAFA=pA["rho"], p_partial_CAFA=pA["p"],
        rho_partial_CAFB=pB["rho"], p_partial_CAFB=pB["p"], stringsAsFactors=FALSE)
  }
}
R <- do.call(rbind, res); rownames(R) <- NULL

## ------------------------------------------------------ heterogeneity test --
het <- list()
for (ax in unique(R$axis)){
  sr <- R[R$axis==ax & R$stratum %in% subtypes, ]
  h  <- hetQ(sr$rho, sr$n)
  h4 <- hetQ(sr$rho[sr$stratum!="Normal-like"], sr$n[sr$stratum!="Normal-like"])
  hp <- hetQ(sr$rho_partial_CAFA, sr$n)
  het[[length(het)+1]] <- data.frame(cohort="TCGA-BRCA", axis=ax,
     k_subtypes=h$k, Q=h$Q, df=h$df, p_heterogeneity=h$p, I2_pct=h$I2,
     pooled_rho=h$pooled,
     Q_4subtypes=h4$Q, p_het_4subtypes=h4$p,
     Q_partialCAFA=hp$Q, p_het_partialCAFA=hp$p,
     rho_min=min(sr$rho,na.rm=TRUE), rho_max=max(sr$rho,na.rm=TRUE),
     rho_Basal=sr$rho[sr$stratum=="Basal"], rho_LumA=sr$rho[sr$stratum=="LumA"],
     stringsAsFactors=FALSE)
}
H <- do.call(rbind, het); rownames(H) <- NULL
write.csv(H, file.path(OUT,"subtype_heterogeneity_v3.csv"), row.names=FALSE)
cat("\n=== HETEROGENEITY ACROSS PAM50 SUBTYPES (TCGA-BRCA) ===\n")
print(H[,c("axis","k_subtypes","Q","df","p_heterogeneity","I2_pct","pooled_rho","rho_Basal","rho_LumA","rho_min","rho_max")], digits=3)

## ------------------------------------------------- CAF score distribution ---
cafd <- do.call(rbind, lapply(c("ALL",subtypes), function(st){
  k <- if (st=="ALL") seq_along(use) else which(sub==st)
  data.frame(cohort="TCGA-BRCA", stratum=st, n=length(k),
    CAFA_mean=mean(CAFA[k]), CAFA_sd=sd(CAFA[k]), CAFA_median=median(CAFA[k]),
    CAFA_q25=quantile(CAFA[k],.25), CAFA_q75=quantile(CAFA[k],.75),
    CAFB_mean=mean(CAFB[k]), CAFB_sd=sd(CAFB[k]), CAFB_median=median(CAFB[k]),
    COL1A1_mean=mean(G["COL1A1",k]), COL3A1_mean=mean(G["COL3A1",k]),
    stringsAsFactors=FALSE)
}))
rownames(cafd) <- NULL
kwA <- kruskal.test(CAFA ~ factor(sub)); kwB <- kruskal.test(CAFB ~ factor(sub))
cat("\n=== CAF score by subtype ===\n"); print(cafd, digits=3)
cat("Kruskal-Wallis CAF-A across 5 subtypes: chi2 =", round(kwA$statistic,3),
    "df =", kwA$parameter, "p =", format.pval(kwA$p.value,digits=3), "\n")
cat("Kruskal-Wallis CAF-B across 5 subtypes: chi2 =", round(kwB$statistic,3),
    "df =", kwB$parameter, "p =", format.pval(kwB$p.value,digits=3), "\n")
cafd$kruskal_p_CAFA <- kwA$p.value; cafd$kruskal_p_CAFB <- kwB$p.value
write.csv(cafd, file.path(OUT,"subtype_caf_distribution_v3.csv"), row.names=FALSE)

## ============================== METABRIC ====================================
cat("\n\n=================== METABRIC ===================\n")
mb <- readRDS("data/metabric.rds"); MB <- mb$M; mcl <- as.data.frame(mb$cl)
cat("METABRIC expression:", dim(MB), "\n")
stopifnot(sum(duplicated(rownames(MB)))==0)
sid <- mcl$SAMPLE_ID; msub <- mcl$CLAUDIN_SUBTYPE
names(msub) <- sid
msub <- msub[colnames(MB)]
cat("CLAUDIN_SUBTYPE:\n"); print(table(msub, useNA="ifany"))
keepM <- which(msub %in% c("LumA","LumB","Her2","Basal","Normal","claudin-low"))
MBs <- MB[, keepM, drop=FALSE]; msub2 <- msub[keepM]
msub2[msub2=="Normal"] <- "Normal-like"
cat("METABRIC analysis set n =", ncol(MBs), "\n"); print(table(msub2))

zM <- function(X){ mu <- rowMeans(X, na.rm=TRUE); s <- rowSds(X, na.rm=TRUE); (X-mu)/s }
cafA_mb <- intersect(cafA_genes, rownames(MBs))
cafB_mb <- intersect(cafB_genes, rownames(MBs))
cat("METABRIC CAF-A genes:", length(cafA_mb), " CAF-B genes:", paste(cafB_mb,collapse=","), "\n")
MCAFA <- colMeans(zM(MBs[cafA_mb,,drop=FALSE]), na.rm=TRUE)
MCAFB <- colMeans(zM(MBs[cafB_mb,,drop=FALSE]), na.rm=TRUE)
cat("METABRIC cor(CAF-A,CAF-B) =", round(cor(MCAFA,MCAFB,method="spearman"),4), "\n")

mb_axes <- rbind(expand.grid(src=tfs, tgt=c("COL1A1","COL3A1"), stringsAsFactors=FALSE),
                 data.frame(src="COL1A1", tgt="COL3A1", stringsAsFactors=FALSE),
                 data.frame(src="CAF_A", tgt=c("COL1A1","COL3A1"), stringsAsFactors=FALSE))
mb_st <- c("LumA","LumB","Her2","Basal","Normal-like","claudin-low")
resM <- list()
for (r in seq_len(nrow(mb_axes))){
  s <- mb_axes$src[r]; t <- mb_axes$tgt[r]
  if (!(t %in% rownames(MBs))) next
  if (s!="CAF_A" && !(s %in% rownames(MBs))) next
  getv <- function(idx) if (s=="CAF_A") MCAFA[idx] else MBs[s, idx]
  for (st in c("ALL", mb_st)){
    k <- if (st=="ALL") seq_len(ncol(MBs)) else which(msub2==st)
    if (length(k) < 10) next
    x <- getv(k); y <- MBs[t,k]
    a <- sp(x,y); pA <- pcor(x,y,MCAFA[k]); pB <- pcor(x,y,MCAFB[k])
    resM[[length(resM)+1]] <- data.frame(cohort="METABRIC", stratum=st, axis=paste(s,"->",t),
      source=s, target=t, source_type=if(s=="CAF_A")"score" else "gene",
      n=a["n"], rho=a["rho"], p=a["p"],
      rho_partial_CAFA=pA["rho"], p_partial_CAFA=pA["p"],
      rho_partial_CAFB=pB["rho"], p_partial_CAFB=pB["p"], stringsAsFactors=FALSE)
  }
}
RM <- do.call(rbind, resM); rownames(RM) <- NULL

hetM <- list()
for (ax in unique(RM$axis)){
  sr <- RM[RM$axis==ax & RM$stratum!="ALL", ]
  h <- hetQ(sr$rho, sr$n); hp <- hetQ(sr$rho_partial_CAFA, sr$n)
  hetM[[length(hetM)+1]] <- data.frame(cohort="METABRIC", axis=ax, k_subtypes=h$k,
    Q=h$Q, df=h$df, p_heterogeneity=h$p, I2_pct=h$I2, pooled_rho=h$pooled,
    Q_4subtypes=NA, p_het_4subtypes=NA,
    Q_partialCAFA=hp$Q, p_het_partialCAFA=hp$p,
    rho_min=min(sr$rho,na.rm=TRUE), rho_max=max(sr$rho,na.rm=TRUE),
    rho_Basal=sr$rho[sr$stratum=="Basal"][1], rho_LumA=sr$rho[sr$stratum=="LumA"][1],
    stringsAsFactors=FALSE)
}
HM <- do.call(rbind, hetM); rownames(HM) <- NULL
write.csv(rbind(H,HM), file.path(OUT,"subtype_heterogeneity_v3.csv"), row.names=FALSE)
cat("\n=== METABRIC heterogeneity ===\n")
print(HM[,c("axis","k_subtypes","Q","p_heterogeneity","I2_pct","pooled_rho","rho_Basal","rho_LumA")], digits=3)

cafdM <- do.call(rbind, lapply(c("ALL", mb_st), function(st){
  k <- if (st=="ALL") seq_len(ncol(MBs)) else which(msub2==st)
  data.frame(cohort="METABRIC", stratum=st, n=length(k),
    CAFA_mean=mean(MCAFA[k]), CAFA_sd=sd(MCAFA[k]), CAFA_median=median(MCAFA[k]),
    CAFA_q25=quantile(MCAFA[k],.25), CAFA_q75=quantile(MCAFA[k],.75),
    CAFB_mean=mean(MCAFB[k]), CAFB_sd=sd(MCAFB[k]), CAFB_median=median(MCAFB[k]),
    COL1A1_mean=mean(MBs["COL1A1",k]), COL3A1_mean=mean(MBs["COL3A1",k]),
    stringsAsFactors=FALSE)
}))
rownames(cafdM) <- NULL
kwAM <- kruskal.test(MCAFA ~ factor(msub2))
cafdM$kruskal_p_CAFA <- kwAM$p.value; cafdM$kruskal_p_CAFB <- kruskal.test(MCAFB~factor(msub2))$p.value
cat("\n=== METABRIC CAF by subtype ===\n"); print(cafdM, digits=3)
cat("METABRIC Kruskal-Wallis CAF-A: chi2 =", round(kwAM$statistic,3), "p =",
    format.pval(kwAM$p.value,digits=3), "\n")
write.csv(rbind(cafd,cafdM), file.path(OUT,"subtype_caf_distribution_v3.csv"), row.names=FALSE)

## ---------------------------------------------------------------- output ---
ALL <- rbind(R, RM)
ALL$q_BH <- p.adjust(ALL$p, method="BH")
write.csv(ALL, file.path(OUT,"subtype_stratified.csv"), row.names=FALSE)
cat("\nwrote", file.path(OUT,"subtype_stratified.csv"), "rows =", nrow(ALL), "\n")

cat("\n=== KEY TABLE: TCGA-BRCA per subtype ===\n")
key <- ALL[ALL$cohort=="TCGA-BRCA" &
  ALL$axis %in% c("hsa-miR-29a -> COL1A1","hsa-miR-29b -> COL1A1","hsa-miR-29c -> COL1A1",
                  "hsa-miR-29a -> COL3A1","hsa-miR-29b -> COL3A1","hsa-miR-29c -> COL3A1",
                  "ETS1 -> COL1A1","NFKB1 -> COL1A1","SP1 -> COL1A1","RELA -> COL1A1",
                  "COL1A1 -> COL3A1"), ]
print(key[,c("stratum","axis","n","rho","p","rho_partial_CAFA")], digits=3)
cat("\n=== KEY TABLE: METABRIC per subtype ===\n")
keyM <- ALL[ALL$cohort=="METABRIC" & ALL$axis %in%
  c("ETS1 -> COL1A1","NFKB1 -> COL1A1","SP1 -> COL1A1","RELA -> COL1A1","COL1A1 -> COL3A1"),]
print(keyM[,c("stratum","axis","n","rho","p","rho_partial_CAFA")], digits=3)
cat("\nDONE\n")
