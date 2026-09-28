## ============================================================================
## TASK A: PAM50 subtype stratification of the miR-29 -> collagen axis
##         and the TF -> collagen axis, in TCGA-BRCA; replication in METABRIC.
## Everything computed fresh in this run. No value carried over.
## ============================================================================
suppressPackageStartupMessages({library(matrixStats)})
setwd("/path/to/revision")
OUT <- "results/v2"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
set.seed(20260908)

cat("=========== PART A : PAM50 SUBTYPE STRATIFICATION ===========\n")

## ------------------------------------------------------------------ data --
gexp  <- readRDS("data/brca_gene_expr.rds")
mexp  <- readRDS("data/brca_mirna_expr_canonical.rds")
pheno <- readRDS("data/brca_pheno.rds")

stopifnot(sum(duplicated(rownames(gexp)))==0, sum(duplicated(rownames(mexp)))==0)
cat("gexp:", dim(gexp), " mexp:", dim(mexp), "\n")

tum  <- pheno$sample[pheno$sample_type=="Primary Tumor"]
samp <- sort(intersect(intersect(colnames(gexp), colnames(mexp)), tum))
cat("paired primary tumours (both assays): n =", length(samp), "\n")

## ------------------------------------------------------------ PAM50 calls -
CLIN <- "/path/to/home/Desktop/DD/R_GPR/ESIA_BRCA/brca_eisa_pilot/clinical/BRCA_clinicalMatrix"
cl <- read.delim(CLIN, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot("PAM50Call_RNAseq" %in% colnames(cl))
cl$PAM50Call_RNAseq[cl$PAM50Call_RNAseq==""] <- NA
cat("clinicalMatrix rows:", nrow(cl), " with a PAM50 call:",
    sum(!is.na(cl$PAM50Call_RNAseq)), "\n")
pam <- setNames(cl$PAM50Call_RNAseq, cl$sampleID)

## harmonise label: Xena writes 'Normal' for the Normal-like PAM50 centroid
pam50 <- pam[samp]
cat("\nPAM50 calls among the", length(samp), "paired primary tumours:\n")
print(table(pam50, useNA="ifany"))

keep <- !is.na(pam50)
sampP <- samp[keep]; pam50 <- pam50[keep]
pam50 <- factor(pam50, levels=c("LumA","LumB","Her2","Basal","Normal"))
levels(pam50)[levels(pam50)=="Normal"] <- "Normal-like"
cat("\nanalysed n =", length(sampP), "\n"); print(table(pam50))

G <- gexp[, sampP, drop=FALSE]; M <- mexp[, sampP, drop=FALSE]

## ------------------------------------------------------------- CAF score --
## reuse ONLY the gene-list definitions from the earlier verification run
## (they are definitions, not results); every score below is recomputed here.
cafdef <- readRDS(file.path(OUT,"v2_caf_scores.rds"))
sigA <- cafdef$sigA; sigB <- cafdef$sigB
stopifnot(length(sigA)==123, length(sigB)==4)
stopifnot(!any(grepl("^COL", sigA)), !any(grepl("^COL", sigB)))
cat("\nCAF-A signature genes:", length(sigA), " CAF-B:", paste(sigB, collapse=","), "\n")

zsc <- function(v) (v - mean(v))/sd(v)
mkscore <- function(sig, X){
  g <- intersect(sig, rownames(X))
  Z <- t(apply(X[g,,drop=FALSE], 1, zsc))
  colMeans(Z)
}
cafA <- mkscore(sigA, G); cafB <- mkscore(sigB, G)
cat("CAF-A genes found:", length(intersect(sigA, rownames(G))),
    " CAF-B genes found:", length(intersect(sigB, rownames(G))), "\n")

## ---------------------------------------------------------- helper stats --
sp  <- function(x,y) suppressWarnings(cor(x,y,method="spearman"))
spp <- function(x,y){ tt <- suppressWarnings(cor.test(x,y,method="spearman", exact=FALSE)); c(tt$estimate, tt$p.value) }
## partial Spearman: rank -> partial Pearson
pcor <- function(x,y,z){
  rx<-rank(x); ry<-rank(y); rz<-rank(z)
  rxy<-cor(rx,ry); rxz<-cor(rx,rz); ryz<-cor(ry,rz)
  (rxy - rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
}
getv <- function(nm){
  if (nm %in% rownames(M)) return(M[nm,])
  if (nm %in% rownames(G)) return(G[nm,])
  return(NULL)
}

## Cochran's Q heterogeneity test on Fisher-z transformed correlations
hetQ <- function(r, n){
  ok <- is.finite(r) & n>3
  r<-r[ok]; n<-n[ok]; k<-length(r)
  if(k<2) return(c(Q=NA,df=NA,p=NA,I2=NA,k=k))
  z <- atanh(r); w <- n-3
  zb <- sum(w*z)/sum(w)
  Q  <- sum(w*(z-zb)^2); df <- k-1
  p  <- pchisq(Q, df, lower.tail=FALSE)
  I2 <- max(0,(Q-df)/Q)*100
  c(Q=Q, df=df, p=p, I2=I2, k=k)
}
## two independent correlations (Fisher z)
zdiff <- function(r1,n1,r2,n2){
  z <- (atanh(r1)-atanh(r2))/sqrt(1/(n1-3)+1/(n2-3))
  c(z=z, p=2*pnorm(-abs(z)))
}

MIRS <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
TFS  <- c("ETS1","NFKB1","SP1","RELA")
COLS <- c("COL1A1","COL3A1")

cat("\nfeature availability: ",
    paste(sapply(c(MIRS,TFS,COLS), function(x) paste0(x,"=",!is.null(getv(x)))), collapse=" "), "\n")

## ------------------------------------------------------- per-subtype table -
rows <- list()
subs <- levels(pam50)
for (ax_s in c(MIRS,TFS)) for (ax_t in COLS){
  xs <- getv(ax_s); ys <- getv(ax_t)
  if (is.null(xs) || is.null(ys)) next
  if (ax_s %in% TFS && ax_t!="COL1A1") next          # TF axes: COL1A1 (per brief) + COL3A1 for completeness
  for (s in subs){
    i <- which(pam50==s); n <- length(i)
    if (n < 10) next
    e <- spp(xs[i], ys[i])
    pa <- pcor(xs[i], ys[i], cafA[i]); pb <- pcor(xs[i], ys[i], cafB[i])
    rows[[length(rows)+1]] <- data.frame(
      cohort="TCGA-BRCA", stratum=s, n=n, axis=paste0(ax_s," -> ",ax_t),
      source=ax_s, target=ax_t,
      rho=unname(e[1]), p=unname(e[2]),
      rho_partial_cafA=pa, rho_partial_cafB=pb,
      stringsAsFactors=FALSE)
  }
  ## ALL
  e <- spp(xs, ys)
  rows[[length(rows)+1]] <- data.frame(cohort="TCGA-BRCA", stratum="ALL", n=length(sampP),
    axis=paste0(ax_s," -> ",ax_t), source=ax_s, target=ax_t,
    rho=unname(e[1]), p=unname(e[2]),
    rho_partial_cafA=pcor(xs,ys,cafA), rho_partial_cafB=pcor(xs,ys,cafB),
    stringsAsFactors=FALSE)
}
## TF -> COL3A1 too
for (ax_s in TFS){
  xs <- getv(ax_s); ys <- getv("COL3A1")
  for (s in subs){
    i <- which(pam50==s); n<-length(i); if(n<10) next
    e <- spp(xs[i],ys[i])
    rows[[length(rows)+1]] <- data.frame(cohort="TCGA-BRCA", stratum=s, n=n,
      axis=paste0(ax_s," -> COL3A1"), source=ax_s, target="COL3A1",
      rho=unname(e[1]), p=unname(e[2]),
      rho_partial_cafA=pcor(xs[i],ys[i],cafA[i]), rho_partial_cafB=pcor(xs[i],ys[i],cafB[i]),
      stringsAsFactors=FALSE)
  }
  e <- spp(xs,ys)
  rows[[length(rows)+1]] <- data.frame(cohort="TCGA-BRCA", stratum="ALL", n=length(sampP),
    axis=paste0(ax_s," -> COL3A1"), source=ax_s, target="COL3A1",
    rho=unname(e[1]), p=unname(e[2]),
    rho_partial_cafA=pcor(xs,ys,cafA), rho_partial_cafB=pcor(xs,ys,cafB), stringsAsFactors=FALSE)
}
## COL1A1 ~ COL3A1
for (s in c(subs,"ALL")){
  i <- if (s=="ALL") seq_along(sampP) else which(pam50==s); n<-length(i); if(n<10) next
  e <- spp(G["COL1A1",i], G["COL3A1",i])
  rows[[length(rows)+1]] <- data.frame(cohort="TCGA-BRCA", stratum=s, n=n,
    axis="COL1A1 ~ COL3A1", source="COL1A1", target="COL3A1",
    rho=unname(e[1]), p=unname(e[2]),
    rho_partial_cafA=pcor(G["COL1A1",i],G["COL3A1",i],cafA[i]),
    rho_partial_cafB=pcor(G["COL1A1",i],G["COL3A1",i],cafB[i]), stringsAsFactors=FALSE)
}
ST <- do.call(rbind, rows)

cat("\n---------------- per-subtype correlations (TCGA-BRCA) ----------------\n")
for (a in unique(ST$axis)){
  d <- ST[ST$axis==a,]
  cat(sprintf("%-22s ", a))
  for (s in c(subs,"ALL")){
    r <- d[d$stratum==s,]
    if(nrow(r)) cat(sprintf("%s n=%d rho=%+.3f(p=%.1e)  ", s, r$n, r$rho, r$p))
  }
  cat("\n")
}

## --------------------------------------------------- heterogeneity tests --
het <- list()
for (a in unique(ST$axis)){
  d <- ST[ST$axis==a & ST$stratum %in% subs,]
  q <- hetQ(d$rho, d$n)
  ## Basal vs non-Basal (pooled complement, recomputed directly)
  parts <- strsplit(a," ~ | -> ")[[1]]
  xs <- getv(trimws(parts[1])); ys <- getv(trimws(parts[2]))
  ib <- which(pam50=="Basal"); io <- which(pam50!="Basal")
  rb <- sp(xs[ib],ys[ib]); ro <- sp(xs[io],ys[io])
  zd <- zdiff(rb,length(ib),ro,length(io))
  ## luminal (LumA+LumB) vs basal
  il <- which(pam50 %in% c("LumA","LumB"))
  rl <- sp(xs[il],ys[il]); zl <- zdiff(rb,length(ib),rl,length(il))
  het[[length(het)+1]] <- data.frame(cohort="TCGA-BRCA", axis=a,
    k_subtypes=unname(q["k"]), Q=unname(q["Q"]), df=unname(q["df"]),
    p_heterogeneity=unname(q["p"]), I2_pct=unname(q["I2"]),
    rho_Basal=rb, n_Basal=length(ib), rho_nonBasal=ro, n_nonBasal=length(io),
    z_Basal_vs_nonBasal=unname(zd["z"]), p_Basal_vs_nonBasal=unname(zd["p"]),
    rho_Luminal=rl, n_Luminal=length(il),
    z_Basal_vs_Luminal=unname(zl["z"]), p_Basal_vs_Luminal=unname(zl["p"]),
    stringsAsFactors=FALSE)
}
HET <- do.call(rbind, het)
HET$p_het_BH <- p.adjust(HET$p_heterogeneity, "BH")
HET$p_BvsNB_BH <- p.adjust(HET$p_Basal_vs_nonBasal, "BH")
cat("\n---------------- heterogeneity across the 5 PAM50 subtypes ----------------\n")
print(format(HET[,c("axis","k_subtypes","Q","df","p_heterogeneity","I2_pct",
                    "rho_Basal","rho_Luminal","p_Basal_vs_Luminal","p_het_BH")], digits=3))

## ------------------------------------------------------ CAF distribution --
cafrows <- list()
for (s in c(subs,"ALL")){
  i <- if(s=="ALL") seq_along(sampP) else which(pam50==s)
  cafrows[[length(cafrows)+1]] <- data.frame(cohort="TCGA-BRCA", stratum=s, n=length(i),
    cafA_mean=mean(cafA[i]), cafA_sd=sd(cafA[i]), cafA_median=median(cafA[i]),
    cafA_q25=quantile(cafA[i],.25), cafA_q75=quantile(cafA[i],.75),
    cafB_mean=mean(cafB[i]), cafB_sd=sd(cafB[i]), cafB_median=median(cafB[i]),
    COL1A1_mean=mean(G["COL1A1",i]), COL3A1_mean=mean(G["COL3A1",i]),
    stringsAsFactors=FALSE)
}
CAFD <- do.call(rbind, cafrows); rownames(CAFD) <- NULL
kwA <- kruskal.test(cafA ~ pam50); kwB <- kruskal.test(cafB ~ pam50)
cat("\n---------------- CAF score by subtype ----------------\n"); print(format(CAFD, digits=3))
cat(sprintf("Kruskal-Wallis CAF-A across subtypes: chi2=%.3f df=%d p=%.3g\n",
            kwA$statistic, kwA$parameter, kwA$p.value))
cat(sprintf("Kruskal-Wallis CAF-B across subtypes: chi2=%.3f df=%d p=%.3g\n",
            kwB$statistic, kwB$parameter, kwB$p.value))
cat("pairwise Wilcoxon CAF-A (BH):\n")
print(pairwise.wilcox.test(cafA, pam50, p.adjust.method="BH")$p.value)

## ================================================================ METABRIC =
cat("\n\n=========== METABRIC REPLICATION ===========\n")
mb <- readRDS("data/metabric.rds")
MB <- mb$M; mcl <- as.data.frame(mb$cl)
cat("METABRIC expr:", dim(MB), " clinical:", dim(mcl), "\n")
idcol <- if("SAMPLE_ID" %in% colnames(mcl)) "SAMPLE_ID" else "PATIENT_ID"
rownames(mcl) <- mcl[[idcol]]
common <- intersect(colnames(MB), rownames(mcl))
cat("samples with expression and clinical:", length(common), "\n")
MB <- MB[, common, drop=FALSE]; mcl <- mcl[common,]
sub_mb <- mcl$CLAUDIN_SUBTYPE
cat("CLAUDIN_SUBTYPE (METABRIC's own PAM50-based call):\n"); print(table(sub_mb, useNA="ifany"))

## any miRNA measured on this platform?
cat("miRNA-like rows in METABRIC matrix:", sum(grepl("^MIR|^hsa-", rownames(MB))), "\n")
cat("MIR29A/B1/B2/C present:",
    paste(c("MIR29A","MIR29B1","MIR29B2","MIR29C") %in% rownames(MB), collapse=","), "\n")

keep_mb <- !is.na(sub_mb) & sub_mb %in% c("LumA","LumB","Her2","Basal","Normal","claudin-low")
MB2 <- MB[, keep_mb, drop=FALSE]; smb <- factor(sub_mb[keep_mb],
        levels=c("LumA","LumB","Her2","Basal","claudin-low","Normal"))
levels(smb)[levels(smb)=="Normal"] <- "Normal-like"
cat("analysed METABRIC n =", ncol(MB2), "\n"); print(table(smb))

cafA_mb <- mkscore(sigA, MB2); cafB_mb <- mkscore(sigB, MB2)
cat("CAF-A genes found in METABRIC:", length(intersect(sigA, rownames(MB2))),
    " CAF-B:", length(intersect(sigB, rownames(MB2))), "\n")

getmb <- function(nm) if(nm %in% rownames(MB2)) MB2[nm,] else NULL
mrows <- list()
mb_axes <- c(paste0(TFS," -> COL1A1"), paste0(TFS," -> COL3A1"), "COL1A1 ~ COL3A1")
for (a in mb_axes){
  parts <- strsplit(a," ~ | -> ")[[1]]
  xs <- getmb(trimws(parts[1])); ys <- getmb(trimws(parts[2]))
  if(is.null(xs)||is.null(ys)){ cat("METABRIC missing feature for", a, "\n"); next }
  for (s in c(levels(smb),"ALL")){
    i <- if(s=="ALL") seq_len(ncol(MB2)) else which(smb==s); n<-length(i); if(n<10) next
    e <- spp(xs[i],ys[i])
    mrows[[length(mrows)+1]] <- data.frame(cohort="METABRIC", stratum=s, n=n, axis=a,
      source=trimws(parts[1]), target=trimws(parts[2]),
      rho=unname(e[1]), p=unname(e[2]),
      rho_partial_cafA=pcor(xs[i],ys[i],cafA_mb[i]),
      rho_partial_cafB=pcor(xs[i],ys[i],cafB_mb[i]), stringsAsFactors=FALSE)
  }
}
MBT <- do.call(rbind, mrows)
cat("\n---------------- METABRIC per-subtype ----------------\n")
for (a in unique(MBT$axis)){
  d <- MBT[MBT$axis==a,]; cat(sprintf("%-22s ", a))
  for (s in c(levels(smb),"ALL")){ r<-d[d$stratum==s,]
    if(nrow(r)) cat(sprintf("%s n=%d rho=%+.3f  ", s, r$n, r$rho)) }
  cat("\n")
}
mhet <- list()
for (a in unique(MBT$axis)){
  d <- MBT[MBT$axis==a & MBT$stratum!="ALL",]
  q <- hetQ(d$rho, d$n)
  parts <- strsplit(a," ~ | -> ")[[1]]
  xs <- getmb(trimws(parts[1])); ys <- getmb(trimws(parts[2]))
  ib <- which(smb=="Basal"); il <- which(smb %in% c("LumA","LumB"))
  rb <- sp(xs[ib],ys[ib]); rl <- sp(xs[il],ys[il])
  zl <- zdiff(rb,length(ib),rl,length(il))
  mhet[[length(mhet)+1]] <- data.frame(cohort="METABRIC", axis=a, k_subtypes=unname(q["k"]),
    Q=unname(q["Q"]), df=unname(q["df"]), p_heterogeneity=unname(q["p"]), I2_pct=unname(q["I2"]),
    rho_Basal=rb, n_Basal=length(ib), rho_nonBasal=sp(xs[-ib],ys[-ib]), n_nonBasal=ncol(MB2)-length(ib),
    z_Basal_vs_nonBasal=NA, p_Basal_vs_nonBasal=NA,
    rho_Luminal=rl, n_Luminal=length(il),
    z_Basal_vs_Luminal=unname(zl["z"]), p_Basal_vs_Luminal=unname(zl["p"]), stringsAsFactors=FALSE)
}
MHET <- do.call(rbind, mhet)
zd2 <- t(mapply(function(rb,nb,ro,no) zdiff(rb,nb,ro,no),
                MHET$rho_Basal, MHET$n_Basal, MHET$rho_nonBasal, MHET$n_nonBasal))
MHET$z_Basal_vs_nonBasal <- zd2[,1]; MHET$p_Basal_vs_nonBasal <- zd2[,2]
MHET$p_het_BH <- p.adjust(MHET$p_heterogeneity,"BH"); MHET$p_BvsNB_BH <- p.adjust(MHET$p_Basal_vs_nonBasal,"BH")
cat("\n---------------- METABRIC heterogeneity ----------------\n")
print(format(MHET[,c("axis","k_subtypes","Q","df","p_heterogeneity","I2_pct","rho_Basal","rho_Luminal","p_Basal_vs_Luminal")], digits=3))

mcaf <- list()
for (s in c(levels(smb),"ALL")){
  i <- if(s=="ALL") seq_len(ncol(MB2)) else which(smb==s)
  mcaf[[length(mcaf)+1]] <- data.frame(cohort="METABRIC", stratum=s, n=length(i),
    cafA_mean=mean(cafA_mb[i]), cafA_sd=sd(cafA_mb[i]), cafA_median=median(cafA_mb[i]),
    cafA_q25=quantile(cafA_mb[i],.25), cafA_q75=quantile(cafA_mb[i],.75),
    cafB_mean=mean(cafB_mb[i]), cafB_sd=sd(cafB_mb[i]), cafB_median=median(cafB_mb[i]),
    COL1A1_mean=if("COL1A1"%in%rownames(MB2)) mean(MB2["COL1A1",i]) else NA,
    COL3A1_mean=if("COL3A1"%in%rownames(MB2)) mean(MB2["COL3A1",i]) else NA,
    stringsAsFactors=FALSE)
}
MCAFD <- do.call(rbind, mcaf); rownames(MCAFD)<-NULL
cat("\n---------------- METABRIC CAF by subtype ----------------\n"); print(format(MCAFD, digits=3))
kwmb <- kruskal.test(cafA_mb ~ smb)
cat(sprintf("Kruskal-Wallis CAF-A across METABRIC subtypes: chi2=%.3f df=%d p=%.3g\n",
            kwmb$statistic, kwmb$parameter, kwmb$p.value))

## ------------------------------------------------------------------ save --
ALL <- rbind(ST, MBT)
ALL$analysis <- "correlation"
CAFALL <- rbind(CAFD, MCAFD)
write.csv(ALL,    file.path(OUT,"subtype_stratified.csv"), row.names=FALSE)
write.csv(rbind(HET,MHET), file.path(OUT,"subtype_heterogeneity.csv"), row.names=FALSE)
write.csv(CAFALL, file.path(OUT,"subtype_caf_distribution.csv"), row.names=FALSE)
cat("\nWROTE:", file.path(OUT,"subtype_stratified.csv"), nrow(ALL), "rows\n")
cat("WROTE:", file.path(OUT,"subtype_heterogeneity.csv"), nrow(HET)+nrow(MHET), "rows\n")
cat("WROTE:", file.path(OUT,"subtype_caf_distribution.csv"), nrow(CAFALL), "rows\n")
cat("DONE PART A\n")
