## ==========================================================================
## M5_analysis_A.R
##  A) DE-direction concordance with TCGA for the network hubs (mRNA + miRNA)
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(survival); library(metafor)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
S  <- readRDS(file.path(CA,"genesets.rds"))
MI <- readRDS(file.path(CA,"mirna_cohorts.rds"))
set.seed(1)

################################################################################
## A1 -- mRNA hub DE-direction concordance vs TCGA (re-derived, not copied)
################################################################################
cat("\n=========== A1  mRNA hub DE direction concordance ===========\n")
DEt <- fread("results/BRCA_DEX_genes.csv")
cat("TCGA DE table cols:", paste(names(DEt), collapse=","), "\n")
fc <- names(DEt)[grepl("^logFC$", names(DEt))][1]
gn <- names(DEt)[1]
cat("TCGA DE rows:", nrow(DEt), " feature column:", gn,
    " looks-like-symbol fraction:", round(mean(grepl("^[A-Za-z]", DEt[[gn]])),4), "\n")
stopifnot(mean(grepl("^[A-Za-z]", DEt[[gn]])) > 0.95)   ## Rule 4 guard
tcga_fc <- setNames(DEt[[fc]], DEt[[gn]])
qcol <- names(DEt)[grepl("adj.P.Val|^FDR$|^q", names(DEt))][1]
tcga_q  <- setNames(DEt[[qcol]], DEt[[gn]])

EG <- fread("results/external_GEO_DE.csv")
A1 <- list()
for(g in unique(EG$gse)){
  e <- EG[gse==g]
  hub <- intersect(S$HUB_PROTEIN, e$feature)
  hub <- hub[hub %in% names(tcga_fc)]
  hub <- hub[is.finite(tcga_q[hub]) & tcga_q[hub] < 0.05]
  ef <- setNames(e$logFC, e$feature)[hub]
  conc <- sum(sign(ef)==sign(tcga_fc[hub]), na.rm=TRUE)
  n <- sum(is.finite(ef))
  bt <- binom.test(conc, n, 0.5)
  ## all network protein-coding nodes as a wider reference
  nod <- intersect(S$ALL_NETWORK_PROTEIN, e$feature); nod <- nod[nod %in% names(tcga_fc)]
  nod <- nod[is.finite(tcga_q[nod]) & tcga_q[nod]<0.05]
  ef2 <- setNames(e$logFC, e$feature)[nod]
  c2 <- sum(sign(ef2)==sign(tcga_fc[nod]), na.rm=TRUE); n2 <- sum(is.finite(ef2))
  A1[[g]] <- data.table(cohort=g, layer="mRNA", feature_set=c("network_hubs","all_network_nodes"),
    n_tested=c(n,n2), n_concordant=c(conc,c2), pct=round(100*c(conc/n,c2/n2),2),
    binom_p=c(bt$p.value, binom.test(c2,n2,0.5)$p.value),
    n_normal=e$n_normal[1], n_tumour=e$n_tumour[1],
    spearman_logFC=c(cor(ef, tcga_fc[hub], method="spearman", use="complete.obs"),
                     cor(ef2, tcga_fc[nod], method="spearman", use="complete.obs")))
  cat(sprintf("  %-10s hubs %3d/%3d (%.1f%%) p=%.3g | all nodes %3d/%3d (%.1f%%) p=%.3g | rho=%.3f\n",
      g, conc, n, 100*conc/n, bt$p.value, c2, n2, 100*c2/n2,
      binom.test(c2,n2,0.5)$p.value, cor(ef2, tcga_fc[nod], method="spearman", use="complete.obs")))
}
A1 <- rbindlist(A1)

################################################################################
## A2 -- miRNA hub DE-direction concordance vs TCGA (NEW)
################################################################################
cat("\n=========== A2  miRNA DE direction concordance ===========\n")
DEm <- fread("results/BRCA_DEX_mirnas.csv")
cat("TCGA miRNA DE cols:", paste(names(DEm), collapse=","), " rows:", nrow(DEm), "\n")
mfc <- setNames(DEm[[grep("^logFC$", names(DEm))[1]]], DEm[[1]])
mq  <- setNames(DEm[[grep("adj.P.Val|^FDR$|^q", names(DEm))[1]]], DEm[[1]])
strip <- function(x) sub("-3p$|-5p$","", x)
A2 <- list()
for(nm in names(MI)){
  o <- MI[[nm]]
  if(is.null(o$normals) || ncol(o$normals)==0) next
  M <- o$M; N <- o$normals
  common <- intersect(rownames(M), rownames(N))
  M <- M[common,,drop=FALSE]; N <- N[common,,drop=FALSE]
  lfc <- rowMeans(M, na.rm=TRUE) - rowMeans(N, na.rm=TRUE)
  pv  <- apply(cbind(M,N), 1, function(r){
    g <- c(rep(1,ncol(M)), rep(0,ncol(N)))
    if(sum(is.finite(r[g==1]))<3 || sum(is.finite(r[g==0]))<3) return(NA_real_)
    tryCatch(wilcox.test(r[g==1], r[g==0])$p.value, error=function(e) NA_real_)})
  qv <- p.adjust(pv,"BH")
  ## match to TCGA miRNA names, allowing the arm suffix to be absent on arrays
  key <- names(mfc); names(key) <- strip(key)
  map <- ifelse(names(lfc) %in% names(mfc), names(lfc), key[names(lfc)])
  ok <- !is.na(map) & is.finite(mfc[map]) & is.finite(lfc)
  cat(sprintf("  %-10s miRNAs measured %4d, matched to the TCGA miRNA DE table %4d\n",
      nm, length(lfc), sum(ok)))
  for(lab in c("all_matched_miRNAs","miRNA_hubs")){
    sel <- ok
    if(lab=="miRNA_hubs") sel <- ok & (map %in% S$HUB_MIRNA | strip(map) %in% strip(S$HUB_MIRNA))
    sel <- sel & is.finite(mq[map]) & mq[map] < 0.05
    n <- sum(sel); if(n<3) next
    conc <- sum(sign(lfc[sel])==sign(mfc[map[sel]]))
    bt <- binom.test(conc,n,0.5)
    A2[[length(A2)+1]] <- data.table(cohort=nm, layer="miRNA", feature_set=lab,
      n_tested=n, n_concordant=conc, pct=round(100*conc/n,2), binom_p=bt$p.value,
      n_normal=ncol(N), n_tumour=ncol(M),
      spearman_logFC=cor(lfc[sel], mfc[map[sel]], method="spearman"))
    cat(sprintf("    %-20s %3d/%3d (%.1f%%) p=%.3g rho=%.3f\n", lab, conc, n,
        100*conc/n, bt$p.value, cor(lfc[sel], mfc[map[sel]], method="spearman")))
  }
}
A2 <- rbindlist(A2)
A <- rbind(A1, A2, fill=TRUE)
fwrite(A, file.path(OUT,"multicohort_A_DE_direction_concordance.csv"))
