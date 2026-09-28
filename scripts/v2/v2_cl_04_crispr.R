#!/usr/bin/env Rscript
## ============================================================================
## STROMA-FREE TEST -- part D: dependency x collagen-expression relationships
## CRISPRGeneEffect.csv = DepMap Public 24Q4 (md5 6edf7ade09b9b34199210b559d4745d3)
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
options(warn=1); set.seed(20260908)
REV <- "/path/to/revision"
DM  <- file.path(REV,"data/depmap24q4"); OUT <- file.path(REV,"results/v2")
CE_PATH <- "/path/to/home/Desktop/DD/R_GPR/ESIA_BRCA/brca_eisa_pilot/depmap/CRISPRGeneEffect.csv"
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
sp <- function(x,y){ ok<-is.finite(x)&is.finite(y); n<-sum(ok)
  if(n<8) return(c(rho=NA,p=NA,n=n))
  ct<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  c(rho=unname(ct$estimate),p=ct$p.value,n=n) }

W <- readRDS(file.path(DM,"cl_workspace.rds")); E<-W$E
mod <- fread(file.path(DM,"Model.csv"))
msg("reading CRISPRGeneEffect.csv (24Q4) ...")
ce <- fread(CE_PATH); setnames(ce,1,"ModelID")
gc_ <- setdiff(names(ce),"ModelID"); gs <- sub(" \\(\\d+\\)$","",gc_)
stopifnot(!any(duplicated(gs)), all(grepl("^[A-Za-z]",gs)))
C <- as.matrix(ce[,gc_,with=FALSE]); rownames(C)<-ce$ModelID; colnames(C)<-gs
msg("CRISPR:", nrow(C),"lines x",ncol(C),"genes")

br_ce <- intersect(W$breast_ca, rownames(C))
all_ce<- intersect(rownames(E), rownames(C))
msg("malignant breast cell lines with BOTH CRISPR and expression:", length(br_ce))
msg("all cell lines with BOTH CRISPR and expression:", length(all_ce))

## ---- D1: does knocking out a TF co-vary with collagen expression? ---------
HUB <- c("NFKB1","RELA","SP1","ETS1","MYC","TP53","EZH2","STAT3","JUN","FOS","HIF1A")
HUB <- intersect(HUB, colnames(C))
TGT <- intersect(c("COL1A1","COL3A1","COL1A2","COL5A1","FN1","SPARC"), colnames(E))
d1 <- rbindlist(lapply(HUB, function(h) rbindlist(lapply(TGT, function(t){
  a <- sp(C[br_ce,h],  E[br_ce,t])
  b <- sp(C[all_ce,h], E[all_ce,t])
  data.table(hub=h, collagen_gene=t,
             n_breast=a["n"], rho_breast=a["rho"], p_breast=a["p"],
             n_all=b["n"],    rho_all=b["rho"],    p_all=b["p"],
             mean_chronos_breast=mean(C[br_ce,h],na.rm=TRUE),
             n_essential_breast=sum(C[br_ce,h] < -0.5, na.rm=TRUE))}))))
d1[, q_breast := p.adjust(p_breast,"BH")][, q_all := p.adjust(p_all,"BH")]
fwrite(d1, file.path(OUT,"celllines_crispr_hub_vs_collagen_expression.csv"))
print(d1[collagen_gene %in% c("COL1A1","COL3A1")][order(hub,collagen_gene)][,
        .(hub,collagen_gene,n_breast,rho_breast,q_breast,rho_all,q_all,
          mean_chronos_breast,n_essential_breast)])

## ---- D2: collagen-high vs collagen-low breast lines: differential deps ----
ECM <- intersect(c("COL1A1","COL3A1","COL1A2","COL5A1","COL5A2","COL6A3","COL11A1",
                   "FN1","SPARC","LOX"), colnames(E))
msg("ECM/collagen module genes used:", paste(ECM, collapse=", "))
score <- rowMeans(scale(E[br_ce, ECM, drop=FALSE]))
msg("collagen module score: n =", length(score), "| range",
    paste(round(range(score),3),collapse=" .. "))
hi <- names(score)[score >  median(score)]; lo <- names(score)[score <= median(score)]
msg("collagen-high lines:", length(hi), "| collagen-low lines:", length(lo))
writeLines(paste(names(score), round(score,4), ifelse(names(score)%in%hi,"high","low")),
           file.path(OUT,"celllines_collagen_module_score.txt"))

nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
net <- intersect(nodes[type %in% c("TF","Gene"), name], colnames(C))
msg("network protein nodes present in the CRISPR matrix:", length(net))

test_set <- function(genes, label){
  rbindlist(lapply(genes, function(g){
    a <- C[hi,g]; b <- C[lo,g]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if(length(a)<5 || length(b)<5) return(NULL)
    w <- suppressWarnings(wilcox.test(a,b))
    s <- sp(score[br_ce], C[br_ce,g])
    data.table(set=label, gene=g, n_hi=length(a), n_lo=length(b),
               mean_hi=mean(a), mean_lo=mean(b), delta=mean(a)-mean(b),
               wilcox_p=w$p.value, rho_score=s["rho"], p_score=s["p"])
  }))
}
dn <- test_set(net, "network_nodes")
dn[, wilcox_q := p.adjust(wilcox_p,"BH")][, q_score := p.adjust(p_score,"BH")]
setorder(dn, wilcox_p)
fwrite(dn, file.path(OUT,"celllines_crispr_collagen_high_vs_low_networknodes.csv"))
msg("network nodes with FDR<0.10 (collagen-high vs low):", sum(dn$wilcox_q<0.10, na.rm=TRUE))
print(head(dn[,.(gene,n_hi,n_lo,mean_hi,mean_lo,delta,wilcox_p,wilcox_q,rho_score,q_score)],12))

## genome-wide exploratory
gw_ok <- colnames(C)[colSums(is.finite(C[br_ce,,drop=FALSE]))>=length(br_ce)*0.9]
msg("genome-wide exploratory genes tested:", length(gw_ok))
dg <- test_set(gw_ok, "genome_wide")
dg[, wilcox_q := p.adjust(wilcox_p,"BH")][, q_score := p.adjust(p_score,"BH")]
setorder(dg, wilcox_p)
fwrite(head(dg,500), file.path(OUT,"celllines_crispr_collagen_high_vs_low_top500.csv"))
msg("genome-wide genes with FDR<0.10:", sum(dg$wilcox_q<0.10, na.rm=TRUE),
    "| with FDR<0.25:", sum(dg$wilcox_q<0.25, na.rm=TRUE))
print(head(dg[,.(gene,mean_hi,mean_lo,delta,wilcox_p,wilcox_q)],10))
msg("done")
