#!/usr/bin/env Rscript
## ============================================================================
## STROMA-FREE TEST -- part B & C
##   B) TF->collagen and miR-29->collagen correlations in breast cancer lines
##   C) absolute collagen expression: cell lines vs bulk tumours vs fibroblasts
## Data: DepMap Public 24Q4 (figshare 27993248; md5 verified) + CCLE miRNA
##       Nanostring 2018-11-03 + TCGA-BRCA (UCSC Xena) from this project.
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
options(warn=1); set.seed(20260908)
REV <- "/path/to/revision"
DM  <- file.path(REV,"data/depmap24q4"); OUT <- file.path(REV,"results/v2")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

sp <- function(x,y){                       # spearman with n and p
  ok <- is.finite(x)&is.finite(y); n <- sum(ok)
  if(n<8) return(list(rho=NA_real_,p=NA_real_,n=n))
  ct <- suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  list(rho=unname(ct$estimate), p=ct$p.value, n=n)
}

## ---------------- load ------------------------------------------------------
mod <- fread(file.path(DM,"Model.csv"))
L   <- readRDS(file.path(DM,"expr_24Q4.rds")); E <- L$E; gsym <- L$genes
msg("expression:", nrow(E),"models x",ncol(E),"genes")

breast_all <- intersect(mod[OncotreeLineage=="Breast", ModelID], rownames(E))
breast_ca  <- intersect(mod[OncotreeLineage=="Breast" & ModelType=="Cell Line" &
                            OncotreePrimaryDisease!="Non-Cancerous", ModelID], rownames(E))
fibro      <- intersect(mod[OncotreeLineage=="Fibroblast", ModelID], rownames(E))
msg("breast models with expression:", length(breast_all),
    "| malignant breast CELL LINES:", length(breast_ca),
    "| fibroblast-lineage models:", length(fibro))

TFs   <- c("ETS1","NFKB1","RELA","SP1")
COLS  <- c("COL1A1","COL3A1")

## ---------------- C) absolute expression -----------------------------------
tpm <- function(v) 2^v - 1
rank_pct <- function(mat) t(apply(mat,1,function(r) rank(r,na.last="keep")/sum(is.finite(r))*100))

## TCGA bulk
X   <- readRDS(file.path(REV,"data/brca_gene_expr.rds"))
ph  <- readRDS(file.path(REV,"data/brca_pheno.rds"))
tum <- ph$sample[ph$sample_type=="Primary Tumor"]; tum <- intersect(tum, colnames(X))
nor <- ph$sample[ph$sample_type=="Solid Tissue Normal"]; nor <- intersect(nor, colnames(X))
msg("TCGA-BRCA: primary tumours", length(tum), "| solid-tissue normals", length(nor))
common <- intersect(rownames(X), colnames(E))
msg("common gene universe CCLE n TCGA:", length(common))

pct_in <- function(mat_genes_by_sample, genes){          # percentile within sample
  sub <- mat_genes_by_sample[common,,drop=FALSE]
  R <- apply(sub, 2, function(v) rank(v,na.last="keep")/sum(is.finite(v))*100)
  R[genes,,drop=FALSE]
}
tcga_pct   <- pct_in(X[,tum,drop=FALSE], COLS)
tcgaN_pct  <- pct_in(X[,nor,drop=FALSE], COLS)
ccle_pct   <- pct_in(t(E[breast_ca,common,drop=FALSE]), COLS)
fib_pct    <- pct_in(t(E[fibro,common,drop=FALSE]), COLS)

absr <- rbindlist(lapply(COLS, function(g){
  bl <- E[breast_ca,g]; fb <- E[fibro,g]
  data.table(gene=g,
    n_breast_lines=length(bl),
    ccle_median_log2TPMp1=median(bl), ccle_mean_log2TPMp1=mean(bl),
    ccle_median_TPM=median(tpm(bl)), ccle_q25_TPM=quantile(tpm(bl),.25),
    ccle_q75_TPM=quantile(tpm(bl),.75), ccle_max_TPM=max(tpm(bl)),
    ccle_frac_TPM_lt1=mean(tpm(bl)<1), ccle_frac_TPM_lt10=mean(tpm(bl)<10),
    ccle_median_pctile=median(ccle_pct[g,]),
    n_fibroblast_lines=length(fb),
    fibro_median_TPM=median(tpm(fb)), fibro_median_pctile=median(fib_pct[g,]),
    fibro_over_breast_TPM_ratio=median(tpm(fb))/pmax(median(tpm(bl)),1e-6),
    tcga_tumour_median_log2=median(X[g,tum]), tcga_tumour_median_pctile=median(tcga_pct[g,]),
    tcga_normal_median_pctile=median(tcgaN_pct[g,]))
}))
fwrite(absr, file.path(OUT,"celllines_absolute_expression.csv"))
print(absr[, .(gene,ccle_median_TPM,ccle_frac_TPM_lt1,ccle_median_pctile,
               fibro_median_TPM,fibro_median_pctile,tcga_tumour_median_pctile)])

## reference genes for scale
refg <- c("ACTB","GAPDH","EPCAM","KRT8","KRT18","CDH1","VIM","DCN","LUM","FAP","THY1",
          "PDGFRB","ACTA2","COL1A2","COL5A1","COL6A3","FN1","SPARC","LOX")
refg <- intersect(refg, common)
refp_ccle <- pct_in(t(E[breast_ca,common,drop=FALSE]), refg)
refp_tcga <- pct_in(X[,tum,drop=FALSE], refg)
refp_fib  <- pct_in(t(E[fibro,common,drop=FALSE]), refg)
ref <- data.table(gene=refg,
  ccle_breast_median_TPM=round(apply(2^E[breast_ca,refg,drop=FALSE]-1,2,median),2),
  ccle_breast_median_pctile=round(apply(refp_ccle,1,median),1),
  fibroblast_median_TPM=round(apply(2^E[fibro,refg,drop=FALSE]-1,2,median),2),
  fibroblast_median_pctile=round(apply(refp_fib,1,median),1),
  tcga_tumour_median_pctile=round(apply(refp_tcga,1,median),1))
fwrite(ref, file.path(OUT,"celllines_reference_gene_context.csv"))
print(ref)

## ---------------- B) TF -> collagen in cell lines ---------------------------
cell_pairs <- CJ(reg=c(TFs), tgt=COLS, sorted=FALSE)
run_set <- function(ids, label){
  rbindlist(lapply(seq_len(nrow(cell_pairs)), function(i){
    r <- cell_pairs$reg[i]; t <- cell_pairs$tgt[i]
    s <- sp(E[ids,r], E[ids,t])
    data.table(set=label, regulator=r, target=t, n=s$n, rho=s$rho, p=s$p)
  }))
}
tfres <- rbind(run_set(breast_ca,"breast_malignant_cell_lines"),
               run_set(breast_all,"all_breast_models"))
## COL1A1 ~ COL3A1
cc <- sp(E[breast_ca,"COL1A1"], E[breast_ca,"COL3A1"])
tfres <- rbind(tfres, data.table(set="breast_malignant_cell_lines",regulator="COL1A1",
                                 target="COL3A1",n=cc$n,rho=cc$rho,p=cc$p))
cc2 <- sp(E[breast_all,"COL1A1"], E[breast_all,"COL3A1"])
tfres <- rbind(tfres, data.table(set="all_breast_models",regulator="COL1A1",
                                 target="COL3A1",n=cc2$n,rho=cc2$rho,p=cc2$p))
## fibroblast-lineage positive control
tfres <- rbind(tfres, rbindlist(lapply(seq_len(nrow(cell_pairs)), function(i){
  s <- sp(E[fibro,cell_pairs$reg[i]], E[fibro,cell_pairs$tgt[i]])
  data.table(set="fibroblast_lines",regulator=cell_pairs$reg[i],target=cell_pairs$tgt[i],
             n=s$n,rho=s$rho,p=s$p)})))
ccf <- sp(E[fibro,"COL1A1"], E[fibro,"COL3A1"])
tfres <- rbind(tfres, data.table(set="fibroblast_lines",regulator="COL1A1",target="COL3A1",
                                 n=ccf$n,rho=ccf$rho,p=ccf$p))

## TCGA bulk re-computed here for a like-for-like column
tc <- rbindlist(lapply(seq_len(nrow(cell_pairs)), function(i){
  s <- sp(X[cell_pairs$reg[i],tum], X[cell_pairs$tgt[i],tum])
  data.table(set="TCGA_BRCA_bulk_tumour",regulator=cell_pairs$reg[i],
             target=cell_pairs$tgt[i],n=s$n,rho=s$rho,p=s$p)}))
cct <- sp(X["COL1A1",tum], X["COL3A1",tum])
tc <- rbind(tc, data.table(set="TCGA_BRCA_bulk_tumour",regulator="COL1A1",target="COL3A1",
                           n=cct$n,rho=cct$rho,p=cct$p))
tfres <- rbind(tfres, tc)
tfres[, q := p.adjust(p, method="BH"), by=set]
fwrite(tfres, file.path(OUT,"celllines_TF_collagen_correlations.csv"))
print(tfres[set %in% c("breast_malignant_cell_lines","TCGA_BRCA_bulk_tumour")][order(regulator,target,set)])

saveRDS(list(E=E, breast_ca=breast_ca, breast_all=breast_all, fibro=fibro,
             common=common, tum=tum), file.path(DM,"cl_workspace.rds"))
msg("part B/C done")
