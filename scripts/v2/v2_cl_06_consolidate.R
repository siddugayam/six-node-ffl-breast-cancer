#!/usr/bin/env Rscript
## ============================================================================
## STROMA-FREE TEST -- consolidation, EMT adjustment, summary table
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
options(warn=1); set.seed(20260908)
REV <- "/path/to/revision"
DM  <- file.path(REV,"data/depmap24q4"); OUT <- file.path(REV,"results/v2")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
sp <- function(x,y){ ok<-is.finite(x)&is.finite(y); n<-sum(ok)
  ct<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  c(rho=unname(ct$estimate),p=ct$p.value,n=n) }
pcor <- function(x,y,Z){ Z<-as.matrix(Z); ok<-is.finite(x)&is.finite(y)&complete.cases(Z)
  xr<-rank(x[ok]); yr<-rank(y[ok]); Zr<-apply(Z[ok,,drop=FALSE],2,rank)
  rx<-resid(lm(xr~Zr)); ry<-resid(lm(yr~Zr)); n<-sum(ok); r<-cor(rx,ry)
  df<-n-2-ncol(Z); t<-r*sqrt(df/(1-r^2)); c(rho=r,p=2*pt(-abs(t),df),n=n) }

W <- readRDS(file.path(DM,"cl_workspace.rds")); E<-W$E
mod<- fread(file.path(DM,"Model.csv"))
ce <- fread("/path/to/home/Desktop/DD/R_GPR/ESIA_BRCA/brca_eisa_pilot/depmap/CRISPRGeneEffect.csv",
            select=1:2); setnames(ce,1,"ModelID")
br_all_ce <- intersect(mod[OncotreeLineage=="Breast" & ModelType=="Cell Line" &
                           OncotreePrimaryDisease!="Non-Cancerous", ModelID], ce$ModelID)
msg("malignant breast CELL LINES with a CRISPR screen (24Q4):", length(br_all_ce))
br_any_ce <- intersect(mod[OncotreeLineage=="Breast", ModelID], ce$ModelID)
msg("all breast models with a CRISPR screen (24Q4):", length(br_any_ce))

ids <- W$breast_ca
det <- function(g,cut) sum(2^E[ids,g]-1 >= cut)
msg("COL1A1 detectable (TPM>=1) in", det("COL1A1",1), "of", length(ids), "breast lines;",
    "TPM>=10 in", det("COL1A1",10))
msg("COL3A1 detectable (TPM>=1) in", det("COL3A1",1), "of", length(ids), "breast lines;",
    "TPM>=10 in", det("COL3A1",10))

## ---- EMT adjustment of the TF arm inside cell lines -----------------------
emtg <- intersect(c("VIM","CDH2","ZEB1","ZEB2","SNAI2","TWIST1","FN1","SERPINE1"), colnames(E))
emt  <- rowMeans(scale(E[ids, emtg, drop=FALSE]))
msg("cell-intrinsic EMT score genes:", paste(emtg, collapse=", "))
pairs <- CJ(reg=c("ETS1","NFKB1","RELA","SP1"), tgt=c("COL1A1","COL3A1"), sorted=FALSE)
emtres <- rbindlist(lapply(seq_len(nrow(pairs)), function(i){
  r <- pairs$reg[i]; t <- pairs$tgt[i]
  a <- sp(E[ids,r], E[ids,t]); b <- pcor(E[ids,r], E[ids,t], cbind(emt))
  data.table(regulator=r, target=t, n=a["n"], rho_raw=a["rho"], p_raw=a["p"],
             rho_EMTadj=b["rho"], p_EMTadj=b["p"])}))
## COL1A1<->COL3A1 too
cc <- sp(E[ids,"COL1A1"],E[ids,"COL3A1"]); ccp <- pcor(E[ids,"COL1A1"],E[ids,"COL3A1"],cbind(emt))
emtres <- rbind(emtres, data.table(regulator="COL1A1",target="COL3A1",n=cc["n"],
                 rho_raw=cc["rho"],p_raw=cc["p"],rho_EMTadj=ccp["rho"],p_EMTadj=ccp["p"]))
emtres[, q_raw := p.adjust(p_raw,"BH")][, q_EMTadj := p.adjust(p_EMTadj,"BH")]
fwrite(emtres, file.path(OUT,"celllines_TF_collagen_EMTadjusted.csv")); print(emtres)

## ---- summary verdict table ------------------------------------------------
abs_ <- fread(file.path(OUT,"celllines_absolute_expression.csv"))
cmp  <- fread(file.path(OUT,"celllines_tcga_vs_cellline_named_axes.csv"))
nul  <- fread(file.path(OUT,"celllines_concordance_vs_matched_null.csv"))
S <- rbindlist(list(
 data.table(item="release", value="DepMap Public 24Q4 (figshare 27993248); md5 of OmicsExpressionProteinCodingGenesTPMLogp1.csv = 71794802b750ce77c422dad0720a40af; CRISPRGeneEffect.csv = 6edf7ade09b9b34199210b559d4745d3"),
 data.table(item="breast_models_total", value=as.character(sum(mod$OncotreeLineage=="Breast"))),
 data.table(item="breast_models_with_expression", value=as.character(length(W$breast_all))),
 data.table(item="malignant_breast_cell_lines_with_expression", value=as.character(length(ids))),
 data.table(item="malignant_breast_cell_lines_with_CRISPR", value=as.character(length(br_all_ce))),
 data.table(item="breast_lines_with_CCLE_miRNA_and_expression", value="50"),
 data.table(item="fibroblast_lineage_models_with_expression", value=as.character(length(W$fibro))),
 data.table(item="COL1A1_median_TPM_breast_lines", value=as.character(round(abs_[gene=="COL1A1",ccle_median_TPM],2))),
 data.table(item="COL3A1_median_TPM_breast_lines", value=as.character(round(abs_[gene=="COL3A1",ccle_median_TPM],2))),
 data.table(item="COL1A1_median_TPM_fibroblast_lines", value=as.character(round(abs_[gene=="COL1A1",fibro_median_TPM],1))),
 data.table(item="COL3A1_median_TPM_fibroblast_lines", value=as.character(round(abs_[gene=="COL3A1",fibro_median_TPM],1))),
 data.table(item="COL1A1_fibroblast_over_breastline_ratio", value=as.character(round(abs_[gene=="COL1A1",fibro_over_breast_TPM_ratio],0))),
 data.table(item="COL3A1_fibroblast_over_breastline_ratio", value=as.character(round(abs_[gene=="COL3A1",fibro_over_breast_TPM_ratio],0))),
 data.table(item="COL3A1_frac_breast_lines_below_1TPM", value=as.character(round(abs_[gene=="COL3A1",ccle_frac_TPM_lt1],3))),
 data.table(item="COL1A1_pctile_TCGA_tumour_vs_breastline", value=paste0(round(abs_[gene=="COL1A1",tcga_tumour_median_pctile],2)," vs ",round(abs_[gene=="COL1A1",ccle_median_pctile],2))),
 data.table(item="COL3A1_pctile_TCGA_tumour_vs_breastline", value=paste0(round(abs_[gene=="COL3A1",tcga_tumour_median_pctile],2)," vs ",round(abs_[gene=="COL3A1",ccle_median_pctile],2)))
))
S <- rbind(S, cmp[, .(item=paste0("axis: ",axis),
      value=sprintf("TCGA rho %.3f (n=%d) -> cell lines rho %.3f [95%% CI %.3f, %.3f] (n=%d), p_diff=%.2g",
                    tcga_rho,tcga_n,cl_rho,cl_ci_lo,cl_ci_hi,cl_n,p_diff))])
S <- rbind(S, nul[, .(item=paste0("panel: ",panel),
      value=sprintf("obs %.3f vs matched null %.3f (n_edges=%d, n_null=%d), p=%.3g",
                    obs_frac,null_frac,n_edges,n_null,p))])
fwrite(S, file.path(OUT,"celllines_summary.csv"))
print(S)
msg("done")
