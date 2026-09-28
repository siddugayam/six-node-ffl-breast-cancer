#!/usr/bin/env Rscript
## Part D: integrated multi-omic summary for the manuscript's featured hubs
suppressPackageStartupMessages({library(data.table)})
RV<-"/path/to/revision"; OUT<-file.path(RV,"results/multiomics")

meth <- fread(file.path(OUT,"brca_methylation_nodes.csv"))
mut  <- fread(file.path(OUT,"brca_mutation_frequency.csv"))
cnv  <- fread(file.path(OUT,"brca_copynumber_nodes.csv"))
cnl  <- fread(file.path(OUT,"brca_copynumber_loci_detail.csv"))
dex  <- fread(file.path(RV,"results/BRCA_DEX_ALL_nodes.csv"))
setnames(dex, 1:3, c("node","logFC","de_fdr"))

## hub definition: manuscript-named hubs (survival_cox_hubs.csv flag) + featured miRNAs in the analysis plan
sc <- fread(file.path(RV,"results/survival_cox_hubs.csv"))
ms_hubs <- unique(sc[hub_named_in_MS==TRUE, .(node=feature, node_type)])
extra <- data.table(node=c("hsa-miR-34a","hsa-miR-200b","hsa-miR-200c","hsa-miR-145","EZH2"),
                    node_type=c("miRNA","miRNA","miRNA","miRNA","TF"))
hubs <- unique(rbind(ms_hubs, extra), by="node")
hubs[, hub_source := ifelse(node %in% ms_hubs$node, "manuscript_named_hub", "manuscript_featured_miRNA/gene")]
cat("featured hubs:", nrow(hubs), "\n")

H <- merge(hubs, dex[,.(node,logFC,de_fdr)], by="node", all.x=TRUE)
H <- merge(H, meth[,.(node,n_meth_probes=n_probes,beta_normal,beta_tumour,delta_beta,
                      meth_fdr=wilcox_fdr,meth_expr_rho,meth_expr_fdr)], by="node", all.x=TRUE)
H <- merge(H, mut[,.(node=gene,mut_n=n_mutated,mut_n_samples=n_samples,mut_freq)], by="node", all.x=TRUE)
H <- merge(H, cnv[,.(node,cn_locus=locus,cn_chrom=chrom,cn_source,cn_method,n_cn_loci=n_loci,
                     frac_amp,frac_del,frac_deepdel,cn_expr_rho,cn_expr_fdr,all_loci)],
           by="node", all.x=TRUE)

## verdict per hub
H[, dir := fifelse(is.na(logFC),"NA", fifelse(de_fdr<0.05 & logFC>0,"UP",
            fifelse(de_fdr<0.05 & logFC<0,"DOWN","not_DE")))]
H[, meth_call := fifelse(is.na(delta_beta),"no_data",
    fifelse(meth_fdr<0.05 & delta_beta>0.05,"promoter_hypermethylated",
    fifelse(meth_fdr<0.05 & delta_beta< -0.05,"promoter_hypomethylated","no_marked_change")))]
H[, cnv_call := fifelse(is.na(frac_amp),"no_data",
    fifelse(frac_amp>=0.30 & frac_amp>frac_del,"recurrently_amplified",
    fifelse(frac_del>=0.30 & frac_del>frac_amp,"recurrently_deleted","no_recurrent_CNV")))]
H[, mut_call := fifelse(is.na(mut_freq),"not_applicable(miRNA)",
    fifelse(mut_freq>=0.05,"recurrently_mutated",
    fifelse(mut_freq>=0.02,"occasionally_mutated","rarely_mutated")))]
H[, basis := {
  b <- character(.N)
  for(i in seq_len(.N)){
    v <- c()
    if(dir[i]=="DOWN" && meth_call[i]=="promoter_hypermethylated") v <- c(v,"epigenetic silencing (consistent)")
    if(dir[i]=="UP"   && meth_call[i]=="promoter_hypomethylated")  v <- c(v,"promoter demethylation (consistent)")
    if(dir[i]=="DOWN" && cnv_call[i]=="recurrently_deleted" && !is.na(cn_expr_fdr[i]) && cn_expr_fdr[i]<0.05 && cn_expr_rho[i]>0) v <- c(v,"copy-number loss (consistent)")
    if(dir[i]=="UP"   && cnv_call[i]=="recurrently_amplified" && !is.na(cn_expr_fdr[i]) && cn_expr_fdr[i]<0.05 && cn_expr_rho[i]>0) v <- c(v,"copy-number gain (consistent)")
    if(mut_call[i]=="recurrently_mutated") v <- c(v,"recurrent somatic mutation")
    b[i] <- if(length(v)) paste(v, collapse="; ") else
      if(dir[i]=="not_DE"||dir[i]=="NA") "not differentially expressed" else "no genomic/epigenomic basis identified"
  }; b }]
H[, cnv_note := fifelse(is.na(frac_amp),"", fifelse(cnv_call!="no_recurrent_CNV" & (is.na(cn_expr_fdr) | cn_expr_fdr>=0.05),
   "CNV present but not coupled to expression (likely arm-level passenger)",""))]
setcolorder(H, c("node","node_type","hub_source","logFC","de_fdr","dir",
  "n_meth_probes","beta_normal","beta_tumour","delta_beta","meth_fdr","meth_call",
  "meth_expr_rho","meth_expr_fdr","mut_n","mut_n_samples","mut_freq","mut_call",
  "cn_locus","cn_chrom","cn_source","cn_method","n_cn_loci","frac_amp","frac_del","frac_deepdel",
  "cn_expr_rho","cn_expr_fdr","cnv_call","cnv_note","basis","all_loci"))
setorder(H, node_type, node)
fwrite(H, file.path(OUT,"hub_multiomic_summary.csv"))
cat("wrote hub_multiomic_summary.csv rows:", nrow(H), "cols:", ncol(H), "\n\n")

print(H[,.(node,type=node_type,logFC=round(logFC,2),dir,
   dBeta=round(delta_beta,3), meth=substr(meth_call,1,26),
   mut=round(mut_freq,3), amp=round(frac_amp,2), del=round(frac_del,2),
   cnv=substr(cnv_call,1,22))], nrows=40)
cat("\n=== BASIS ===\n")
for(i in seq_len(nrow(H))) cat(sprintf("%-14s %-6s %-8s : %s\n", H$node[i], H$node_type[i], H$dir[i], H$basis[i]))
cat("\n=== counts ===\n"); print(H[,.N,by=basis][order(-N)])

## miR-29 cluster loci detail
cat("\n=== miR-29 cluster / featured loci (all loci) ===\n")
print(cnl[node %in% c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-130a","COL1A1","COL3A1"),
  .(node,locus,chrom,cn_source,cn_method,amp=round(frac_amp,3),del=round(frac_del,3),
    rho=round(cn_expr_rho,3),fdr=signif(cn_expr_fdr,2))])
