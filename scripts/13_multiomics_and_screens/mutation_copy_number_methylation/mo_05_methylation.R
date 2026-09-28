#!/usr/bin/env Rscript
## Part A: TCGA-BRCA Illumina 450k promoter methylation of network nodes
suppressPackageStartupMessages({library(data.table)})
RV<-"/path/to/revision"; CA<-file.path(RV,"cache/multiomics")
OUT<-file.path(RV,"results/multiomics"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

M <- fread(cmd=paste0("zcat ",CA,"/meth_subset.tsv.gz"))
setnames(M,1,"probe")
bm <- as.matrix(M[,-1]); rownames(bm) <- M$probe
cat("beta matrix:", nrow(bm),"probes x",ncol(bm),"samples\n")
nna <- rowSums(!is.na(bm))
cat("probes with 0 non-NA values (masked):", sum(nna==0),"\n")
bm <- bm[nna >= 100, , drop=FALSE]
cat("probes retained (>=100 non-NA):", nrow(bm),"\n")
sfx <- sub(".*-","",colnames(bm))
tum <- colnames(bm)[sfx=="01"]; nor <- colnames(bm)[sfx=="11"]
cat("tumour samples:",length(tum)," normal samples:",length(nor),
    " (excluded metastatic:",sum(sfx=="06"),")\n")

nodes <- fread(file.path(RV,"data/canonical_nodes.tsv"))
pg <- fread(file.path(CA,"promoter_probes_network_genes.tsv"))
pb <- fread(file.path(CA,"promoter_probes_background_genes.tsv"))
pm <- fread(file.path(CA,"promoter_probes_mirnas.tsv"))

ge <- readRDS(file.path(RV,"data/brca_gene_expr.rds"))
me <- readRDS(file.path(RV,"data/brca_mirna_expr_canonical.rds"))
dex <- fread(file.path(RV,"results/BRCA_DEX_ALL_nodes.csv"))
setnames(dex, c("Gene","logFC","adj.P.Val","class"), c("node","logFC","de_fdr","de_class"), skip_absent=TRUE)

analyse <- function(name, probes, type, exprmat){
  probes <- intersect(unique(probes), rownames(bm))
  if(!length(probes)) return(NULL)
  sub <- bm[probes,,drop=FALSE]
  ps  <- colMeans(sub, na.rm=TRUE)                # promoter score per sample
  t_v <- ps[tum]; n_v <- ps[nor]
  t_v <- t_v[is.finite(t_v)]; n_v <- n_v[is.finite(n_v)]
  if(length(t_v)<30 || length(n_v)<10) return(NULL)
  wt <- suppressWarnings(wilcox.test(t_v, n_v))
  rho<-NA_real_; pv<-NA_real_; nn<-NA_integer_
  if(!is.null(exprmat) && name %in% rownames(exprmat)){
    ex <- exprmat[name,]; cs <- intersect(names(t_v), names(ex))
    if(length(cs)>=30){ ct <- suppressWarnings(cor.test(ps[cs], ex[cs], method="spearman"))
      rho<-unname(ct$estimate); pv<-ct$p.value; nn<-length(cs) }
  }
  data.table(node=name, type=type, n_probes=length(probes),
             beta_tumour=mean(t_v), beta_normal=mean(n_v),
             delta_beta=mean(t_v)-mean(n_v), wilcox_p=wt$p.value,
             n_tumour=length(t_v), n_normal=length(n_v),
             meth_expr_rho=rho, meth_expr_p=pv, n_meth_expr=nn)
}

## network protein-coding
r1 <- rbindlist(lapply(unique(pg$gene), function(g)
        analyse(g, pg[gene==g, Name], nodes$type[match(g,nodes$name)], ge)))
## network miRNAs
r2 <- rbindlist(lapply(unique(pm$mature), function(m)
        analyse(m, pm[mature==m, Name], "miRNA", me)))
## background genes
r3 <- rbindlist(lapply(unique(pb$gene), function(g) analyse(g, pb[gene==g, Name], "background", ge)))

res <- rbind(r1,r2)
res[, wilcox_fdr := p.adjust(wilcox_p,"BH")]
res[!is.na(meth_expr_p), meth_expr_fdr := p.adjust(meth_expr_p,"BH")]
res <- merge(res, dex[,.(node,logFC,de_fdr,de_class)], by="node", all.x=TRUE)
fwrite(res[order(-delta_beta)], file.path(OUT,"brca_methylation_nodes.csv"))
cat("\nwrote brca_methylation_nodes.csv rows:", nrow(res),"\n")
print(res[, .N, by=type])
r3[, wilcox_fdr := p.adjust(wilcox_p,"BH")]
r3[!is.na(meth_expr_p), meth_expr_fdr := p.adjust(meth_expr_p,"BH")]
fwrite(r3, file.path(OUT,"brca_methylation_background_genes.csv"))
cat("wrote brca_methylation_background_genes.csv rows:", nrow(r3),"\n")

cat("\n=== differential promoter methylation, tumour vs normal ===\n")
for(t in c("Gene","TF","miRNA")){ s<-res[type==t]
  cat(sprintf("%-6s n=%3d | median delta=%+.4f | hyper(FDR<.05 & d>0.05)=%d | hypo(FDR<.05 & d<-0.05)=%d\n",
      t,nrow(s),median(s$delta_beta), sum(s$wilcox_fdr<0.05 & s$delta_beta>0.05),
      sum(s$wilcox_fdr<0.05 & s$delta_beta< -0.05))) }
cat(sprintf("%-6s n=%3d | median delta=%+.4f | hyper=%d | hypo=%d\n","bkgd",nrow(r3),
    median(r3$delta_beta), sum(r3$wilcox_fdr<0.05 & r3$delta_beta>0.05),
    sum(r3$wilcox_fdr<0.05 & r3$delta_beta< -0.05)))

cat("\n=== KEY TEST: promoter methylation change vs expression logFC ===\n")
kt <- function(d,lab){ d<-d[!is.na(logFC) & !is.na(delta_beta)]
  if(nrow(d)<10){cat(lab,": too few\n"); return(invisible())}
  ct<-suppressWarnings(cor.test(d$delta_beta,d$logFC,method="spearman"))
  cp<-suppressWarnings(cor.test(d$delta_beta,d$logFC,method="pearson"))
  cat(sprintf("%-28s n=%3d  Spearman rho=%+.3f p=%-10.3g  Pearson r=%+.3f p=%.3g\n",
      lab,nrow(d),unname(ct$estimate),ct$p.value,unname(cp$estimate),cp$p.value)) }
kt(res,"ALL network nodes"); kt(res[type=="Gene"],"network Genes")
kt(res[type=="TF"],"network TFs"); kt(res[type=="miRNA"],"network miRNAs")
bgm <- merge(r3, dex[,.(node,logFC)], by="node"); kt(bgm,"background genes")

cat("\n=== within-tumour promoter-methylation vs expression correlation ===\n")
for(t in c("Gene","TF","miRNA")){ s<-res[type==t & !is.na(meth_expr_rho)]
  cat(sprintf("%-6s n=%3d median rho=%+.3f  %%negative=%.1f  %%FDR<.05 & rho<0=%.1f\n",
      t,nrow(s),median(s$meth_expr_rho),100*mean(s$meth_expr_rho<0),
      100*mean(s$meth_expr_fdr<0.05 & s$meth_expr_rho<0))) }
s<-r3[!is.na(meth_expr_rho)]
cat(sprintf("%-6s n=%3d median rho=%+.3f  %%negative=%.1f  %%FDR<.05 & rho<0=%.1f\n","bkgd",
    nrow(s),median(s$meth_expr_rho),100*mean(s$meth_expr_rho<0),
    100*mean(s$meth_expr_fdr<0.05 & s$meth_expr_rho<0)))

cat("\n=== FEATURED miRNA PROMOTERS ===\n")
feat <- c("hsa-miR-130a","hsa-miR-124","hsa-miR-101","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c",
          "hsa-let-7b","hsa-let-7e","hsa-miR-34a","hsa-miR-200b","hsa-miR-200c","hsa-miR-145")
print(res[node %in% feat, .(node,n_probes,beta_normal=round(beta_normal,3),
   beta_tumour=round(beta_tumour,3), delta=round(delta_beta,4), fdr=signif(wilcox_fdr,2),
   logFC=round(logFC,3), de_fdr=signif(de_fdr,2), meth_expr_rho=round(meth_expr_rho,3),
   me_fdr=signif(meth_expr_fdr,2))][match(feat,node)])

cat("\n=== featured genes/TFs ===\n")
fg <- c("COL1A1","COL3A1","EZH2","RELA","NFKB1","SP1","TP53","MYC","ESR1")
print(res[node %in% fg, .(node,type,n_probes,beta_normal=round(beta_normal,3),
   beta_tumour=round(beta_tumour,3), delta=round(delta_beta,4), fdr=signif(wilcox_fdr,2),
   logFC=round(logFC,3), meth_expr_rho=round(meth_expr_rho,3), me_fdr=signif(meth_expr_fdr,2))])

cat("\n=== most hypermethylated network nodes (FDR<0.05) ===\n")
print(res[wilcox_fdr<0.05][order(-delta_beta)][1:15,.(node,type,delta=round(delta_beta,3),
   beta_normal=round(beta_normal,3),beta_tumour=round(beta_tumour,3),logFC=round(logFC,2),
   meth_expr_rho=round(meth_expr_rho,2))])
cat("\n=== most hypomethylated network nodes (FDR<0.05) ===\n")
print(res[wilcox_fdr<0.05][order(delta_beta)][1:15,.(node,type,delta=round(delta_beta,3),
   beta_normal=round(beta_normal,3),beta_tumour=round(beta_tumour,3),logFC=round(logFC,2),
   meth_expr_rho=round(meth_expr_rho,2))])

## per-probe detail for featured miRNAs
det <- rbindlist(lapply(feat, function(m){
  pr <- intersect(pm[mature==m, Name], rownames(bm)); if(!length(pr)) return(NULL)
  data.table(mature=m, probe=pr,
    beta_normal=rowMeans(bm[pr,nor,drop=FALSE],na.rm=TRUE),
    beta_tumour=rowMeans(bm[pr,tum,drop=FALSE],na.rm=TRUE))}))
det[, delta := beta_tumour-beta_normal]
det <- merge(det, unique(pm[,.(mature,Name,src,precursor)]), by.x=c("mature","probe"),
             by.y=c("mature","Name"), all.x=TRUE)
fwrite(det[order(mature,-delta)], file.path(OUT,"brca_methylation_featured_mirna_probes.csv"))
cat("\nwrote brca_methylation_featured_mirna_probes.csv rows:", nrow(det),"\n")
saveRDS(list(res=res,bg=r3), file.path(CA,"meth_results.rds"))
