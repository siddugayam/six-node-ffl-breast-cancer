#!/usr/bin/env Rscript
## Length-matched control for the "network genes are more mutated" claim
suppressPackageStartupMessages({library(data.table)})
RV<-"/path/to/revision"; CA<-file.path(RV,"cache/multiomics")
OUT<-file.path(RV,"results/multiomics")
mf <- fread(file.path(OUT,"brca_mutation_frequency.csv"))
nodes <- fread(file.path(RV,"data/canonical_nodes.tsv"))

## recompute full frequency table (network + background)
ph <- fread(cmd="zcat /path/to/home/Desktop/DD/R_GPR/ML/TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
setnames(ph, make.names(names(ph)))
brca <- ph[grepl("breast",X_primary_disease,ignore.case=TRUE) & sample_type=="Primary Tumor", sample]
mc3 <- fread(cmd="zcat /path/to/home/Desktop/DD/R_GPR/ML/external_cohorts/MC3/mc3_nonsilentGene.xena.gz")
setnames(mc3,1,"gene")
keep <- intersect(brca, names(mc3))
m <- as.matrix(mc3[,..keep]); rownames(m) <- mc3$gene; storage.mode(m)<-"numeric"
all <- data.table(gene=rownames(m), mut_freq=rowMeans(m>0,na.rm=TRUE))
all[, in_network := gene %in% nodes[type %in% c("Gene","TF"), name]]

## CDS length from refGene hg19 (max across transcripts)
rg <- fread(cmd=paste0("zcat ",CA,"/refGene_hg19.txt.gz"), header=FALSE, select=c(3,7,8,10,11,13))
setnames(rg, c("chrom","cdsStart","cdsEnd","exonStarts","exonEnds","sym"))
rg <- rg[chrom %in% paste0("chr",c(1:22,"X","Y")) & cdsEnd > cdsStart]
cdslen <- function(es,ee,cs,ce){
  a <- as.integer(strsplit(es,",")[[1]]); b <- as.integer(strsplit(ee,",")[[1]])
  a2 <- pmax(a,cs); b2 <- pmin(b,ce); sum(pmax(0,b2-a2)) }
rg[, len := mapply(cdslen, exonStarts, exonEnds, cdsStart, cdsEnd)]
gl <- rg[, .(cds_len=max(len)), by=sym]
cat("genes with CDS length:", nrow(gl), "\n")

all <- merge(all, gl, by.x="gene", by.y="sym")
cat("MC3 genes with length:", nrow(all), "| network:", sum(all$in_network), "\n")
net <- all[in_network==TRUE]; oth <- all[in_network==FALSE]
cat("median CDS length: network", median(net$cds_len), " other", median(oth$cds_len),
    " (Wilcoxon p =", format.pval(wilcox.test(net$cds_len,oth$cds_len)$p.value), ")\n")

## length-stratified permutation: sample background genes matched on CDS-length decile
set.seed(1)
brk <- unique(quantile(all$cds_len, probs=seq(0,1,0.05)))
all[, bin := cut(cds_len, brk, include.lowest=TRUE)]
tab <- net[, .N, by=.(bin=cut(cds_len, brk, include.lowest=TRUE))]
pool <- split(oth$mut_freq, cut(oth$cds_len, brk, include.lowest=TRUE))
obs <- mean(net$mut_freq)
nperm <- 10000
perm <- replicate(nperm, {
  sum(sapply(seq_len(nrow(tab)), function(i){
    p <- pool[[as.character(tab$bin[i])]]
    if(is.null(p)||!length(p)) return(0)
    sum(sample(p, tab$N[i], replace=TRUE)) })) / sum(tab$N) })
cat("\n=== LENGTH-MATCHED PERMUTATION (", nperm, "draws ) ===\n")
cat("observed mean mutation freq of network genes:", signif(obs,4), "\n")
cat("length-matched null mean:", signif(mean(perm),4),
    " sd:", signif(sd(perm),4), " 95% CI:", signif(quantile(perm,c(.025,.975)),4), "\n")
pv <- (sum(perm >= obs)+1)/(nperm+1)
cat("empirical one-sided p =", signif(pv,4), " | fold over length-matched null =",
    signif(obs/mean(perm),3), "\n")
cat("naive (unmatched) background mean:", signif(mean(oth$mut_freq),4),
    " => naive fold =", signif(obs/mean(oth$mut_freq),3), "\n")
fwrite(data.table(metric=c("obs_mean_freq","lengthmatched_null_mean","lengthmatched_null_sd",
   "lengthmatched_null_lo95","lengthmatched_null_hi95","perm_p_onesided","fold_vs_lengthmatched",
   "naive_background_mean","fold_vs_naive","n_perm","median_cds_len_network","median_cds_len_other"),
  value=c(obs,mean(perm),sd(perm),quantile(perm,.025),quantile(perm,.975),pv,obs/mean(perm),
   mean(oth$mut_freq),obs/mean(oth$mut_freq),nperm,median(net$cds_len),median(oth$cds_len))),
  file.path(OUT,"brca_mutation_lengthmatched_permutation.csv"))
cat("wrote brca_mutation_lengthmatched_permutation.csv\n")
