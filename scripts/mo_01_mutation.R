#!/usr/bin/env Rscript
## Part B: somatic non-silent mutation frequency in TCGA-BRCA for network genes/TFs
suppressPackageStartupMessages({library(data.table)})
RV  <- "/path/to/revision"
OUT <- file.path(RV,"results/multiomics")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

nodes <- fread(file.path(RV,"data/canonical_nodes.tsv"))
cat("nodes:", nrow(nodes), "| by type:\n"); print(table(nodes$type))
prot <- nodes[type %in% c("Gene","TF")]
cat("protein-coding network nodes (Gene+TF):", nrow(prot), "\n")

## phenotype -> BRCA primary tumour samples
ph <- fread(cmd=paste0("zcat /path/to/home/Desktop/DD/R_GPR/ML/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"))
setnames(ph, make.names(names(ph)))
cat("pheno cols:", paste(names(ph),collapse=","), "\n")
brca <- ph[grepl("breast", X_primary_disease, ignore.case=TRUE) & sample_type=="Primary Tumor", sample]
cat("BRCA primary tumour samples in pheno:", length(brca), "\n")

mc3 <- fread(cmd="zcat /path/to/home/Desktop/DD/R_GPR/ML/external_cohorts/MC3/mc3_nonsilentGene.xena.gz")
setnames(mc3, 1, "gene")
cat("MC3 matrix:", nrow(mc3), "genes x", ncol(mc3)-1, "samples\n")

keep <- intersect(brca, names(mc3))
cat("BRCA samples present in MC3:", length(keep), "\n")
stopifnot(length(keep) > 500)

m <- as.matrix(mc3[, ..keep])
rownames(m) <- mc3$gene
storage.mode(m) <- "numeric"
nmiss <- sum(is.na(m)); cat("NA cells:", nmiss, "\n")
nsamp <- ncol(m)
nmut  <- rowSums(m > 0, na.rm=TRUE)
freq  <- nmut / nsamp

res <- data.table(gene=rownames(m), n_mutated=nmut, n_samples=nsamp, mut_freq=freq)
res[, in_network := gene %in% prot$name]
res[, node_type  := nodes$type[match(gene, nodes$name)]]
res[is.na(node_type), node_type := "not_in_network"]

## network vs non-network comparison
netf <- res[in_network==TRUE, mut_freq]; othf <- res[in_network==FALSE, mut_freq]
cat("\n=== network vs background ===\n")
cat("network genes measured:", length(netf), " median freq:", median(netf),
    " mean:", mean(netf), "\n")
cat("other genes:", length(othf), " median freq:", median(othf), " mean:", mean(othf), "\n")
wt <- wilcox.test(netf, othf)
cat("Wilcoxon rank-sum p =", format.pval(wt$p.value), " W =", wt$statistic, "\n")
## proportion recurrently mutated (>=2%)
tb <- matrix(c(sum(netf>=0.02), sum(netf<0.02), sum(othf>=0.02), sum(othf<0.02)), nrow=2)
ft <- fisher.test(tb)
cat("genes with freq>=2%: network", sum(netf>=0.02),"/",length(netf),
    sprintf(" (%.1f%%)", 100*mean(netf>=0.02)),
    " | other", sum(othf>=0.02),"/",length(othf),
    sprintf(" (%.1f%%)", 100*mean(othf>=0.02)), "\n")
cat("Fisher OR =", ft$estimate, " p =", format.pval(ft$p.value), "\n")

## which network genes are missing from MC3
missing <- setdiff(prot$name, res$gene)
cat("\nnetwork protein-coding nodes absent from MC3 matrix:", length(missing), "\n")
if(length(missing)) cat(paste(missing, collapse=", "), "\n")

## save: network genes only, plus flag
outdt <- res[order(-mut_freq)]
fwrite(outdt[in_network==TRUE, .(gene, node_type, n_mutated, n_samples, mut_freq,
                                 background_median_freq=median(othf),
                                 background_mean_freq=mean(othf))],
       file.path(OUT,"brca_mutation_frequency.csv"))
cat("\nwrote brca_mutation_frequency.csv rows:", sum(res$in_network), "\n")

cat("\n=== TOP 20 most mutated network genes ===\n")
print(outdt[in_network==TRUE][1:20, .(gene,node_type,n_mutated,n_samples,mut_freq=round(mut_freq,4))])

## stats summary file
stats <- data.table(
  metric=c("n_brca_samples","n_network_genes_in_mc3","n_other_genes",
           "network_median_freq","other_median_freq","network_mean_freq","other_mean_freq",
           "wilcoxon_p","frac_network_ge2pct","frac_other_ge2pct","fisher_OR","fisher_p"),
  value=c(nsamp,length(netf),length(othf),median(netf),median(othf),mean(netf),mean(othf),
          wt$p.value, mean(netf>=0.02), mean(othf>=0.02), unname(ft$estimate), ft$p.value))
fwrite(stats, file.path(OUT,"brca_mutation_network_vs_background.csv"))
print(stats)
