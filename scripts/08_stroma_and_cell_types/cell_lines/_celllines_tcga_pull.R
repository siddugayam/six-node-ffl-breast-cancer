# Written and run by celllines_02_corr.py, which holds this code; running that script regenerates this file.

g <- readRDS("/path/to/revision/data/brca_gene_expr.rds")
m <- readRDS("/path/to/revision/data/brca_mirna_expr_canonical.rds")
p <- readRDS("/path/to/revision/data/brca_pheno.rds")
tum <- p$sample[p$sample_type == "Primary Tumor"]
gt <- intersect(colnames(g), tum); mt <- intersect(colnames(m), tum)
both <- intersect(gt, mt)
cat("N_gene_tumours", length(gt), "\n"); cat("N_paired", length(both), "\n")
genes <- c("COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1","MYB","MKL1","STAT6","TFAP2A",
           "EZH2","DCN","LUM","FAP","THY1")
genes <- genes[genes %in% rownames(g)]
write.table(g[genes, gt], "/path/to/revision/results/v2/_tcga_genes_tumour.tsv", sep="\t", quote=FALSE)
mirs <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-let-7b","hsa-let-7e",
          "hsa-miR-143","hsa-miR-218","hsa-miR-133a","hsa-miR-133b","hsa-miR-101",
          "hsa-miR-21")
mirs <- mirs[mirs %in% rownames(m)]
write.table(m[mirs, both], "/path/to/revision/results/v2/_tcga_mirs_tumour.tsv", sep="\t", quote=FALSE)
write.table(g[genes, both], "/path/to/revision/results/v2/_tcga_genes_paired.tsv", sep="\t", quote=FALSE)
# within-sample percentile rank of the collagens across ALL genes
r <- apply(g[, gt], 2, function(v) rank(v)/length(v))
write.table(t(r[c("COL1A1","COL3A1"), ]), "/path/to/revision/results/v2/_tcga_collagen_pct.tsv", sep="\t", quote=FALSE)
cat("N_genes_matrix", nrow(g), "\n")
