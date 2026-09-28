
g <- readRDS("/path/to/revision/data/brca_gene_expr.rds"); m <- readRDS("/path/to/revision/data/brca_mirna_expr_canonical.rds")
p <- readRDS("/path/to/revision/data/brca_pheno.rds")
tum <- p$sample[p$sample_type == "Primary Tumor"]
both <- intersect(intersect(colnames(g), tum), intersect(colnames(m), tum))
cat("N_paired", length(both), "\n")
saveRDS(list(g = g[, both], m = m[, both]), "/path/to/revision/results/v2/_tcga_full_paired.rds")
