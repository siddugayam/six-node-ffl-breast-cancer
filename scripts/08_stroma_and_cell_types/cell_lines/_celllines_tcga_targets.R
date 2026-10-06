# Written and run by celllines_09_network_sign_tcga.py, which holds this code; running that script regenerates
# this file.

g <- readRDS("/path/to/revision/data/brca_gene_expr.rds"); m <- readRDS("/path/to/revision/data/brca_mirna_expr_canonical.rds")
p <- readRDS("/path/to/revision/data/brca_pheno.rds")
tum <- p$sample[p$sample_type == "Primary Tumor"]
both <- intersect(intersect(colnames(g), tum), intersect(colnames(m), tum))
cat("N_paired", length(both), "\n")
gn <- readLines("/path/to/revision/results/v2/_need_genes.txt"); gn <- gn[gn %in% rownames(g)]
cat("N_genes_found", length(gn), "\n")
write.table(g[gn, both], "/path/to/revision/results/v2/_tcga_targets.tsv", sep="\t", quote=FALSE)
write.table(m[, both], "/path/to/revision/results/v2/_tcga_allmirs.tsv", sep="\t", quote=FALSE)
