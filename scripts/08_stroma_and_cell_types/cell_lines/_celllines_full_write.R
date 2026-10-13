# Written and run by celllines_11_purity_axis.py, which holds this code; running that script regenerates this
# file.

x <- readRDS("/path/to/six-node-ffl-breast-cancer/results/expression_validation/_tcga_full_paired.rds")
write.table(x$g, gzfile("/path/to/six-node-ffl-breast-cancer/results/expression_validation/_tcga_g.tsv.gz"), sep="\t", quote=FALSE)
write.table(x$m, "/path/to/six-node-ffl-breast-cancer/results/expression_validation/_tcga_m.tsv", sep="\t", quote=FALSE)
