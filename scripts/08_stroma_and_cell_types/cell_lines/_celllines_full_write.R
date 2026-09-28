
x <- readRDS("/path/to/revision/results/v2/_tcga_full_paired.rds")
write.table(x$g, gzfile("/path/to/revision/results/v2/_tcga_g.tsv.gz"), sep="\t", quote=FALSE)
write.table(x$m, "/path/to/revision/results/v2/_tcga_m.tsv", sep="\t", quote=FALSE)
