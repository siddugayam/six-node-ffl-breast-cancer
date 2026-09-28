## Independent re-computation (this run) of edge-level Spearman rho in TCGA-BRCA
## for every miRNA_target and TF_target edge of the canonical network.
suppressWarnings(suppressMessages({library(data.table)}))
BASE <- "/path/to/revision"
g <- readRDS(file.path(BASE,"data/brca_gene_expr.rds"))
m <- readRDS(file.path(BASE,"data/brca_mirna_expr_canonical.rds"))
ph <- readRDS(file.path(BASE,"data/brca_pheno.rds"))
stopifnot(sum(duplicated(rownames(g)))==0, sum(duplicated(rownames(m)))==0)

tum <- ph$sample[ph$sample_type=="Primary Tumor"]
sg  <- intersect(colnames(g), tum); sm <- intersect(colnames(m), tum)
com <- intersect(sg, sm)
cat("primary tumours with gene assay:", length(sg), "\n")
cat("primary tumours with miRNA assay:", length(sm), "\n")
cat("paired (both assays):", length(com), "\n")

G <- g[, com, drop=FALSE]; M <- m[, com, drop=FALSE]
## drop zero-variance rows
G <- G[apply(G,1,function(x) sd(x,na.rm=TRUE)>0), , drop=FALSE]
M <- M[apply(M,1,function(x) sd(x,na.rm=TRUE)>0), , drop=FALSE]
cat("genes with variance:", nrow(G), " miRNAs with variance:", nrow(M), "\n")
Gr <- t(apply(G,1,rank)); Mr <- t(apply(M,1,rank))

E <- fread(file.path(BASE,"data/canonical_edges.tsv"))
cat("edges:", nrow(E), "\n"); print(table(E$edge_type))

sub <- E[edge_type %in% c("miRNA_target","TF_target","TF_miRNA","gene_gene","miRNA_miRNA")]
res <- vector("list", nrow(sub))
n <- length(com)
for (i in seq_len(nrow(sub))) {
  s <- sub$source[i]; t <- sub$target[i]; et <- sub$edge_type[i]
  xs <- if (s %in% rownames(Mr)) Mr[s,] else if (s %in% rownames(Gr)) Gr[s,] else NULL
  ys <- if (t %in% rownames(Mr)) Mr[t,] else if (t %in% rownames(Gr)) Gr[t,] else NULL
  if (is.null(xs) || is.null(ys)) { res[[i]] <- list(rho=NA_real_, p=NA_real_); next }
  ct <- suppressWarnings(cor.test(xs, ys, method="pearson"))
  res[[i]] <- list(rho=unname(ct$estimate), p=ct$p.value)
}
sub[, rho := sapply(res, function(z) z$rho)]
sub[, p   := sapply(res, function(z) z$p)]
sub[, n_samples := n]
cat("measurable edges:", sum(!is.na(sub$rho)), "of", nrow(sub), "\n")
cat("measurable by type:\n"); print(sub[!is.na(rho), .N, by=edge_type])
## concordance definition: miRNA_target expected negative
mt <- sub[edge_type=="miRNA_target" & !is.na(rho)]
cat("miRNA_target measurable:", nrow(mt), " frac rho<0:", round(mean(mt$rho<0),4), "\n")
tt <- sub[edge_type=="TF_target" & !is.na(rho)]
cat("TF_target measurable:", nrow(tt), " frac rho>0:", round(mean(tt$rho>0),4), "\n")
fwrite(sub, file.path(BASE,"results/v2/v2d_edge_rho_recomputed.csv"))
cat("WROTE results/v2/v2d_edge_rho_recomputed.csv\n")
