## 21_metabolic_scores.R ---------------------------------------------------
## Per-sample ssGSEA + GSVA scores for the metabolic gene-set space,
## over all 1,211 TCGA-BRCA samples (1,097 primary tumours + 114 normals).
suppressPackageStartupMessages({library(data.table); library(GSVA); library(BiocParallel)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
E <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds"))
stopifnot(!any(duplicated(rownames(E))), is.matrix(E))
sdv <- apply(E, 1, sd); Ef <- E[sdv > 0, , drop=FALSE]
cat("expr:", dim(E), "-> non-constant:", dim(Ef), "\n")
sets <- readRDS(file.path(RES,"metabolic_genesets.rds"))
sets <- lapply(sets, function(g) intersect(g, rownames(Ef)))
sets <- sets[lengths(sets) >= 5]
cat("sets scored:", length(sets), "\n")
bp <- MulticoreParam(workers=12, RNGseed=1)

t0 <- Sys.time()
ss <- gsva(ssgseaParam(Ef, sets, alpha=0.25, normalize=TRUE, minSize=5, maxSize=Inf), BPPARAM=bp)
cat("ssGSEA done", format(Sys.time()-t0), " dim:", dim(ss), "\n")
t0 <- Sys.time()
gv <- gsva(gsvaParam(Ef, sets, kcdf="Gaussian", minSize=5, maxSize=Inf), BPPARAM=bp)
cat("GSVA done", format(Sys.time()-t0), " dim:", dim(gv), "\n")

stopifnot(identical(colnames(ss), colnames(E)), identical(colnames(gv), colnames(E)))
common <- intersect(rownames(ss), rownames(gv))
agree <- sapply(common, function(s) suppressWarnings(cor(ss[s,], gv[s,], method="spearman")))
cat("ssGSEA vs GSVA per-set Spearman: median", round(median(agree,na.rm=TRUE),3),
    " range", paste(round(range(agree,na.rm=TRUE),3), collapse=" to "), "\n")

saveRDS(ss, file.path(RES,"metabolic_ssgsea_scores.rds"))
saveRDS(gv, file.path(RES,"metabolic_gsva_scores.rds"))
fwrite(data.table(set_id=rownames(ss), as.data.table(ss)), file.path(RES,"metabolic_ssgsea_scores.csv"))
fwrite(data.table(set_id=rownames(gv), as.data.table(gv)), file.path(RES,"metabolic_gsva_scores.csv"))
fwrite(data.table(set_id=common, spearman_ssgsea_vs_gsva=agree[common],
                  n_genes=lengths(sets)[common]),
       file.path(RES,"metabolic_score_method_agreement.csv"))
cat("DONE 21\n")
