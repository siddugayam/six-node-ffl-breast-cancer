#!/usr/bin/env Rscript
# 14_external_met500.R
# EXTERNAL COHORT 1: MET500 metastatic breast cancer (local files).
# normal DE in primary TCGA tumours cannot support a statement about metastasis.
# MET500 provides real metastatic breast biopsies, so hub expression and module
# scores can be compared: metastatic breast (MET500) vs primary breast (TCGA).
# Output: results/external_MET500_breast.csv

suppressPackageStartupMessages({
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
})
BASE <- "/path/to/revision"
MET  <- "/path/to/home/Desktop/DD/R_GPR/ML/external_cohorts/MET500"
LOG  <- file.path(BASE, "logs", "14_external_met500.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")

## ------------------------------------------------------------- metadata -----
meta <- fread(file.path(MET, "M.meta.plus.txt"))
logf("MET500 metadata rows=", nrow(meta), " cols=", paste(names(meta), collapse = ","))
br <- meta[tolower(cohort) == "brca"]
logf("MET500 breast (cohort==BRCA) samples=", nrow(br),
     " distinct patients (sample_source)=", uniqueN(br$sample_source))
logf("library prep: ", paste(sprintf("%s=%d",
      names(table(sub(".*?-(capt|poly)-.*", "\\1", br$Sample_id))),
      table(sub(".*?-(capt|poly)-.*", "\\1", br$Sample_id))), collapse = " "))
logf("biopsy sites: ", paste(sprintf("%s=%d", names(table(br$biopsy_tissue)),
                                     table(br$biopsy_tissue)), collapse = " "))

## ---------------------------------------------------- expression (ENSG) -----
logf("reading M.mx.log2.txt.gz (only the breast columns) ...")
hdr <- names(fread(cmd = sprintf("zcat %s | head -1", file.path(MET, "M.mx.log2.txt.gz"))))
keep <- c(hdr[1], intersect(hdr, br$Sample_id))
logf("columns requested: ", length(keep) - 1, " breast samples of ", length(hdr) - 1, " total")
E <- fread(cmd = sprintf("zcat %s", file.path(MET, "M.mx.log2.txt.gz")), select = keep)
setnames(E, 1, "ensg")
logf("MET500 breast expression matrix: ", nrow(E), " genes x ", ncol(E) - 1, " samples")

E[, ensg_base := sub("\\..*$", "", ensg)]
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(E$ensg_base),
                                              keytype = "ENSEMBL", columns = "SYMBOL"))
map <- as.data.table(map)[!is.na(SYMBOL)]
map <- map[, .SD[1], by = ENSEMBL]                       # first symbol per ENSG
E <- merge(E, map, by.x = "ensg_base", by.y = "ENSEMBL")
logf("ENSG -> SYMBOL mapped rows: ", nrow(E), " ; distinct symbols=", uniqueN(E$SYMBOL))

scols <- setdiff(names(E), c("ensg", "ensg_base", "SYMBOL"))
M <- as.matrix(E[, ..scols]); rownames(M) <- E$SYMBOL
# collapse duplicate symbols by max mean expression
ord <- order(rowMeans(M, na.rm = TRUE), decreasing = TRUE)
M <- M[ord, , drop = FALSE]; M <- M[!duplicated(rownames(M)), , drop = FALSE]
logf("MET500 symbol-level matrix: ", paste(dim(M), collapse = " x "),
     " ; value range ", paste(round(range(M, na.rm = TRUE), 2), collapse = " to "))

# one library per patient: prefer poly-A, else capture (avoids pseudo-replication)
libtype <- ifelse(grepl("-poly-", colnames(M)), "poly", "capt")
pat <- br$sample_source[match(colnames(M), br$Sample_id)]
sel <- data.table(col = colnames(M), pat = pat, lib = libtype)
sel[, rank := ifelse(lib == "poly", 1L, 2L)]
setorder(sel, pat, rank)
sel1 <- sel[, .SD[1], by = pat]
M1 <- M[, sel1$col, drop = FALSE]
logf("one library per patient: ", ncol(M1), " metastatic breast tumours (",
     sum(sel1$lib == "poly"), " polyA, ", sum(sel1$lib == "capt"), " capture)")
saveRDS(list(M = M1, meta = sel1, br = br), file.path(BASE, "data/met500_breast.rds"))

## ------------------------------------------ TCGA primary breast for compare --
gex <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
pheno <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
tum <- intersect(colnames(gex), pheno$sample[pheno$sample_type == "Primary Tumor"])
G <- gex[, tum, drop = FALSE]
logf("TCGA primary breast tumours for comparison: ", ncol(G))

## MET500 is log2(FPKM)-scale RNA-seq, TCGA is log2(RSEM norm_count+1): the two are
## NOT on a common absolute scale. Cross-cohort comparison therefore uses the
## WITHIN-SAMPLE PERCENTILE RANK of each gene among the genes shared by both
## platforms (0 = lowest expressed in that sample, 1 = highest). This is invariant
## to any monotone per-sample transformation, so it is comparable across platforms.
## (Per-gene z-scoring within cohort would force every mean difference to exactly 0.)
common_genes <- intersect(rownames(M1), rownames(G))
logf("genes shared MET500 & TCGA: ", length(common_genes))

pct_rank <- function(m) apply(m, 2, function(col) (rank(col, ties.method = "average") - 0.5) / length(col))
Rm <- pct_rank(M1[common_genes, , drop = FALSE]); rownames(Rm) <- common_genes
Rg <- pct_rank(G[common_genes, , drop = FALSE]);  rownames(Rg) <- common_genes
logf("within-sample percentile-rank matrices: MET500 ", paste(dim(Rm), collapse = "x"),
     " ; TCGA ", paste(dim(Rg), collapse = "x"))
saveRDS(list(Rm = Rm, Rg = Rg, common_genes = common_genes),
        file.path(BASE, "data/met500_vs_tcga_ranks.rds"))

nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
prot <- nodes[type != "miRNA", name]
tested <- intersect(prot, common_genes)
logf("network protein-coding nodes testable in MET500: ", length(tested), " / ", length(prot))

res <- rbindlist(lapply(tested, function(gn) {
  za <- as.numeric(Rm[gn, ]); zb <- as.numeric(Rg[gn, ])
  w <- suppressWarnings(stats::wilcox.test(za, zb))
  data.table(gene = gn, n_met = length(za), n_prim = length(zb),
             mean_pctrank_MET500_metastasis = mean(za, na.rm = TRUE),
             mean_pctrank_TCGA_primary = mean(zb, na.rm = TRUE),
             delta_pctrank_met_minus_prim = mean(za, na.rm = TRUE) - mean(zb, na.rm = TRUE),
             wilcox_p = w$p.value)
}))
res[, wilcox_q_BH := p.adjust(wilcox_p, method = "BH")]
setorder(res, wilcox_p)
fwrite(res, file.path(BASE, "results/external_MET500_breast.csv"))
logf("WROTE results/external_MET500_breast.csv rows=", nrow(res))
logf("MET500 metastasis vs TCGA primary: q<0.05 in ", sum(res$wilcox_q_BH < 0.05),
     " / ", nrow(res), " network protein-coding nodes")

named <- c("NFKB1", "RELA", "SP1", "ETS1", "COL1A1", "COL3A1", "VEGFA", "CCND2",
           "MYC", "E2F1", "TP53")
logf("---- manuscript's named protein-coding hubs, metastasis vs primary ----")
for (gn in named) {
  x <- res[gene == gn]
  if (nrow(x) == 0) { logf(sprintf("%-8s NOT MEASURED in MET500", gn)); next }
  logf(sprintf("%-8s delta pct-rank (met - prim) = %+.4f  Wilcoxon p=%.3g q=%.3g",
               gn, x$delta_pctrank_met_minus_prim, x$wilcox_p, x$wilcox_q_BH))
}
logf("DONE")
