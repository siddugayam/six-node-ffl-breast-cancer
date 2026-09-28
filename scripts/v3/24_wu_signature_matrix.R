#!/usr/bin/env Rscript
## 24_wu_signature_matrix.R
## Build a BREAST-TUMOUR-SPECIFIC signature matrix for 9 major cell types from the
## Wu et al. 2021 (GSE176078) atlas, following the CIBERSORT signature-matrix
## construction recipe: (i) per-cell-type differential expression across patients,
## (ii) rank candidate markers by fold change, (iii) sweep the number of markers per
## cell type and keep the matrix with the smallest 2-norm condition number.
## Two matrices are produced: the full one and a COLLAGEN-FREE one (all COL* genes
## removed) so the mediation analysis can be run with a fibroblast estimate that
## cannot contain the outcome variable.
suppressPackageStartupMessages({library(matrixStats); library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/24_wu_signature_matrix.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")

meancpm <- as.matrix(read.delim(file.path(D, "wu2021_meancpm_by_celltype.tsv"), row.names = 1, check.names = FALSE))
pc      <- as.matrix(read.delim(file.path(D, "wu2021_sum_umi_by_patient_celltype.tsv"), row.names = 1, check.names = FALSE))
logf("mean-CPM reference: ", nrow(meancpm), " genes x ", ncol(meancpm), " cell types")
logf("patient x celltype pseudobulks (>=20 cells): ", ncol(pc))

bulk <- readRDS(file.path(BASE, "cache/v3/deconv/bulk_log2.rds"))
common <- intersect(rownames(meancpm), rownames(bulk))
logf("genes shared with the TCGA bulk matrix: ", length(common))
meancpm <- meancpm[common, ]; pc <- pc[common, , drop = FALSE]

## CPM-normalise each patient x celltype pseudobulk, log for testing
pcc <- sweep(pc, 2, colSums(pc), "/") * 1e6
lpc <- log2(pcc + 1)
grp <- sub("^.*\\|", "", colnames(pc))
cts <- colnames(meancpm)
logf("pseudobulk groups per cell type: ", paste(sprintf("%s=%d", cts, table(grp)[cts]), collapse = ", "))

## per cell type: two-sided Wilcoxon of that type's patient pseudobulks vs all others
pv <- matrix(NA_real_, nrow(lpc), length(cts), dimnames = list(rownames(lpc), cts))
for (ct in cts) {
  i1 <- grp == ct; i2 <- !i1
  pv[, ct] <- sapply(seq_len(nrow(lpc)), function(g)
    tryCatch(wilcox.test(lpc[g, i1], lpc[g, i2])$p.value, error = function(e) NA_real_))
}
qv <- apply(pv, 2, p.adjust, method = "BH")
logf("genes with BH q<0.05 per cell type: ", paste(sprintf("%s=%d", cts, colSums(qv < 0.05, na.rm = TRUE)), collapse = ", "))

## fold change of the mean-CPM profile: own type vs the largest other type
fc <- sapply(cts, function(ct) {
  own <- meancpm[, ct]; oth <- rowMaxs(meancpm[, setdiff(cts, ct), drop = FALSE])
  (own + 1) / (oth + 1)
})
dimnames(fc) <- list(rownames(meancpm), cts)

build_sig <- function(G, exclude = character(0)) {
  sel <- character(0)
  for (ct in cts) {
    ok <- which(qv[, ct] < 0.05 & fc[, ct] > 2 & meancpm[, ct] >= 10 &
                !(rownames(meancpm) %in% exclude))
    ok <- ok[order(-fc[ok, ct])]
    sel <- union(sel, rownames(meancpm)[head(ok, G)])
  }
  meancpm[sel, , drop = FALSE]
}
Gs <- seq(25, 400, by = 25)
kap <- sapply(Gs, function(G) { S <- build_sig(G); kappa(S, exact = TRUE) })
ngn <- sapply(Gs, function(G) nrow(build_sig(G)))
for (i in seq_along(Gs)) logf(sprintf("G=%3d  genes=%4d  condition number = %8.2f", Gs[i], ngn[i], kap[i]))
Gbest <- Gs[which.min(kap)]
SIG <- build_sig(Gbest)
logf("CHOSEN G = ", Gbest, "  -> signature matrix ", nrow(SIG), " genes x ", ncol(SIG),
     " cell types, condition number ", round(min(kap), 2))

colg <- grep("^COL[0-9]", rownames(meancpm), value = TRUE)
kap2 <- sapply(Gs, function(G) { S <- build_sig(G, exclude = colg); kappa(S, exact = TRUE) })
Gbest2 <- Gs[which.min(kap2)]
SIGC <- build_sig(Gbest2, exclude = colg)
logf("COLLAGEN-FREE: excluded ", length(colg), " COL* genes; chosen G = ", Gbest2,
     " -> ", nrow(SIGC), " genes, condition number ", round(min(kap2), 2))
logf("collagens in full signature: ", paste(intersect(rownames(SIG), colg), collapse = ", "))
logf("COL1A1 in full signature: ", "COL1A1" %in% rownames(SIG),
     " ; COL3A1 in full signature: ", "COL3A1" %in% rownames(SIG))

write.table(data.frame(Gene = rownames(SIG), SIG, check.names = FALSE),
            file.path(D, "wu2021_signature_matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = rownames(SIGC), SIGC, check.names = FALSE),
            file.path(D, "wu2021_signature_matrix_nocollagen.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(data.frame(G = Gs, n_genes = ngn, condition_number = kap,
                     n_genes_nocol = sapply(Gs, function(G) nrow(build_sig(G, exclude = colg))),
                     condition_number_nocol = kap2),
          file.path(BASE, "results/v3/deconv_wu_signature_condition_numbers.csv"), row.names = FALSE)
mk <- do.call(rbind, lapply(cts, function(ct) {
  gs <- rownames(SIG)[apply(SIG, 1, which.max) == match(ct, cts)]
  data.frame(celltype = ct, n_markers_maxin_type = length(gs),
             top_markers = paste(head(gs[order(-fc[gs, ct])], 15), collapse = ", "))
}))
write.csv(mk, file.path(BASE, "results/v3/deconv_wu_signature_markers.csv"), row.names = FALSE)
for (i in seq_len(nrow(mk))) logf(sprintf("%-18s n=%3d  %s", mk$celltype[i], mk$n_markers_maxin_type[i], mk$top_markers[i]))
logf("DONE 24")
