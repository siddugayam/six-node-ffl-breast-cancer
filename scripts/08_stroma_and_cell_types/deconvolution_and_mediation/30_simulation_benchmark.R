#!/usr/bin/env Rscript
## 30_simulation_benchmark.R
## Which fibroblast estimate is actually accurate? Benchmark every method against
## pseudo-bulk mixtures with KNOWN composition, built from the Wu et al. 2021 atlas.
## Patient-level hold-out: the Wu-derived signature matrix / InstaPrism reference are
## rebuilt from half the patients, and mixtures are simulated from the OTHER half, so
## the reference-based methods are not scored on their own training data.
## CAVEAT recorded in the log: simulated mixtures inherit no bulk-tissue technical
## noise, no ambient RNA and no cell-type mRNA-content differences (each cell type
## contributes CPM-normalised mRNA in proportion to its fraction). Absolute accuracy
## here is therefore an upper bound; the ranking between methods is the usable output.
suppressPackageStartupMessages({library(matrixStats); library(data.table); library(nnls)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/30_simulation.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")
set.seed(20260909)

pc <- as.matrix(read.delim(file.path(D, "wu2021_sum_umi_by_patient_celltype.tsv"), row.names = 1, check.names = FALSE))
bulkgenes <- rownames(readRDS(file.path(D, "bulk_log2.rds")))
pat <- sub("\\|.*$", "", colnames(pc)); ct <- sub("^.*\\|", "", colnames(pc))
cts <- sort(unique(ct)); pats <- sort(unique(pat))
logf("pseudobulk profiles: ", ncol(pc), " from ", length(pats), " patients, ", length(cts), " cell types")
trainP <- sort(sample(pats, floor(length(pats) / 2)))
testP  <- setdiff(pats, trainP)
logf("TRAIN patients (", length(trainP), "): ", paste(trainP, collapse = ", "))
logf("TEST  patients (", length(testP), "): ", paste(testP, collapse = ", "))

cpm <- sweep(pc, 2, colSums(pc), "/") * 1e6

## ------------------------------------------------------- simulate test mixtures ----
NSIM <- 300
P <- matrix(rgamma(NSIM * length(cts), 0.7, 1), NSIM, length(cts)); P <- P / rowSums(P)
colnames(P) <- cts
MIX <- matrix(0, nrow(cpm), NSIM, dimnames = list(rownames(cpm), paste0("sim", seq_len(NSIM))))
for (i in seq_len(NSIM)) for (k in seq_along(cts)) {
  cand <- which(pat %in% testP & ct == cts[k])
  if (!length(cand)) { P[i, k] <- 0; next }
  MIX[, i] <- MIX[, i] + P[i, k] * cpm[, sample(cand, 1)]
}
P <- P / rowSums(P)
MIX <- sweep(MIX, 2, colSums(MIX), "/") * 1e6
MIX <- MIX[intersect(rownames(MIX), bulkgenes), , drop = FALSE]
LMIX <- log2(MIX + 1)
logf("simulated ", NSIM, " mixtures on ", nrow(MIX), " genes; true CAF fraction median ",
     round(median(P[, "CAFs"]), 3), " range ", round(min(P[, "CAFs"]), 3), " - ", round(max(P[, "CAFs"]), 3))
write.table(data.frame(Gene = rownames(LMIX), LMIX, check.names = FALSE),
            file.path(D, "sim_mix_log2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## --------------------------------------------- train-half Wu signature + reference ---
mean_by <- function(cols) rowMeans(cpm[, cols, drop = FALSE])
trIdx <- which(pat %in% trainP)
refTr <- sapply(cts, function(k) mean_by(intersect(trIdx, which(ct == k))))
refTr <- refTr[intersect(rownames(refTr), bulkgenes), ]
fc <- sapply(cts, function(k) (refTr[, k] + 1) / (rowMaxs(refTr[, setdiff(cts, k), drop = FALSE]) + 1))
dimnames(fc) <- dimnames(refTr)
build <- function(G) { sel <- character(0)
  for (k in cts) { ok <- which(fc[, k] > 2 & refTr[, k] >= 10); ok <- ok[order(-fc[ok, k])]
                   sel <- union(sel, rownames(refTr)[head(ok, G)]) }
  refTr[sel, , drop = FALSE] }
Gs <- seq(25, 300, 25); kp <- sapply(Gs, function(G) kappa(build(G), exact = TRUE))
SIGtr <- build(Gs[which.min(kp)])
logf("train-half signature: G = ", Gs[which.min(kp)], " -> ", nrow(SIGtr), " genes, kappa ", round(min(kp), 2))

## ------------------------------------------------------------------ run methods ----
EST <- list()
EST$TRUE_CAF <- P[, "CAFs"]
try({ library(estimate)
  f1 <- file.path(D, "sim_filtered.gct"); f2 <- file.path(D, "sim_score.gct")
  filterCommonGenes(input.f = file.path(D, "sim_mix_log2.txt"), output.f = f1, id = "GeneSymbol")
  estimateScore(input.ds = f1, output.ds = f2, platform = "illumina")
  e <- t(as.matrix(read.delim(f2, skip = 2, row.names = 1, check.names = FALSE)[, -1]))
  EST$ESTIMATE_StromalScore <- e[, "StromalScore"] }, silent = FALSE)
try({ library(xCell); x <- t(xCellAnalysis(LMIX, rnaseq = TRUE, parallel.sz = 8))
  EST$xCell_Fibroblasts <- x[, "Fibroblasts"]; EST$xCell_StromaScore <- x[, "StromaScore"] }, silent = FALSE)
try({ library(MCPcounter); m <- t(MCPcounter.estimate(LMIX, featuresType = "HUGO_symbols"))
  EST$MCPcounter_Fibroblasts <- m[, "Fibroblasts"] }, silent = FALSE)
try({ library(EPIC); e <- EPIC(bulk = MIX, reference = "TRef", withOtherCells = TRUE)$cellFractions
  EST$EPIC_CAFs <- e[, "CAFs"] }, silent = FALSE)
try({ library(ConsensusTME); c1 <- t(consensusTMEAnalysis(as.matrix(LMIX), cancer = "BRCA", statMethod = "ssgsea"))
  EST$ConsensusTME_Fibroblasts <- c1[, "Fibroblasts"] }, silent = FALSE)
try({ source(file.path(BASE, "scripts/08_stroma_and_cell_types/deconvolution_and_mediation/23_cibersort_impl.R"))
  r <- run_cibersort(SIGtr, MIX, QN = FALSE, cores = 8)
  EST$CBSX_Wu_CAFs <- r$fractions[, "CAFs"] }, silent = FALSE)
try({ gi <- intersect(rownames(SIGtr), rownames(MIX))
  w <- t(sapply(seq_len(ncol(MIX)), function(j) { z <- nnls(SIGtr[gi, ], MIX[gi, j])$x; z / sum(z) }))
  colnames(w) <- cts; EST$NNLS_Wu_CAFs <- w[, "CAFs"] }, silent = FALSE)
try({ library(InstaPrism)
  csTr <- colnames(pc)[trIdx]
  rp <- refPrepare(sc_Expr = pc[intersect(rownames(pc), rownames(MIX)), csTr],
                   cell.type.labels = ct[trIdx], cell.state.labels = csTr)
  gi <- intersect(rownames(MIX), rownames(rp@phi.cs))
  ip <- InstaPrism(bulk_Expr = MIX[gi, ], refPhi_cs = rp, n.iter = 100, n.core = 8, verbose = FALSE)
  EST$InstaPrism_Wu_CAFs <- t(ip@Post.ini.ct@theta)[, "CAFs"] }, silent = FALSE)
try({ suppressPackageStartupMessages(library(GSVA))
  sets <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
  nodes <- read.delim(file.path(BASE, "data/canonical_nodes.tsv"), stringsAsFactors = FALSE)$name
  gmt <- readLines(system.file("extdata", "SI_geneset.gmt", package = "estimate"))
  strom <- strsplit(gmt[grep("StromalSignature", gmt)], "\t")[[1]][-c(1, 2)]
  SL <- list(CAF_A_estimate_noCOL_noNodes = intersect(setdiff(setdiff(strom, grep("^COL", strom, value = TRUE)), nodes), rownames(LMIX)),
             CAF_scRNA = intersect(setdiff(sets$CAF_scRNA_50, nodes), rownames(LMIX)),
             FARMER_clean = intersect(setdiff(setdiff(sets$FARMER_STROMAL, grep("^COL", sets$FARMER_STROMAL, value = TRUE)), nodes), rownames(LMIX)))
  zz <- sapply(SL, function(g) { X <- LMIX[g, , drop = FALSE]; colMeans((X - rowMeans(X)) / rowSds(X)) })
  EST$meanZ_CAF_A_estimate <- zz[, 1]; EST$meanZ_CAF_scRNA <- zz[, 2]; EST$meanZ_FARMER_clean <- zz[, 3]
  sg <- t(gsva(ssgseaParam(as.matrix(LMIX), SL, normalize = TRUE)))
  EST$ssGSEA_CAF_scRNA <- sg[, "CAF_scRNA"]
  EST$ssGSEA_ESTIMATE_stromal_clean <- sg[, "CAF_A_estimate_noCOL_noNodes"] }, silent = FALSE)

M <- do.call(cbind, EST)
saveRDS(list(P = P, EST = M), file.path(D, "sim_benchmark.rds"))
tr <- P[, "CAFs"]
B <- rbindlist(lapply(setdiff(colnames(M), "TRUE_CAF"), function(j) {
  v <- M[, j]
  data.table(estimate = j, spearman_vs_true_CAF = cor(v, tr, method = "spearman"),
             pearson_vs_true_CAF = cor(v, tr),
             spearman_vs_true_PVL = cor(v, P[, "PVL"], method = "spearman"),
             spearman_vs_true_Endo = cor(v, P[, "Endothelial"], method = "spearman"),
             spearman_vs_true_CancerEpi = cor(v, P[, "Cancer Epithelial"], method = "spearman"),
             spearman_vs_true_Myeloid = cor(v, P[, "Myeloid"], method = "spearman"),
             spearman_vs_true_stroma_CAFplusPVL = cor(v, P[, "CAFs"] + P[, "PVL"], method = "spearman"),
             rmse_if_fraction = if (all(v >= 0 & v <= 1)) sqrt(mean((v - tr)^2)) else NA_real_)
}))
setorder(B, -spearman_vs_true_CAF)
fwrite(B, file.path(BASE, "results/v3/deconv_simulation_benchmark.csv"))
fwrite(data.table(sim = rownames(P), as.data.table(P)), file.path(BASE, "results/v3/deconv_simulation_true_fractions.csv"))
logf("\n=== ACCURACY AGAINST KNOWN CAF FRACTION (", NSIM, " held-out-patient pseudobulk mixtures) ===")
logf(sprintf("%-32s %8s %8s %8s %8s %8s %8s", "estimate", "rho_CAF", "r_CAF", "rho_PVL", "rho_Endo", "rho_Epi", "rho_CAF+PVL"))
for (i in seq_len(nrow(B))) logf(sprintf("%-32s %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f",
  B$estimate[i], B$spearman_vs_true_CAF[i], B$pearson_vs_true_CAF[i], B$spearman_vs_true_PVL[i],
  B$spearman_vs_true_Endo[i], B$spearman_vs_true_CancerEpi[i], B$spearman_vs_true_stroma_CAFplusPVL[i]))
logf("DONE 30")
