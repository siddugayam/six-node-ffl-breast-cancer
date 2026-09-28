#!/usr/bin/env Rscript
## 27_assemble_scores.R -- assemble every per-sample cell-type estimate into one
## table, add the signature-score methods the project already used (so the new
## methods are compared with the incumbent on identical samples), and write:
##   results/v3/deconv_all_celltype_estimates.csv   long table, every method x feature
##   results/v3/deconv_fibroblast_estimates.csv     wide table, one column per fibroblast/stromal estimate
suppressPackageStartupMessages({library(matrixStats); library(GSVA); library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/27_assemble.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")

LOG2 <- readRDS(file.path(D, "bulk_log2.rds")); SS <- colnames(LOG2)
nodes <- read.delim(file.path(BASE, "data/canonical_nodes.tsv"), stringsAsFactors = FALSE)$name
sets  <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))

fs <- list.files(file.path(D, "out"), pattern = "\\.rds$", full.names = TRUE)
fs <- fs[!grepl("_stats\\.rds$", fs)]
mats <- list()
for (f in fs) {
  nm <- sub("\\.rds$", "", basename(f)); M <- readRDS(f)
  stopifnot(identical(rownames(M), SS))
  mats[[nm]] <- M
  logf(sprintf("%-26s %4d features", nm, ncol(M)))
}

## ---- signature-score methods used elsewhere in the project, recomputed here ----
zsc <- function(gs) { X <- LOG2[intersect(gs, rownames(LOG2)), , drop = FALSE]
                      colMeans((X - rowMeans(X)) / rowSds(X)) }
gmt   <- readLines(system.file("extdata", "SI_geneset.gmt", package = "estimate"))
strom <- strsplit(gmt[grep("StromalSignature", gmt)], "\t")[[1]][-c(1, 2)]
sigA  <- intersect(setdiff(setdiff(strom, grep("^COL", strom, value = TRUE)), nodes), rownames(LOG2))
farmer_clean <- intersect(setdiff(setdiff(sets$FARMER_STROMAL, grep("^COL", sets$FARMER_STROMAL, value = TRUE)), nodes), rownames(LOG2))
cafsc <- intersect(setdiff(sets$CAF_scRNA_50, nodes), rownames(LOG2))
cafmk <- intersect(setdiff(c("DCN","LUM","FAP","PDGFRB","THY1","POSTN"), nodes), rownames(LOG2))
logf("signature sizes: CAF_A_estimate_noCOL_noNodes=", length(sigA), " FARMER_clean=", length(farmer_clean),
     " CAF_scRNA=", length(cafsc), " CAF_B_markers=", length(cafmk))
SIGSETS <- list(CAF_A_estimate_noCOL_noNodes = sigA, FARMER_clean = farmer_clean,
                CAF_scRNA = cafsc, CAF_B_markers = cafmk)
mz <- sapply(SIGSETS, zsc)
mats[["meanZ"]] <- mz

gp <- gsvaParam(as.matrix(LOG2), SIGSETS, kcdf = "Gaussian")
ss <- t(GSVA::gsva(ssgseaParam(as.matrix(LOG2), SIGSETS, normalize = TRUE)))
colnames(ss) <- paste0("ssGSEA_", colnames(ss))
mats[["ssGSEA"]] <- ss[SS, , drop = FALSE]
logf("ssGSEA scores computed for ", ncol(ss), " signatures")

## ---------------------------------------------------------------- long table ----
long <- rbindlist(lapply(names(mats), function(nm) {
  M <- mats[[nm]]
  data.table(method = nm, feature = rep(colnames(M), each = nrow(M)),
             sample = rep(rownames(M), ncol(M)), value = as.vector(M))
}))
fwrite(long, file.path(BASE, "results/v3/deconv_all_celltype_estimates.csv"))
logf("WROTE deconv_all_celltype_estimates.csv rows=", nrow(long),
     " methods=", length(unique(long$method)), " features=", uniqueN(paste(long$method, long$feature)))

## ------------------------------------------------ fibroblast / stromal estimates ----
FIB <- list(
  ESTIMATE_StromalScore        = mats$ESTIMATE[, "StromalScore"],
  ESTIMATE_ESTIMATEScore       = mats$ESTIMATE[, "ESTIMATEScore"],
  xCell_Fibroblasts            = mats$xCell[, "Fibroblasts"],
  xCell_StromaScore            = mats$xCell[, "StromaScore"],
  xCell_MicroenvironmentScore  = mats$xCell[, "MicroenvironmentScore"],
  MCPcounter_Fibroblasts       = mats$MCPcounter[, "Fibroblasts"],
  EPIC_CAFs                    = mats$EPIC[, "CAFs"],
  ConsensusTME_Fibroblasts     = mats$ConsensusTME[, "Fibroblasts"],
  CBSX_Wu_CAFs                 = mats$CBSX_Wu[, "CAFs"],
  CBSX_Wu_nocol_CAFs           = mats$CBSX_Wu_nocollagen[, "CAFs"],
  NNLS_Wu_CAFs                 = mats$NNLS_Wu[, "CAFs"],
  ssGSEA_CAF_scRNA             = mats$ssGSEA[, "ssGSEA_CAF_scRNA"],
  ssGSEA_ESTIMATE_stromal_clean= mats$ssGSEA[, "ssGSEA_CAF_A_estimate_noCOL_noNodes"],
  meanZ_CAF_A_estimate         = mats$meanZ[, "CAF_A_estimate_noCOL_noNodes"],
  meanZ_FARMER_clean           = mats$meanZ[, "FARMER_clean"],
  meanZ_CAF_scRNA              = mats$meanZ[, "CAF_scRNA"],
  meanZ_CAF_B_markers          = mats$meanZ[, "CAF_B_markers"])
for (nm in c("InstaPrism_Wu", "InstaPrism_Wu_nocollagen"))
  if (!is.null(mats[[nm]])) FIB[[paste0(nm, "_CAFs")]] <- mats[[nm]][, "CAFs"]
FIBM <- do.call(cbind, FIB); rownames(FIBM) <- SS
saveRDS(FIBM, file.path(D, "fibroblast_estimates.rds"))
fwrite(data.table(sample = SS, as.data.table(FIBM)), file.path(BASE, "results/v3/deconv_fibroblast_estimates.csv"))
logf("WROTE deconv_fibroblast_estimates.csv : ", ncol(FIBM), " fibroblast/stromal estimates x ", nrow(FIBM), " tumours")
logf("columns: ", paste(colnames(FIBM), collapse = ", "))
for (j in colnames(FIBM)) logf(sprintf("  %-32s median %10.4f  IQR %10.4f - %10.4f  n_zero %d",
    j, median(FIBM[, j]), quantile(FIBM[, j], .25), quantile(FIBM[, j], .75), sum(FIBM[, j] == 0)))
logf("DONE 27")
