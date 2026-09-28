#!/usr/bin/env Rscript
## 45_deconv2_assemble.R
## (1) run the collagen-free variants of the new least-squares / dtangle methods, so
##     that every new fibroblast estimate also exists in a form whose signature cannot
##     contain COL1A1 or COL3A1;
## (2) assemble EVERY fibroblast/stromal estimate available for the 1,097 TCGA-BRCA
##     primary tumours -- the 26 from the first pass, the new algorithmic families, and
##     the two NON-TRANSCRIPTOMIC estimates (ABSOLUTE DNA purity, pathologist slide
##     scoring) -- into one matrix;
## (3) write the full pairwise correlation matrices (Spearman and Pearson) -- part (B).
suppressPackageStartupMessages({library(data.table); library(limSolve); library(quadprog); library(MASS)})
.libPaths(c(path.expand("~/Rlib_deconv2"), .libPaths()))
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/45_deconv2_assemble.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
OUT  <- file.path(BASE, "cache/v3/deconv2/out")

TPM  <- readRDS(file.path(BASE, "cache/v3/deconv/bulk_tpm.rds"))
LOG2 <- readRDS(file.path(BASE, "cache/v3/deconv/bulk_log2.rds"))
SS   <- colnames(TPM)

## ------------------------------------------------- (1) collagen-free new methods ----
SIGNC <- as.matrix(read.delim(file.path(BASE, "cache/v3/deconv/wu2021_signature_matrix_nocollagen.txt"),
                              row.names = 1, check.names = FALSE))
logf("collagen-free Wu signature: ", nrow(SIGNC), " genes x ", ncol(SIGNC), " types; ",
     "COL* genes present = ", sum(grepl("^COL", rownames(SIGNC))))
gi <- intersect(rownames(SIGNC), rownames(TPM))
S  <- SIGNC[gi, ]; Y <- TPM[gi, ]; K <- ncol(S)
logf("genes used = ", length(gi))
## the collagen-free variants are deterministic; skip recomputation if already cached
SKIP <- all(file.exists(file.path(OUT, paste0(c("QPROG_Wu_nocol", "RLS_Wu_nocol",
          "NNLS_Wu_nocol", "DTANGLE_Wu_nocol"), ".rds"))))
logf("collagen-free variants already cached, skipping recomputation: ", SKIP)
if (!SKIP) {

sc <- max(S); S2 <- S / sc; Y2 <- Y / sc
qp <- t(simplify2array(lapply(seq_len(ncol(Y2)), function(j) {
  s <- tryCatch(limSolve::lsei(A = S2, B = Y2[, j], E = matrix(1, 1, K), F = 1,
                               G = diag(K), H = rep(0, K), type = 2)$X,
                error = function(e) rep(NA_real_, K)); s / sum(s) })))
dimnames(qp) <- list(SS, colnames(S)); saveRDS(qp, file.path(OUT, "QPROG_Wu_nocol.rds"))
logf("QPROG_Wu_nocol done; failures = ", sum(!complete.cases(qp)))

rl <- t(simplify2array(lapply(seq_len(ncol(Y)), function(j) {
  co <- tryCatch(coef(MASS::rlm(S, Y[, j], maxit = 100)), error = function(e) rep(NA_real_, K))
  co[!is.finite(co) | co < 0] <- 0
  if (sum(co) <= 0) rep(NA_real_, K) else co / sum(co) })))
dimnames(rl) <- list(SS, colnames(S)); saveRDS(rl, file.path(OUT, "RLS_Wu_nocol.rds"))
logf("RLS_Wu_nocol done; failures = ", sum(!complete.cases(rl)))

nn <- t(simplify2array(lapply(seq_len(ncol(Y)), function(j) {
  w <- nnls::nnls(S, Y[, j])$x; w / sum(w) })))
dimnames(nn) <- list(SS, colnames(S)); saveRDS(nn, file.path(OUT, "NNLS_Wu_nocol.rds"))
logf("NNLS_Wu_nocol done")

## dtangle on the full Wu reference with collagens removed from the marker pool
REFWU <- as.matrix(read.delim(file.path(BASE, "cache/v3/deconv/wu2021_meancpm_by_celltype.tsv"),
                              row.names = 1, check.names = FALSE))
{
  library(dtangle)
  ref <- log2(REFWU + 1)
  ref <- ref[!grepl("^COL[0-9]", rownames(ref)), , drop = FALSE]
  g2 <- intersect(rownames(ref), rownames(LOG2))
  Yd <- t(LOG2[g2, , drop = FALSE]); Rd <- t(ref[g2, , drop = FALSE])
  ps <- lapply(seq_len(nrow(Rd)), function(i) nrow(Yd) + i)
  mk <- dtangle::find_markers(Y = rbind(Yd, Rd), pure_samples = ps, data_type = "rna-seq",
                              marker_method = "ratio")
  nmk <- unlist(lapply(mk$L, function(v) max(5L, floor(length(v) * 0.05))))
  dd <- dtangle::dtangle(Y = rbind(Yd, Rd), pure_samples = ps, n_markers = nmk,
                         markers = mk$L, data_type = "rna-seq")
  P <- dd$estimates[seq_len(nrow(Yd)), , drop = FALSE]
  dimnames(P) <- list(rownames(Yd), colnames(REFWU))
  saveRDS(P, file.path(OUT, "DTANGLE_Wu_nocol.rds"))
  logf("DTANGLE_Wu_nocol done; genes = ", length(g2))
  ## record which genes actually became CAF markers in the collagen-containing run
  mk0 <- dtangle::find_markers(Y = rbind(t(LOG2[intersect(rownames(REFWU), rownames(LOG2)), ]),
                                          t(log2(REFWU[intersect(rownames(REFWU), rownames(LOG2)), ] + 1))),
                               pure_samples = ps, data_type = "rna-seq", marker_method = "ratio")
  cafidx <- which(colnames(REFWU) == "CAFs")
  gAll <- intersect(rownames(REFWU), rownames(LOG2))
  ntop <- max(5L, floor(length(mk0$L[[cafidx]]) * 0.05))
  cafmk <- gAll[mk0$L[[cafidx]][seq_len(ntop)]]
  logf("dtangle CAF markers used (n = ", length(cafmk), "); COL* among them: ",
       paste(grep("^COL", cafmk, value = TRUE), collapse = ", "))
  writeLines(cafmk, file.path(BASE, "cache/v3/deconv2/dtangle_CAF_markers.txt"))
}
}   # end if (!SKIP)

## --------------------------------------------------------- (2) assemble everything ----
FIB0 <- readRDS(file.path(BASE, "cache/v3/deconv/fibroblast_estimates_all.rds"))
stopifnot(identical(rownames(FIB0), SS))
logf("\nfirst-pass estimates: ", ncol(FIB0))

grab <- function(f, col, nm) {
  M <- readRDS(file.path(OUT, paste0(f, ".rds")))
  stopifnot(identical(rownames(M), SS))
  v <- M[, col]; names(v) <- SS
  setNames(list(v), nm)
}
newl <- c(
  grab("QPROG_Wu",         "CAFs", "QPROG_Wu_CAFs"),
  grab("RLS_Wu",           "CAFs", "RLS_Wu_CAFs"),
  grab("OLS_Wu",           "CAFs", "OLS_Wu_CAFs"),
  grab("DWLS_Wu",          "CAFs", "DWLS_Wu_CAFs"),
  grab("DTANGLE_Wu",       "CAFs", "DTANGLE_Wu_CAFs"),
  grab("QPROG_Wu_nocol",   "CAFs", "QPROG_Wu_nocol_CAFs"),
  grab("RLS_Wu_nocol",     "CAFs", "RLS_Wu_nocol_CAFs"),
  grab("NNLS_Wu_nocol",    "CAFs", "NNLS_Wu_nocol_CAFs"),
  grab("DTANGLE_Wu_nocol", "CAFs", "DTANGLE_Wu_nocol_CAFs"),
  ## strictly matched collagen controls from script 50 (same 641-gene signature,
  ## only the collagens deleted) -- these isolate the circularity of the signature
  grab("NNLS_Wu641_full",         "CAFs", "NNLS_Wu641_full_CAFs"),
  grab("NNLS_Wu641_minusAllCOL",  "CAFs", "NNLS_Wu641_minusAllCOL_CAFs"),
  grab("QPROG_Wu641_full",        "CAFs", "QPROG_Wu641_full_CAFs"),
  grab("QPROG_Wu641_minusAllCOL", "CAFs", "QPROG_Wu641_minusAllCOL_CAFs"))
NEW <- do.call(cbind, newl)
logf("new algorithmic estimates: ", ncol(NEW), " -> ", paste(colnames(NEW), collapse = ", "))

## ---- non-transcriptomic mediators ----
EXT <- readRDS(file.path(BASE, "cache/v3/deconv2/external_estimates.rds"))
pur <- EXT$absolute_purity[SS]
nontum <- 1 - pur
logf("ABSOLUTE non-tumour fraction: n = ", sum(!is.na(nontum)), "; median = ",
     round(median(nontum, na.rm = TRUE), 3))

PA <- fread(file.path(BASE, "data/deconv/tcga_brca_slide_pathology.tsv"))
PA[, short := substr(sample, 1, 15)]
PA[, kind := fifelse(grepl("-BS", slide), "BS", fifelse(grepl("-TS", slide), "TS",
              fifelse(grepl("-MS", slide), "MS", "DX")))]
frozen <- PA[kind %in% c("BS", "TS", "MS") & !is.na(percent_stromal_cells),
             .(v = mean(percent_stromal_cells)), by = short]
allsl  <- PA[!is.na(percent_stromal_cells), .(v = mean(percent_stromal_cells)), by = short]
lymph  <- PA[kind %in% c("BS", "TS", "MS") & !is.na(percent_lymphocyte_infiltration),
             .(v = mean(percent_lymphocyte_infiltration)), by = short]
tumnuc <- PA[kind %in% c("BS", "TS", "MS") & !is.na(percent_tumor_nuclei),
             .(v = mean(percent_tumor_nuclei)), by = short]
path_frozen <- setNames(frozen$v[match(SS, frozen$short)], SS)
path_all    <- setNames(allsl$v[match(SS, allsl$short)], SS)
path_lymph  <- setNames(lymph$v[match(SS, lymph$short)], SS)
path_tumnuc <- setNames(tumnuc$v[match(SS, tumnuc$short)], SS)
logf("pathology percent_stromal_cells (frozen BS/TS/MS slides): n = ", sum(!is.na(path_frozen)),
     "; median = ", round(median(path_frozen, na.rm = TRUE), 1), "%")
logf("pathology percent_stromal_cells (all slides): n = ", sum(!is.na(path_all)),
     "; median = ", round(median(path_all, na.rm = TRUE), 1), "%")

## DNA-methylation-based fibroblast fraction (EpiDISH, script 52) -- a different ASSAY
EPI <- tryCatch(readRDS(file.path(BASE, "cache/v3/deconv2/epidish_fractions.rds")),
                error = function(e) NULL)
if (!is.null(EPI)) {
  stopifnot(identical(rownames(EPI), SS))
  logf("EpiDISH 450k fractions available for ", sum(!is.na(EPI[, "Fib"])), " tumours")
} else logf("EpiDISH fractions not available yet")

ORTH <- cbind(ABSOLUTE_nontumour = nontum,
              PATH_stroma_frozen = path_frozen,
              PATH_stroma_allslides = path_all)
if (!is.null(EPI)) ORTH <- cbind(ORTH, METH_EpiDISH_Fib = EPI[, "Fib"],
                                       METH_EpiDISH_Epi = EPI[, "Epi"],
                                       METH_EpiDISH_IC  = EPI[, "IC"])
logf("orthogonal (non-RNA) mediators: ", ncol(ORTH))

ALL <- cbind(FIB0, NEW, ORTH)
saveRDS(ALL, file.path(BASE, "cache/v3/deconv2/fibroblast_estimates_extended.rds"))
logf("EXTENDED matrix: ", nrow(ALL), " samples x ", ncol(ALL), " estimates")
fwrite(data.table(sample = rownames(ALL), as.data.frame(ALL)),
       file.path(BASE, "results/v3/deconv_fibroblast_estimates_extended.csv"))

## extra covariates used later
saveRDS(list(path_lymph = path_lymph, path_tumnuc = path_tumnuc, purity = pur),
        file.path(BASE, "cache/v3/deconv2/extra_covariates.rds"))

## ------------------------------------------------------- (3) correlation matrices ----
for (mth in c("spearman", "pearson")) {
  C <- cor(ALL, method = mth, use = "pairwise.complete.obs")
  fwrite(data.table(estimate = rownames(C), as.data.frame(C)),
         file.path(BASE, paste0("results/v3/deconv_fibcorr_extended_", mth, ".csv")))
  ut <- C[upper.tri(C)]
  logf("\n", toupper(mth), " across ", ncol(C), " estimates: min ", round(min(ut, na.rm = TRUE), 3),
       " median ", round(median(ut, na.rm = TRUE), 3), " max ", round(max(ut, na.rm = TRUE), 3))
}
C <- cor(ALL, method = "spearman", use = "pairwise.complete.obs")
n_ok <- sapply(seq_len(ncol(ALL)), function(i) sum(!is.na(ALL[, i])))
pairs <- rbindlist(lapply(seq_len(ncol(C) - 1), function(i) {
  data.table(a = colnames(C)[i], b = colnames(C)[(i + 1):ncol(C)], rho = C[i, (i + 1):ncol(C)])
}))
fwrite(pairs[order(rho)], file.path(BASE, "results/v3/deconv_fibcorr_extended_pairs.csv"))
logf("\n--- 15 LOWEST pairwise Spearman correlations among all stromal estimates ---")
p15 <- pairs[order(rho)][1:15]
for (i in 1:15) logf(sprintf("  %+0.3f   %s  vs  %s", p15$rho[i], p15$a[i], p15$b[i]))
logf("\n--- correlations of the two NON-RNA estimates with the RNA-based ones ---")
for (o in c("ABSOLUTE_nontumour", "PATH_stroma_frozen", "PATH_stroma_allslides")) {
  v <- pairs[a == o | b == o]
  logf(sprintf("  %-24s n = %4d  rho range %+0.3f .. %+0.3f  median %+0.3f",
               o, n_ok[match(o, colnames(ALL))], min(v$rho), max(v$rho), median(v$rho)))
}
logf("DONE 45")
