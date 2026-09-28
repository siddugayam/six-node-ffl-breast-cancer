#!/usr/bin/env Rscript
## 48_deconv2_celltype_correlations.R -- part (D).
## Correlate EVERY cell-type fraction from EVERY method with COL1A1, COL3A1, the miR-29
## family and the two TFs, so the paper can state which compartment each component of
## the circuit tracks. Also correlates against the two non-transcriptomic composition
## measures (ABSOLUTE purity, pathologist stroma score) for calibration.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/48_deconv2_celltype_corr.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

gexp <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mexp <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
stopifnot(!any(duplicated(rownames(gexp))))
SS <- readLines(file.path(BASE, "cache/v3/deconv/samples_primary_tumour.txt"))

D1 <- file.path(BASE, "cache/v3/deconv/out")
D2 <- file.path(BASE, "cache/v3/deconv2/out")
mats <- list()
for (f in list.files(D1, pattern = "\\.rds$", full.names = TRUE)) {
  nm <- sub("\\.rds$", "", basename(f))
  if (grepl("_stats$", nm)) next
  M <- readRDS(f)
  if (!is.matrix(M) || is.null(rownames(M)) || !all(SS %in% rownames(M))) next
  mats[[nm]] <- M[SS, , drop = FALSE]
}
for (f in list.files(D2, pattern = "\\.rds$", full.names = TRUE)) {
  nm <- sub("\\.rds$", "", basename(f))
  M <- readRDS(f)
  if (!is.matrix(M) || is.null(rownames(M)) || !all(SS %in% rownames(M))) {
    logf("skipping ", nm, " (not on the full sample set)"); next }
  mats[[nm]] <- M[SS, , drop = FALSE]
}
EXT <- readRDS(file.path(BASE, "cache/v3/deconv2/external_estimates.rds"))
TH <- EXT$thorsson_cibersort
TH2 <- matrix(NA_real_, length(SS), ncol(TH), dimnames = list(SS, colnames(TH)))
TH2[rownames(TH), ] <- TH
mats[["THORSSON_CIBERSORT"]] <- TH2
EP <- tryCatch(readRDS(file.path(BASE, "cache/v3/deconv2/epidish_fractions.rds")), error = function(e) NULL)
if (!is.null(EP)) { colnames(EP) <- paste0("EpiDISH_", colnames(EP)); mats[["METH_EpiDISH"]] <- EP[SS, , drop = FALSE] }
logf("method matrices loaded: ", length(mats), " -> ", paste(names(mats), collapse = ", "))
logf("total cell-type features: ", sum(sapply(mats, ncol)))

## targets
mir29 <- mexp[c("hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"), , drop = FALSE]
mir29fam <- colMeans(mir29)
tg <- list(COL1A1 = gexp["COL1A1", ], COL3A1 = gexp["COL3A1", ],
           ETS1 = gexp["ETS1", ], NFKB1 = gexp["NFKB1", ],
           `hsa-miR-29a` = mexp["hsa-miR-29a", ], `hsa-miR-29b` = mexp["hsa-miR-29b", ],
           `hsa-miR-29c` = mexp["hsa-miR-29c", ], miR29_family_mean = mir29fam)
cov2 <- readRDS(file.path(BASE, "cache/v3/deconv2/extra_covariates.rds"))
tg$ABSOLUTE_purity <- cov2$purity
tg$PATH_lymphocytes <- cov2$path_lymph
tg$PATH_tumour_nuclei <- cov2$path_tumnuc

res <- rbindlist(lapply(names(mats), function(mn) {
  M <- mats[[mn]]
  rbindlist(lapply(colnames(M), function(ct) {
    v <- M[, ct]
    rbindlist(lapply(names(tg), function(tn) {
      t0 <- tg[[tn]]
      ss <- intersect(SS, names(t0))
      a <- v[ss]; b <- as.numeric(t0[ss])
      ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 50 || sd(a[ok]) == 0 || sd(b[ok]) == 0)
        return(data.table(method = mn, cell_type = ct, target = tn, n = sum(ok),
                          rho = NA_real_, p = NA_real_))
      ct1 <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = FALSE))
      data.table(method = mn, cell_type = ct, target = tn, n = sum(ok),
                 rho = unname(ct1$estimate), p = ct1$p.value)
    }))
  }))
}))
res[, fdr := p.adjust(p, "BH"), by = target]
fwrite(res, file.path(BASE, "results/v3/deconv_celltype_target_correlations_extended.csv"))
logf("WROTE results/v3/deconv_celltype_target_correlations_extended.csv rows = ", nrow(res))

show <- function(tn, k = 12) {
  s <- res[target == tn & !is.na(rho)][order(-rho)]
  logf("\n=== ", tn, " : ", k, " most POSITIVELY correlated cell-type estimates ===")
  for (i in seq_len(min(k, nrow(s)))) logf(sprintf("  %+0.3f  (FDR %8.2g)  %s :: %s",
    s$rho[i], s$fdr[i], s$method[i], s$cell_type[i]))
  logf("--- ", k, " most NEGATIVELY correlated ---")
  s2 <- s[order(rho)]
  for (i in seq_len(min(k, nrow(s2)))) logf(sprintf("  %+0.3f  (FDR %8.2g)  %s :: %s",
    s2$rho[i], s2$fdr[i], s2$method[i], s2$cell_type[i]))
}
for (tn in c("COL1A1", "COL3A1", "miR29_family_mean", "hsa-miR-29a", "ETS1", "NFKB1")) show(tn)

## harmonised view: median rho per cell-type LABEL across methods
lab <- function(x) {
  x <- tolower(x)
  fifelse(grepl("caf|fibroblast", x), "Fibroblast/CAF",
   fifelse(grepl("endothel", x), "Endothelial",
    fifelse(grepl("pvl|pericyt", x), "Perivascular",
     fifelse(grepl("cancer epithel|tumou?r|uncharacter|othercell|other cell", x), "Malignant/other",
      fifelse(grepl("normal epithel", x), "Normal epithelial",
       fifelse(grepl("macrophage|monocyt|myeloid|dendritic|neutrophil|eosinophil|mast|granulocyt", x), "Myeloid",
        fifelse(grepl("t cell|t-cell|tregs|cd8|cd4|th1|th2|tgd|nk", x), "T/NK",
         fifelse(grepl("b cell|b-cell|plasma", x), "B/plasma", "other"))))))))
}
res[, compartment := lab(cell_type)]
H <- res[!is.na(rho) & target %in% c("COL1A1", "COL3A1", "miR29_family_mean"),
         .(n_features = .N, median_rho = median(rho), min_rho = min(rho), max_rho = max(rho),
           n_pos_sig = sum(rho > 0 & fdr < 0.05), n_neg_sig = sum(rho < 0 & fdr < 0.05)),
         by = .(target, compartment)][order(target, -median_rho)]
fwrite(H, file.path(BASE, "results/v3/deconv_compartment_target_summary.csv"))
logf("\n=== COMPARTMENT-LEVEL SUMMARY (median Spearman across all methods' features) ===")
for (i in seq_len(nrow(H))) { r <- H[i]
  logf(sprintf("%-20s %-18s n=%3d  median rho %+0.3f  [%+0.3f,%+0.3f]  sig+ %d  sig- %d",
    r$target, r$compartment, r$n_features, r$median_rho, r$min_rho, r$max_rho, r$n_pos_sig, r$n_neg_sig)) }
logf("DONE 48")
