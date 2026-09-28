#!/usr/bin/env Rscript
## 51_deconv2_contradictions.R -- part (E).
## Which methods disagree, and why. For every mediation axis, classify each mediator's
## verdict, find the minority verdicts, and test whether the disagreement is explained by
##   (i)  how strongly the mediator itself tracks the outcome collagen (circularity),
##   (ii) whether the mediator's signature contains collagen genes,
##   (iii) the mediator's measured accuracy on held-out simulated mixtures,
##   (iv) whether the mediator is transcriptomic at all.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/51_deconv2_contradictions.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

R <- fread(file.path(BASE, "results/v3/deconv_mediation_extended.csv"))
FIB <- readRDS(file.path(BASE, "cache/v3/deconv2/fibroblast_estimates_extended.rds"))
gexp <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
SS <- rownames(FIB)

meta <- data.table(mediator = colnames(FIB))
meta[, rho_med_COL1A1 := sapply(mediator, function(m)
  cor(FIB[, m], gexp["COL1A1", SS], method = "spearman", use = "complete.obs"))]
meta[, rho_med_COL3A1 := sapply(mediator, function(m)
  cor(FIB[, m], gexp["COL3A1", SS], method = "spearman", use = "complete.obs"))]
NONRNA <- c("ABSOLUTE_nontumour", "PATH_stroma_frozen", "PATH_stroma_allslides",
            "METH_EpiDISH_Fib", "METH_EpiDISH_Epi", "METH_EpiDISH_IC")
meta[, transcriptomic := !mediator %in% NONRNA]
## collagen status of each mediator's signature, as far as it is knowable
colfree <- c("xCell_Fibroblasts", "ssGSEA_CAF_scRNA", "meanZ_CAF_scRNA", "meanZ_CAF_B_markers",
             "CBSX_Wu_nocol_CAFs", "InstaPrism_Wu_nocollagen_CAFs", "MCPcounter_Fibroblasts_nocol",
             "EPIC_CAFs_nocol", "ConsensusTME_Fibroblasts_nocol", "QPROG_Wu_nocol_CAFs",
             "RLS_Wu_nocol_CAFs", "NNLS_Wu_nocol_CAFs", "DTANGLE_Wu_nocol_CAFs",
             "NNLS_Wu641_minusAllCOL_CAFs", "QPROG_Wu641_minusAllCOL_CAFs",
             "meanZ_FARMER_clean", "ssGSEA_ESTIMATE_stromal_clean", "CONSENSUS_rankmean_collagensafe",
             "ABSOLUTE_nontumour", "PATH_stroma_frozen", "PATH_stroma_allslides",
             "METH_EpiDISH_Fib", "METH_EpiDISH_Epi", "METH_EpiDISH_IC")
meta[, signature_collagen_free := mediator %in% colfree]

## mediator FAMILY: a compositional fraction (sums to 1 across cell types) vs an
## enrichment / marker score with no compositional constraint vs a non-RNA measure
fracpat <- "^(CBSX_Wu|NNLS_Wu|EPIC_CAFs|InstaPrism|QPROG_|RLS_|OLS_|DWLS_|DTANGLE_)"
meta[, family := fifelse(!transcriptomic, "non-RNA",
                  fifelse(grepl(fracpat, mediator), "compositional fraction",
                          "enrichment / marker score"))]

## measured accuracy on held-out simulated mixtures, where available
acc <- rbindlist(list(
  tryCatch(fread(file.path(BASE, "results/v3/deconv_simulation_benchmark.csv"))[
             , .(mediator = estimate, acc = spearman_vs_true_CAF)], error = function(e) NULL),
  tryCatch(fread(file.path(BASE, "results/v3/deconv_benchmark_newmethods.csv"))[
             , .(mediator = paste0(method, "_Wu_CAFs"), acc = spearman_vs_true_CAF)],
           error = function(e) NULL)), fill = TRUE)
acc <- acc[!duplicated(mediator)]
meta <- merge(meta, acc, by = "mediator", all.x = TRUE)
typ <- tryCatch(fread(file.path(BASE, "results/v3/deconv_fibroblast_typicality.csv"))[
                  , .(mediator = estimate, typicality = mean_rho_with_others)],
                error = function(e) NULL)
if (!is.null(typ)) meta <- merge(meta, typ, by = "mediator", all.x = TRUE)
fwrite(meta, file.path(BASE, "results/v3/deconv_mediator_metadata.csv"))
logf("mediators with a simulated-mixture accuracy value: ", sum(!is.na(meta$acc)))
logf("mediators: ", nrow(meta), "; declared collagen-free: ", sum(meta$signature_collagen_free),
     "; non-transcriptomic: ", sum(!meta$transcriptomic))

R <- merge(R, meta, by = "mediator")
R[, axis := paste(x, y, sep = " -> ")]
R[, verdict := fifelse(sobel_p >= 0.05, "no mediation",
                fifelse(inconsistent, "INCONSISTENT",
                 fifelse(ade_p >= 0.05, "complete mediation", "partial same-signed")))]

logf("\n=== VERDICT TABLE BY AXIS ===")
V <- dcast(R[, .N, by = .(axis, verdict)], axis ~ verdict, value.var = "N", fill = 0)
print(V); cat(capture.output(print(V)), sep = "\n", file = LOG, append = TRUE)
fwrite(V, file.path(BASE, "results/v3/deconv_verdict_table_by_axis.csv"))

logf("\n=== MINORITY VERDICTS (the methods that disagree) ===")
minor <- rbindlist(lapply(unique(R$axis), function(a) {
  s <- R[axis == a]
  tb <- sort(table(s$verdict), decreasing = TRUE)
  maj <- names(tb)[1]
  s[verdict != maj, .(axis = a, majority = maj, mediator, verdict, n, ACME, ADE,
                      prop_mediated, sobel_p, rho_raw, rho_partial,
                      rho_med_COL1A1, transcriptomic, signature_collagen_free)]
}))
fwrite(minor, file.path(BASE, "results/v3/deconv_minority_verdicts.csv"))
for (a in unique(minor$axis)) {
  s <- minor[axis == a]
  logf("\n--- ", a, "  (majority = ", s$majority[1], "; ", nrow(s), " dissenting mediators) ---")
  for (i in seq_len(nrow(s))) logf(sprintf("   %-32s %-20s ACME %+0.3f  ADE %+0.3f  sobel_p %8.2g  rho(med,COL1A1) %+0.2f  %s",
    s$mediator[i], s$verdict[i], s$ACME[i], s$ADE[i], s$sobel_p[i], s$rho_med_COL1A1[i],
    ifelse(s$transcriptomic[i], "", "[NON-RNA]")))
}

## ---------------- is the size of the mediated effect explained by circularity? -------
logf("\n=== ACME vs how strongly the mediator itself tracks the outcome ===")
fitrows <- rbindlist(lapply(unique(R$axis), function(a) {
  s <- R[axis == a]
  rr <- cor(s$ACME, s$rho_med_COL1A1, method = "spearman")
  f <- lm(ACME ~ rho_med_COL1A1, data = s)
  data.table(axis = a, n_mediators = nrow(s),
             spearman_ACME_vs_rho_med_COL1A1 = rr,
             slope = coef(f)[2], slope_p = summary(f)$coef[2, 4], r2 = summary(f)$r.squared,
             ACME_collagenfree_median = median(s[signature_collagen_free == TRUE]$ACME),
             ACME_collagencontaining_median = median(s[signature_collagen_free == FALSE]$ACME),
             wilcox_p = tryCatch(wilcox.test(ACME ~ signature_collagen_free, data = s)$p.value,
                                 error = function(e) NA_real_))
}))
fwrite(fitrows, file.path(BASE, "results/v3/deconv_acme_vs_circularity.csv"))
for (i in seq_len(nrow(fitrows))) { r <- fitrows[i]
  logf(sprintf("%-24s rho(ACME, rho_med_COL1A1) = %+0.3f | R2 = %0.3f | median ACME collagen-free %+0.3f vs collagen-containing %+0.3f (Wilcoxon p = %.3g)",
    r$axis, r$spearman_ACME_vs_rho_med_COL1A1, r$r2, r$ACME_collagenfree_median,
    r$ACME_collagencontaining_median, r$wilcox_p)) }

## ---------------- restrict to the defensible mediators only --------------------------
logf("\n=== HEADLINE RESTRICTED TO COLLAGEN-FREE MEDIATORS ONLY ===")
S2 <- R[signature_collagen_free == TRUE,
        .(n = .N, ACME_med = median(ACME), ACME_min = min(ACME), ACME_max = max(ACME),
          ADE_med = median(ADE), ADE_min = min(ADE), ADE_max = max(ADE),
          n_sobel_sig = sum(sobel_p < 0.05), n_inconsistent = sum(inconsistent),
          propmed_med = median(prop_mediated)), by = axis]
fwrite(S2, file.path(BASE, "results/v3/deconv_mediation_collagenfree_only.csv"))
for (i in seq_len(nrow(S2))) { r <- S2[i]
  logf(sprintf("%-24s n=%2d  ACME %+0.3f [%+0.3f,%+0.3f]  ADE %+0.3f [%+0.3f,%+0.3f]  propmed %+0.2f  Sobel sig %d/%d  INCONSISTENT %d/%d",
    r$axis, r$n, r$ACME_med, r$ACME_min, r$ACME_max, r$ADE_med, r$ADE_min, r$ADE_max,
    r$propmed_med, r$n_sobel_sig, r$n, r$n_inconsistent, r$n)) }

## ------------- does the disagreement track mediator FAMILY or mediator ACCURACY? ----
logf("\n=== ACME BY MEDIATOR FAMILY ===")
FAM <- R[, .(n = .N, ACME_med = median(ACME), ACME_min = min(ACME), ACME_max = max(ACME),
             ADE_med = median(ADE), propmed_med = median(prop_mediated),
             n_sobel_sig = sum(sobel_p < 0.05), n_inconsistent = sum(inconsistent)),
         by = .(axis, family)][order(axis, family)]
fwrite(FAM, file.path(BASE, "results/v3/deconv_acme_by_mediator_family.csv"))
for (i in seq_len(nrow(FAM))) { r <- FAM[i]
  logf(sprintf("%-24s %-26s n=%2d  ACME %+0.3f [%+0.3f,%+0.3f]  ADE %+0.3f  propmed %+0.2f  sig %d  incons %d",
    r$axis, r$family, r$n, r$ACME_med, r$ACME_min, r$ACME_max, r$ADE_med, r$propmed_med,
    r$n_sobel_sig, r$n_inconsistent)) }

logf("\n=== ACME vs MEASURED ACCURACY of the mediator (held-out simulated mixtures) ===")
ACC <- rbindlist(lapply(unique(R$axis), function(a) {
  s <- R[axis == a & !is.na(acc)]
  if (nrow(s) < 5) return(NULL)
  ct1 <- suppressWarnings(cor.test(s$acc, s$ACME, method = "spearman", exact = FALSE))
  data.table(axis = a, n_mediators_with_acc = nrow(s),
             spearman_acc_vs_ACME = unname(ct1$estimate), p = ct1$p.value,
             spearman_acc_vs_absACME = cor(s$acc, abs(s$ACME), method = "spearman"))
}))
fwrite(ACC, file.path(BASE, "results/v3/deconv_acme_vs_accuracy.csv"))
for (i in seq_len(nrow(ACC))) { r <- ACC[i]
  logf(sprintf("%-24s n=%2d  rho(accuracy, ACME) = %+0.3f (p = %.3g); rho(accuracy, |ACME|) = %+0.3f",
    r$axis, r$n_mediators_with_acc, r$spearman_acc_vs_ACME, r$p, r$spearman_acc_vs_absACME)) }

logf("\n=== THE TWO NON-TRANSCRIPTOMIC MEDIATORS, ALL AXES ===")
NT <- R[transcriptomic == FALSE][order(axis, mediator)]
for (i in seq_len(nrow(NT))) { r <- NT[i]
  logf(sprintf("%-24s %-22s n=%4d  rho %+0.3f  rho|M %+0.3f  ACME %+0.3f [%+0.3f,%+0.3f]  ADE %+0.3f  propmed %+0.3f  sobel_p %.2g  %s",
    r$axis, r$mediator, r$n, r$rho_raw, r$rho_partial, r$ACME, r$acme_lo, r$acme_hi,
    r$ADE, r$prop_mediated, r$sobel_p, r$verdict)) }
fwrite(NT, file.path(BASE, "results/v3/deconv_mediation_nonRNA_mediators.csv"))
logf("DONE 51")
