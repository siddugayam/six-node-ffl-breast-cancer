#!/usr/bin/env Rscript
## 53_deconv2_consolidate.R -- one compact, manuscript-ready summary of the multi-method
## deconvolution re-test. Every number here is copied from a file written by scripts
## 41-52; nothing is recomputed by hand.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/53_deconv2_consolidate.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
add <- function(section, item, value) data.table(section = section, item = item, value = as.character(value))
OUT <- list()

FIB <- readRDS(file.path(BASE, "cache/v3/deconv2/fibroblast_estimates_extended.rds"))
C <- cor(FIB, method = "spearman", use = "pairwise.complete.obs")
nonfib <- c("METH_EpiDISH_Epi", "METH_EpiDISH_IC")
fib <- setdiff(colnames(FIB), nonfib)
rna <- setdiff(fib, c("ABSOLUTE_nontumour", "PATH_stroma_frozen", "PATH_stroma_allslides", "METH_EpiDISH_Fib"))
u_all <- C[fib, fib][upper.tri(diag(length(fib)))]
u_rna <- C[rna, rna][upper.tri(diag(length(rna)))]
OUT <- c(OUT, list(
  add("A methods", "stromal/fibroblast estimates assembled", length(fib)),
  add("A methods", "of which RNA-based", length(rna)),
  add("A methods", "non-RNA (DNA purity, pathology, DNA methylation)", length(fib) - length(rna)),
  add("B correlation", "pairwise Spearman, all fibroblast estimates: min",  round(min(u_all), 3)),
  add("B correlation", "pairwise Spearman, all fibroblast estimates: median", round(median(u_all), 3)),
  add("B correlation", "pairwise Spearman, RNA-based only: min", round(min(u_rna), 3)),
  add("B correlation", "pairwise Spearman, RNA-based only: median", round(median(u_rna), 3))))

R <- fread(file.path(BASE, "results/v3/deconv_mediation_extended.csv"))
S <- fread(file.path(BASE, "results/v3/deconv_mediation_extended_summary.csv"))
for (i in seq_len(nrow(S))) { r <- S[i]
  ax <- paste(r$x, "->", r$y)
  OUT <- c(OUT, list(
    add("C mediation", paste(ax, "| n mediators"), r$n_methods),
    add("C mediation", paste(ax, "| ACME median [min,max]"),
        sprintf("%+.3f [%+.3f, %+.3f]", r$ACME_med, r$ACME_min, r$ACME_max)),
    add("C mediation", paste(ax, "| ADE median [min,max]"),
        sprintf("%+.3f [%+.3f, %+.3f]", r$ADE_med, r$ADE_min, r$ADE_max)),
    add("C mediation", paste(ax, "| proportion mediated median"),
        sprintf("%+.2f", r$propmed_med)),
    add("C mediation", paste(ax, "| Sobel p<0.05"), sprintf("%d/%d", r$n_sobel_sig, r$n_methods)),
    add("C mediation", paste(ax, "| sign-inconsistent (ACME vs ADE)"),
        sprintf("%d/%d", r$n_inconsistent, r$n_methods))))
}
CF <- fread(file.path(BASE, "results/v3/deconv_mediation_collagenfree_only.csv"))
for (i in seq_len(nrow(CF))) { r <- CF[i]
  OUT <- c(OUT, list(add("C mediation, collagen-free mediators only",
    paste(r$axis, "| ACME median, Sobel sig, inconsistent"),
    sprintf("%+.3f ; %d/%d ; %d/%d", r$ACME_med, r$n_sobel_sig, r$n, r$n_inconsistent, r$n))))
}

H <- fread(file.path(BASE, "results/v3/deconv_compartment_target_summary.csv"))
for (i in seq_len(nrow(H))) { r <- H[i]
  if (!r$compartment %in% c("Fibroblast/CAF", "Malignant/other", "Myeloid", "T/NK")) next
  OUT <- c(OUT, list(add("D compartment", paste(r$target, "|", r$compartment),
    sprintf("median rho %+.3f over %d features [%+.3f,%+.3f]; %d sig+, %d sig-",
            r$median_rho, r$n_features, r$min_rho, r$max_rho, r$n_pos_sig, r$n_neg_sig))))
}

AC <- fread(file.path(BASE, "results/v3/deconv_acme_vs_accuracy.csv"))
for (i in seq_len(nrow(AC))) { r <- AC[i]
  OUT <- c(OUT, list(add("E why methods differ", paste(r$axis, "| rho(simulated-mixture accuracy, |ACME|)"),
    sprintf("%+.3f (n=%d mediators, p=%.3g)", r$spearman_acc_vs_absACME, r$n_mediators_with_acc, r$p))))
}
CI <- fread(file.path(BASE, "results/v3/deconv_acme_vs_circularity.csv"))
for (i in seq_len(nrow(CI))) { r <- CI[i]
  OUT <- c(OUT, list(add("E why methods differ", paste(r$axis, "| rho(ACME, rho[mediator,COL1A1])"),
    sprintf("%+.3f ; median ACME collagen-free %+.3f vs collagen-containing %+.3f (Wilcoxon p=%.3g)",
            r$spearman_ACME_vs_rho_med_COL1A1, r$ACME_collagenfree_median,
            r$ACME_collagencontaining_median, r$wilcox_p))))
}

EXTA <- fread(file.path(BASE, "results/v3/deconv_external_reference_agreement.csv"))
OUT <- c(OUT, list(
  add("E external checks", "in-house nu-SVR CIBERSORT vs official Thorsson CIBERSORT (22 LM22 types, n=1095)",
      sprintf("median Spearman %.3f; abundance-weighted %.3f; %d/22 types above 0.6",
              median(EXTA$spearman), sum(EXTA$spearman * EXTA$mean_thorsson) / sum(EXTA$mean_thorsson),
              sum(EXTA$spearman > 0.6)))))
MA <- fread(file.path(BASE, "results/v3/deconv_methylation_vs_rna_agreement.csv"))
OUT <- c(OUT, list(
  add("E external checks", "DNA-methylation (EpiDISH) fibroblast fraction vs RNA-based estimates",
      sprintf("median Spearman %.3f (range %.3f to %.3f, n=783 tumours)",
              median(MA$rho), min(MA$rho), max(MA$rho)))))
BN <- tryCatch(fread(file.path(BASE, "results/v3/deconv_benchmark_newmethods.csv")), error = function(e) NULL)
if (!is.null(BN)) for (i in seq_len(nrow(BN))) OUT <- c(OUT, list(
  add("A benchmark (held-out simulated mixtures, n=300)", BN$method[i],
      sprintf("Spearman vs true CAF fraction %.3f", BN$spearman_vs_true_CAF[i]))))
BP <- tryCatch(fread(file.path(BASE, "results/v3/deconv_bayesprism_vs_instaprism.csv")), error = function(e) NULL)
if (!is.null(BP)) OUT <- c(OUT, list(
  add("E external checks", "BayesPrism (Gibbs) vs InstaPrism (analytic), same Wu reference",
      sprintf("median Spearman across cell types %.3f; CAFs %.3f",
              median(BP$spearman), BP[cell_type == "CAFs"]$spearman))))

TAB <- rbindlist(OUT)
fwrite(TAB, file.path(BASE, "results/v3/deconv_FINAL_summary.csv"))
logf("WROTE results/v3/deconv_FINAL_summary.csv  rows = ", nrow(TAB))
for (i in seq_len(nrow(TAB))) logf(sprintf("[%-42s] %-62s %s", TAB$section[i], TAB$item[i], TAB$value[i]))
logf("DONE 53")
