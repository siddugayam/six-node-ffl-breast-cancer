#!/usr/bin/env Rscript
## 47_neoadj_survival_and_agreement.R
##  (a) GSE25066 distant relapse-free survival (the only neoadjuvant cohort with outcome time)
##  (b) duplicate-patient detection across the GPL96 neoadjuvant series (they come from
##      overlapping MDACC/USO/LAB FNA collections, so the meta-analysis must be caveated)
##  (c) Task F: agreement between the Farmer signature and every other stromal / ECM /
##      wound / EMT signature, in TCGA-BRCA, METABRIC and the pooled neoadjuvant cohorts.
suppressPackageStartupMessages({library(data.table); library(survival); library(matrixStats)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v5/47_neoadj_survival_and_agreement.log")
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); set.seed(3)

store <- readRDS(file.path(BASE, "cache/v5/farmer_neoadjuvant_store.rds"))

## ------------------------------------------------------- (a) GSE25066 DRFS -------------
s <- store$GSE25066
ev <- suppressWarnings(as.numeric(sapply(s$meta, function(m) m[["drfs_1_event_0_censored"]])))
ti <- suppressWarnings(as.numeric(sapply(s$meta, function(m) m[["drfs_even_time_years"]])))
er <- s$er
logf("GSE25066 DRFS: evaluable n = ", sum(is.finite(ev) & is.finite(ti) & ti > 0),
     " | events = ", sum(ev == 1, na.rm = TRUE))
rows <- list(); k <- 0L
for (nm in rownames(s$SC)) {
  x <- as.numeric(scale(s$SC[nm, ]))
  ok <- is.finite(ev) & is.finite(ti) & ti > 0 & is.finite(x)
  if (sum(ok) < 100) next
  d1 <- data.frame(t = ti[ok], e = ev[ok], x = x[ok])
  f1 <- coxph(Surv(t, e) ~ x, data = d1)
  c1 <- summary(f1)$coefficients          # columns: coef, exp(coef), se(coef), z, Pr(>|z|)
  ci1 <- summary(f1)$conf.int
  d2 <- data.frame(t = ti[ok], e = ev[ok], x = x[ok], er = er[ok])
  d2 <- d2[complete.cases(d2), ]
  f2 <- coxph(Surv(t, e) ~ x + factor(er), data = d2); c2 <- summary(f2)$coefficients
  k <- k + 1L
  rows[[k]] <- data.table(cohort = "GSE25066", endpoint = "DRFS", score = nm, n = sum(ok),
    n_event = sum(ev[ok] == 1), HR_per_SD = c1["x", "exp(coef)"],
    CI_low = ci1["x", "lower .95"], CI_high = ci1["x", "upper .95"],
    p = c1["x", "Pr(>|z|)"], C_index = summary(f1)$concordance[1],
    HR_ER_adjusted = c2["x", "exp(coef)"], p_ER_adjusted = c2["x", "Pr(>|z|)"],
    n_ER_adjusted = nrow(d2))
}
drfs <- rbindlist(rows); drfs[, q_BH := p.adjust(p, "BH")]
setorder(drfs, p)
fwrite(drfs, file.path(BASE, "results/v5/farmer_gse25066_drfs.csv"))
logf("\n=== GSE25066 distant relapse-free survival, HR per SD ===")
for (i in seq_len(nrow(drfs))) { r <- drfs[i]
  logf(sprintf("  %-26s HR=%.3f (%.3f-%.3f) p=%.3g q=%.3g C=%.3f | ER-adj HR=%.3f p=%.3g",
    r$score, r$HR_per_SD, r$CI_low, r$CI_high, r$p, r$q_BH, r$C_index, r$HR_ER_adjusted, r$p_ER_adjusted)) }

## ------------------------------------------- (b) duplicate patients across cohorts ------
logf("\n=== duplicate-sample screen across the GPL96 neoadjuvant series ===")
gpl96 <- c("GSE25066","GSE20194","GSE22093","GSE23988","GSE42822")
common <- Reduce(intersect, lapply(store[gpl96], function(z) z$genes))
logf("genes common to all five GPL96 series: ", length(common))
gmf <- file.path(BASE, "cache/v5/farmer_gpl96_genematrix.rds")
if (file.exists(gmf)) { store$.genemat <- readRDS(gmf) } else {
  source(file.path(BASE, "scripts/10_survival_and_clinical/reactive_stroma_and_neoadjuvant/46_helpers.R"))
  mats <- lapply(gpl96, function(g) {
    s <- read_series(g); X <- s$X
    if (max(X, na.rm = TRUE) > 60) { X[X < 1] <- 1; X <- log2(X) }
    X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
    Y <- collapse_to_symbol(X, annot_map("GPL96"))[common, , drop = FALSE]
    colnames(Y) <- paste0(g, "|", colnames(Y)); Y })
  store$.genemat <- do.call(cbind, mats)
  saveRDS(store$.genemat, gmf)
}
logf("gene-level matrix for duplicate screen: ", paste(dim(store$.genemat), collapse = " x "))
## gene-level fingerprint: Spearman correlation of the full expression profile restricted
## to the 13,101 genes common to all five GPL96 series. Technical replicates of the same
## FNA hybridise to rho > 0.99; unrelated tumours on the same platform sit near 0.90-0.95.
gm <- store$.genemat
CC <- cor(gm, method = "spearman")
lab <- sub("\\|.*$", "", colnames(gm))
THR <- 0.99
dup <- which(CC > THR & upper.tri(CC), arr.ind = TRUE)
dup <- dup[lab[dup[,1]] != lab[dup[,2]], , drop = FALSE]
dt <- data.table(a = colnames(gm)[dup[,1]], b = colnames(gm)[dup[,2]], rho = CC[dup])
fwrite(dt, file.path(BASE, "results/v5/farmer_neoadjuvant_possible_duplicates.csv"))
offdiag <- CC[upper.tri(CC)][lab[which(upper.tri(CC), arr.ind = TRUE)[,1]] !=
                             lab[which(upper.tri(CC), arr.ind = TRUE)[,2]]]
logf("cross-series gene-level Spearman rho: median ", round(median(offdiag), 3),
     ", 99th pct ", round(quantile(offdiag, 0.99), 3))
logf("candidate cross-series duplicate pairs (gene-level Spearman rho > ", THR, "): ", nrow(dt))
logf("distinct samples involved: ", length(unique(c(dt$a, dt$b))))
tab <- dt[, .N, by = .(pair = paste(sub("\\|.*","",a), sub("\\|.*","",b), sep = " vs "))]
for (i in seq_len(nrow(tab))) logf("   ", tab$pair[i], ": ", tab$N[i], " candidate pairs")
logf("NOTE: these series are drawn from overlapping MDACC / USO / LAB FNA collections; the ",
     "meta-analysis is therefore reported as a summary of consistency, not as ",
     "an analysis of independent patients. GSE25066 alone (n=488 evaluable) is the ",
     "single largest and is reported separately.")

## ------------------------------------------------ (c) agreement between signatures ------
logf("\n=== Task F: agreement between the Farmer signature and the other published signatures ===")
agr <- list(); k2 <- 0L
add_cohort <- function(name, M) {
  M <- M[complete.cases(M), , drop = FALSE]
  f <- "FARMER_STROMAL"
  if (!f %in% rownames(M)) return(NULL)
  for (nm in setdiff(rownames(M), f)) {
    h <- cor.test(as.numeric(M[f, ]), as.numeric(M[nm, ]), method = "spearman", exact = FALSE)
    k2 <<- k2 + 1L
    agr[[k2]] <<- data.table(cohort = name, n = ncol(M), signature = nm,
                             rho_vs_Farmer = unname(h$estimate), p = h$p.value)
  }
}
tsc <- readRDS(file.path(BASE, "cache/v5/farmer_tcga_scores.rds"))
add_cohort("TCGA-BRCA", tsc$ssgsea)
mbf <- file.path(BASE, "cache/v5/farmer_metabric_scores.rds")
if (file.exists(mbf)) add_cohort("METABRIC", readRDS(mbf))
for (g in setdiff(names(store), ".genemat")) add_cohort(g, store[[g]]$SC)
A <- rbindlist(agr)
A[, q_BH := p.adjust(p, "BH"), by = cohort]
fwrite(A, file.path(BASE, "results/v5/farmer_signature_agreement.csv"))
W <- dcast(A, signature ~ cohort, value.var = "rho_vs_Farmer")
num <- setdiff(names(W), "signature")
W[, median_rho := apply(.SD, 1, median, na.rm = TRUE), .SDcols = num]
W[, min_rho := apply(.SD, 1, min, na.rm = TRUE), .SDcols = num]
W[, max_rho := apply(.SD, 1, max, na.rm = TRUE), .SDcols = num]
setorder(W, -median_rho)
fwrite(W, file.path(BASE, "results/v5/farmer_signature_agreement_wide.csv"))
logf("\nSpearman rho with the Farmer stromal score, per cohort (ssGSEA scores):")
print(W)
for (i in seq_len(nrow(W)))
  logf(sprintf("  %-26s median rho = %+.3f  (range %+.3f to %+.3f across %d cohorts)",
       W$signature[i], W$median_rho[i], W$min_rho[i], W$max_rho[i], length(num)))
logf("DONE 47")
