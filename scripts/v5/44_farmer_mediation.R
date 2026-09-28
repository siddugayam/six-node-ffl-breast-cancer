#!/usr/bin/env Rscript
## 44_farmer_mediation.R
## Task D: repeat the causal-mediation analysis of the TF -> collagen and
## miR-29 -| collagen edges using the *independently published, clinically validated*
## Farmer stromal signature as the mediator, alongside our own CAF estimates.
## Circularity guards: the Farmer signature contains 6 collagens and overlaps the
## network, so the primary mediator is FARMER_clean = Farmer 50 minus all COL* genes
## minus every gene that is a node of our network.
suppressPackageStartupMessages({library(mediation); library(matrixStats); library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v5/44_farmer_mediation.log")
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG)
SIMS <- 1000L
logf("mediation bootstrap sims = ", SIMS, " (boot=TRUE, percentile CI)")

sets  <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
gexp  <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mexp  <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
ph    <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
nodes <- read.delim(file.path(BASE, "data/canonical_nodes.tsv"), stringsAsFactors = FALSE)$name
tum   <- ph$sample[ph$sample_type == "Primary Tumor"]
sG <- sort(intersect(colnames(gexp), tum)); logf("gene-assay primary tumours n = ", length(sG))
sB <- sort(intersect(sG, colnames(mexp))); logf("both-assay primary tumours  n = ", length(sB))

## ------------------------------------------------------------------ mediator sets ------
farmer      <- intersect(sets$FARMER_STROMAL, rownames(gexp))
farmer_noC  <- setdiff(farmer, grep("^COL", farmer, value = TRUE))
farmer_clean<- setdiff(farmer_noC, nodes)
gmt  <- readLines(system.file("extdata", "SI_geneset.gmt", package = "estimate"))
strom <- strsplit(gmt[grep("StromalSignature", gmt)], "\t")[[1]][-c(1, 2)]
sigA <- intersect(setdiff(setdiff(strom, grep("^COL", strom, value = TRUE)), nodes), rownames(gexp))
sigB <- intersect(setdiff(c("DCN","LUM","FAP","PDGFRB","THY1","POSTN"), nodes), rownames(gexp))
cafsc<- intersect(setdiff(sets$CAF_scRNA_50, nodes), rownames(gexp))
MED <- list(FARMER_full = farmer, FARMER_noCOL = farmer_noC, FARMER_clean = farmer_clean,
            CAF_A_estimate_noCOL_noNodes = sigA, CAF_B_markers = sigB, CAF_scRNA = cafsc)
for (nm in names(MED))
  logf(sprintf("mediator %-30s n=%3d  contains COL*: %s  contains network nodes: %s", nm,
       length(MED[[nm]]), any(grepl("^COL", MED[[nm]])), any(MED[[nm]] %in% nodes)))
logf("FARMER_clean members: ", paste(sort(farmer_clean), collapse = ", "))

zscore <- function(S, ss) { X <- gexp[S, ss, drop = FALSE]; colMeans((X - rowMeans(X)) / rowSds(X)) }
ns <- function(v) { n <- length(v); qnorm((rank(v) - 0.5) / n) }

run_med <- function(xname, yname, medname, ss, xassay = "gene") {
  X <- if (xassay == "miRNA") mexp[xname, ss] else gexp[xname, ss]
  Y <- gexp[yname, ss]
  gs <- setdiff(MED[[medname]], yname)          # never let the outcome be in the mediator
  M  <- zscore(gs, ss)
  d  <- data.frame(x = ns(X), y = ns(Y), m = ns(M))
  mm <- lm(m ~ x, data = d); om <- lm(y ~ x + m, data = d)
  set.seed(7)
  fit <- mediate(mm, om, treat = "x", mediator = "m", boot = TRUE, sims = SIMS,
                 boot.ci.type = "perc")
  a <- coef(mm)["x"]; b <- coef(om)["m"]
  sa <- summary(mm)$coef["x", 2]; sb <- summary(om)$coef["m", 2]
  z  <- a * b / sqrt(b^2 * sa^2 + a^2 * sb^2)
  data.frame(x = xname, y = yname, mediator = medname, n_mediator_genes = length(gs),
    n = length(ss), rho_raw = cor(X, Y, method = "spearman"),
    rho_partial_given_mediator = suppressWarnings(
      cor(residuals(lm(ns(X) ~ ns(M))), residuals(lm(ns(Y) ~ ns(M))), method = "spearman")),
    total = fit$tau.coef, total_p = fit$tau.p,
    ACME = fit$d0, acme_lo = fit$d0.ci[1], acme_hi = fit$d0.ci[2], acme_p = fit$d0.p,
    ADE = fit$z0, ade_lo = fit$z0.ci[1], ade_hi = fit$z0.ci[2], ade_p = fit$z0.p,
    prop_mediated = fit$n0, prop_lo = fit$n0.ci[1], prop_hi = fit$n0.ci[2], prop_p = fit$n0.p,
    a_path = a, b_path = b, cprime = coef(om)["x"], sobel_z = z, sobel_p = 2 * pnorm(-abs(z)),
    stringsAsFactors = FALSE)
}

gene_jobs  <- list(c("ETS1","COL1A1"), c("NFKB1","COL1A1"), c("SP1","COL1A1"),
                   c("RELA","COL1A1"), c("ETS1","COL3A1"), c("NFKB1","COL3A1"))
mirna_jobs <- list(c("hsa-miR-29a","COL1A1"), c("hsa-miR-29a","COL3A1"), c("hsa-miR-29b","COL3A1"))
res <- list(); k <- 0L
for (md in names(MED)) {
  for (j in gene_jobs)  { k <- k + 1L; res[[k]] <- run_med(j[1], j[2], md, sG) }
  for (j in mirna_jobs) { k <- k + 1L; res[[k]] <- run_med(j[1], j[2], md, sB, xassay = "miRNA") }
  logf("finished mediator ", md)
}
R <- as.data.table(rbindlist(res))
fwrite(R, file.path(BASE, "results/v5/farmer_mediation.csv"))
logf("WROTE results/v5/farmer_mediation.csv rows=", nrow(R))
logf("\n=== mediation: does an independently published stromal signature give the same answer? ===")
logf(sprintf("%-12s %-8s %-30s %6s %8s %8s %8s %8s %10s", "exposure","outcome","mediator",
             "rho","rho|M","total","ACME","ADE","prop_med"))
for (i in seq_len(nrow(R))) { r <- R[i]
  logf(sprintf("%-12s %-8s %-30s %+6.3f %+8.3f %+8.3f %+8.3f %+8.3f %9.3f  (ACME p=%.3g, prop 95%% CI %.2f-%.2f)",
    r$x, r$y, r$mediator, r$rho_raw, r$rho_partial_given_mediator, r$total, r$ACME, r$ADE,
    r$prop_mediated, r$acme_p, r$prop_lo, r$prop_hi)) }

## side-by-side agreement between the published Farmer mediator and our own CAF estimates
w <- dcast(R, x + y ~ mediator, value.var = "prop_mediated")
fwrite(w, file.path(BASE, "results/v5/farmer_mediation_propmediated_wide.csv"))
logf("\n=== proportion mediated, side by side ===")
print(w)
cc <- R[mediator %in% c("FARMER_clean","CAF_A_estimate_noCOL_noNodes")]
a <- cc[mediator == "FARMER_clean"][order(x, y)]; b <- cc[mediator == "CAF_A_estimate_noCOL_noNodes"][order(x, y)]
logf("Pearson r between the two mediators' proportion-mediated estimates across the 9 axes: ",
     round(cor(a$prop_mediated, b$prop_mediated), 3),
     " | Spearman ", round(cor(a$prop_mediated, b$prop_mediated, method = "spearman"), 3))
logf("DONE 44")
