## E6: why Table S1 (published as Table S2) gives rho = -0.246 for miR-29a -| COL1A1 (n = 1,066) and
## results/v6/pancancer_axes.csv gives -0.160 (BRCA, n = 1,065).  Both are Spearman correlations on
## PanCanAtlas matrices; this script recomputes each and swaps the miRNA feature and the sample set one at
## a time.  Only summary numbers are written.  Read-only on the project.
suppressPackageStartupMessages(library(data.table))
REV <- "/path/to/revision"
G0 <- readRDS(file.path(REV, "data/brca_gene_expr.rds")); M0 <- readRDS(file.path(REV, "data/brca_mirna_expr.rds"))
MC <- readRDS(file.path(REV, "data/brca_mirna_expr_canonical.rds")); ph <- readRDS(file.path(REV, "data/brca_pheno.rds"))
tum <- ph$sample[ph$sample_type == "Primary Tumor"]
s1066 <- sort(intersect(intersect(colnames(G0), tum), intersect(colnames(M0), tum)))
cat("08_expression_validation.R sample set (primary tumour, both assays):", length(s1066), "\n")
arm3 <- M0["hsa-miR-29a-3p", s1066]; arm5 <- M0["hsa-miR-29a-5p", s1066]; can <- MC["hsa-miR-29a", s1066]
cat("canonical hsa-miR-29a equals the mean of the 3p and 5p rows:", isTRUE(all.equal(unname(can), unname((arm3 + arm5) / 2))), "\n")
col <- G0["COL1A1", s1066]
sp <- function(a, b) { ok <- is.finite(a) & is.finite(b); c(rho = unname(cor(a[ok], b[ok], method = "spearman")), n = sum(ok)) }
r <- rbind(`canonical (mean of 3p and 5p), 1,066 samples` = sp(can, col),
           `3p arm only, same samples` = sp(arm3, col), `5p arm only, same samples` = sp(arm5, col))
PC <- fread(file.path(REV, "results/v6/pancancer_axes.csv"))[cohort == "BRCA" & axis == "miR-29a -| COL1A1"]
cat("stored pancancer_axes.csv BRCA row: rho", PC$rho, "n", PC$n, "\n")
print(r)
out <- data.table(comparison = rownames(r), rho = r[, "rho"], n = r[, "n"])
out <- rbind(out, data.table(comparison = "stored results/v6/pancancer_axes.csv (3p arm, pan-cancer matrices)", rho = PC$rho, n = PC$n))
fwrite(out, "e6_two_correlations.csv")
