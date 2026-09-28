#!/usr/bin/env Rscript
# 04a_extract_mir29_data.R -- Part D data extraction.
# Pull the measured quantities that constrain the miR-29 / collagen circuit:
#   TCGA-BRCA  : log2 expression of miR-29a/b/c, NFKB1, SP1, RELA, COL1A1, COL3A1
#   CPTAC-BRCA : matched mRNA and protein for the same genes, plus miR-29 arms
# Nothing is modelled here; this script only writes the measurement tables.
suppressMessages(library(stats))
REV <- "/path/to/revision"
RES <- file.path(REV, "results", "v3")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf("%s | %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

G  <- readRDS(file.path(REV, "data", "brca_gene_expr.rds"))
Mi <- readRDS(file.path(REV, "data", "brca_mirna_expr_canonical.rds"))
ph <- readRDS(file.path(REV, "data", "brca_pheno.rds"))
stopifnot(!any(duplicated(rownames(G))), !any(duplicated(rownames(Mi))))
say("gene matrix ", paste(dim(G), collapse = " x "), " ; miRNA ", paste(dim(Mi), collapse = " x "))

tum <- ph$sample[ph$sample_type == "Primary Tumor"]
s   <- intersect(intersect(colnames(G), colnames(Mi)), tum)
say("paired primary tumours: ", length(s))

genes  <- c("COL1A1", "COL3A1", "COL1A2", "NFKB1", "RELA", "SP1", "ETS1")
mirs   <- c("hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c")
stopifnot(all(genes %in% rownames(G)), all(mirs %in% rownames(Mi)))

D <- data.frame(sample = s, stringsAsFactors = FALSE)
for (g in genes) D[[g]] <- as.numeric(G[g, s])
for (m in mirs)  D[[gsub("-", "_", m)]] <- as.numeric(Mi[m, s])
# miR-29 family aggregate: log2 of the summed linear abundance of the three members
D$miR29_family <- log2(rowSums(2 ^ D[, gsub("-", "_", mirs)]))
write.csv(D, file.path(RES, "dynamics_mir29_tcga_measurements.csv"), row.names = FALSE)
say("wrote TCGA table: ", nrow(D), " tumours")

## ---- Spearman correlations actually used to constrain the circuit -----------
cr <- do.call(rbind, lapply(c(gsub("-", "_", mirs), "miR29_family"), function(m)
  do.call(rbind, lapply(c("COL1A1", "COL3A1"), function(g) {
    ct <- suppressWarnings(cor.test(D[[m]], D[[g]], method = "spearman"))
    lm1 <- lm(D[[g]] ~ D[[m]])
    data.frame(cohort = "TCGA-BRCA", layer = "mRNA", miRNA = m, target = g, n = nrow(D),
               rho = unname(ct$estimate), p = ct$p.value,
               ols_slope_log2_per_log2 = unname(coef(lm1)[2]),
               ols_slope_se = summary(lm1)$coefficients[2, 2],
               r2 = summary(lm1)$r.squared)
  }))))
say("TCGA correlations computed")

## ---- CPTAC: mRNA and protein ------------------------------------------------
b <- readRDS(file.path(REV, "results", "multiomics", "cptac_bundle.rds"))
sc <- b$s_all
Rg <- b$Rg; Pg <- b$Pg; MI <- b$MI; MIc <- b$MIc
arms <- c("hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p")
say("CPTAC samples with RNA+protein+miRNA: ", length(sc))

C <- data.frame(sample = sc, stringsAsFactors = FALSE)
for (g in c("COL1A1", "COL3A1", "NFKB1", "SP1", "RELA")) {
  C[[paste0(g, "_mRNA")]]    <- if (g %in% rownames(Rg)) as.numeric(Rg[g, sc]) else NA_real_
  C[[paste0(g, "_protein")]] <- if (g %in% rownames(Pg)) as.numeric(Pg[g, sc]) else NA_real_
}
for (a in arms) C[[gsub("-", "_", a)]] <- if (a %in% rownames(MI)) as.numeric(MI[a, sc]) else NA_real_
C$miR29_family_3p <- log2(rowSums(2 ^ C[, gsub("-", "_", arms)], na.rm = TRUE))
write.csv(C, file.path(RES, "dynamics_mir29_cptac_measurements.csv"), row.names = FALSE)

cc <- do.call(rbind, lapply(c(gsub("-", "_", arms), "miR29_family_3p"), function(m)
  do.call(rbind, lapply(c("COL1A1", "COL3A1"), function(g) {
    do.call(rbind, lapply(c("mRNA", "protein"), function(ly) {
      yv <- C[[paste0(g, "_", ly)]]; xv <- C[[m]]
      ok <- is.finite(yv) & is.finite(xv)
      if (sum(ok) < 20) return(NULL)
      ct <- suppressWarnings(cor.test(xv[ok], yv[ok], method = "spearman"))
      lm1 <- lm(yv[ok] ~ xv[ok])
      data.frame(cohort = "CPTAC-BRCA", layer = ly, miRNA = m, target = g, n = sum(ok),
                 rho = unname(ct$estimate), p = ct$p.value,
                 ols_slope_log2_per_log2 = unname(coef(lm1)[2]),
                 ols_slope_se = summary(lm1)$coefficients[2, 2],
                 r2 = summary(lm1)$r.squared)
    }))
  }))))

out <- rbind(cr, cc)
write.csv(out, file.path(RES, "dynamics_mir29_effect_sizes.csv"), row.names = FALSE)
print(out, row.names = FALSE)
say("wrote effect sizes: ", nrow(out), " rows")
