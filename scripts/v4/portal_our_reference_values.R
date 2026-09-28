# Our own TCGA-BRCA reference values for the 16 portal-cross-check genes, computed here with the
# SAME conventions as results/survival_cox_hubs.csv (Cox on z-scored log2 expression, HR per SD)
# and results/BRCA_DEX_genes.csv (limma moderated t, Tumor vs Normal).
suppressPackageStartupMessages({library(survival)})
DAT <- "/path/to/revision/data/"
RES <- "/path/to/revision/results/"
V4  <- paste0(RES, "v4/")
GENES <- c("PTEN","SMAD4","TGFBR2","DICER1","KLF4","XIAP","COL1A1","COL3A1",
           "NFKB1","RELA","SP1","ETS1","VEGFA","CCND2","MYC","EZH2")

gex <- readRDS(paste0(DAT, "brca_gene_expr.rds"))
stopifnot(!any(duplicated(rownames(gex))))            # rule 4 guard
cat("missing from expression matrix:",
    paste(setdiff(GENES, rownames(gex)), collapse=", "), "\n")

## ---- differential expression: reuse the manuscript's own limma table ----
dex <- read.csv(paste0(RES, "BRCA_DEX_genes.csv"), stringsAsFactors = FALSE)
stopifnot(is.character(dex$feature), !all(grepl("^[0-9]+$", dex$feature)))   # rule 4 assert
de <- dex[match(GENES, dex$feature), ]
de$gene <- GENES

## ---- median TPM-ish level by sample type, for a like-for-like read against UALCAN ----
ph <- readRDS(paste0(DAT, "brca_pheno.rds"))
tum <- ph$sample[ph$sample_type == "Primary Tumor"]
nor <- ph$sample[ph$sample_type == "Solid Tissue Normal"]
tum <- intersect(tum, colnames(gex)); nor <- intersect(nor, colnames(gex))
cat("tumours:", length(tum), " normals:", length(nor), "\n")

## ---- Cox, HR per SD, OS / PFI / DSS ----
sv <- readRDS(paste0(DAT, "brca_survival.rds"))
rownames(sv) <- sv$sample
cox_rows <- list()
for (ep in c("OS","PFI","DSS")) {
  ev <- sv[[ep]]; tt <- sv[[paste0(ep, ".time")]]
  keep <- !is.na(ev) & !is.na(tt) & tt > 0
  smp <- intersect(sv$sample[keep], tum)
  for (g in GENES) {
    if (!(g %in% rownames(gex))) next
    x <- as.numeric(gex[g, smp]); z <- as.numeric(scale(x))
    ok <- is.finite(z)
    fit <- coxph(Surv(sv[smp[ok], paste0(ep, ".time")], sv[smp[ok], ep]) ~ z[ok])
    s <- summary(fit)
    cox_rows[[length(cox_rows)+1]] <- data.frame(
      gene = g, endpoint = ep, n = sum(ok), n_event = sum(sv[smp[ok], ep]),
      HR_per_SD = s$coefficients[1, "exp(coef)"],
      CI_low = s$conf.int[1, "lower .95"], CI_high = s$conf.int[1, "upper .95"],
      p_value = s$coefficients[1, "Pr(>|z|)"], C_index = unname(s$concordance[1]))
  }
}
cox <- do.call(rbind, cox_rows)
cox$FDR_within_endpoint <- ave(cox$p_value, cox$endpoint, FUN = function(p) p.adjust(p, "BH"))
write.csv(cox, paste0(V4, "portals_our_tcga_cox_16genes.csv"), row.names = FALSE)

out <- data.frame(
  gene = GENES,
  our_logFC_TvsN = de$logFC, our_t_TvsN = de$t, our_p_TvsN = de$P.Value,
  our_FDR_TvsN = de$adj.P.Val,
  our_direction = ifelse(is.na(de$logFC), NA,
                         ifelse(de$logFC > 0, "UP in tumour", "DOWN in tumour")),
  our_median_tumour = sapply(GENES, function(g) if (g %in% rownames(gex))
      round(median(2^as.numeric(gex[g, tum]) - 1, na.rm = TRUE), 3) else NA),
  our_median_normal = sapply(GENES, function(g) if (g %in% rownames(gex))
      round(median(2^as.numeric(gex[g, nor]) - 1, na.rm = TRUE), 3) else NA),
  n_tumour = length(tum), n_normal = length(nor))
for (ep in c("OS","PFI","DSS")) {
  m <- cox[cox$endpoint == ep, ]
  out[[paste0("our_HRperSD_", ep)]]  <- round(m$HR_per_SD[match(GENES, m$gene)], 4)
  out[[paste0("our_CIlow_", ep)]]    <- round(m$CI_low[match(GENES, m$gene)], 4)
  out[[paste0("our_CIhigh_", ep)]]   <- round(m$CI_high[match(GENES, m$gene)], 4)
  out[[paste0("our_p_", ep)]]        <- signif(m$p_value[match(GENES, m$gene)], 4)
}
write.csv(out, paste0(V4, "portals_our_tcga_reference_16genes.csv"), row.names = FALSE)
print(out[, c("gene","our_logFC_TvsN","our_FDR_TvsN","our_direction",
              "our_median_tumour","our_median_normal","our_HRperSD_OS","our_p_OS",
              "our_HRperSD_PFI","our_p_PFI")], row.names = FALSE)

## ---- miR-130a reference (already established, re-derived here for the table) ----
mir <- readRDS(paste0(DAT, "brca_mirna_expr.rds"))
dm  <- read.csv(paste0(RES, "BRCA_DEX_mirnas.csv"), stringsAsFactors = FALSE)
mrows <- list()
for (m in c("hsa-miR-130a-3p","hsa-miR-130a-5p","hsa-miR-29a-3p","hsa-miR-29b-3p",
            "hsa-miR-29c-3p","hsa-let-7b-5p","hsa-let-7e-5p")) {
  d <- dm[dm$feature == m, ]
  for (ep in c("OS","PFI","DSS")) {
    ev <- sv[[ep]]; tt <- sv[[paste0(ep, ".time")]]
    keep <- !is.na(ev) & !is.na(tt) & tt > 0
    smp <- intersect(sv$sample[keep], intersect(colnames(mir), tum))
    if (!(m %in% rownames(mir))) next
    z <- as.numeric(scale(as.numeric(mir[m, smp]))); ok <- is.finite(z)
    fit <- coxph(Surv(sv[smp[ok], paste0(ep,".time")], sv[smp[ok], ep]) ~ z[ok])
    s <- summary(fit)
    mrows[[length(mrows)+1]] <- data.frame(
      feature = m, endpoint = ep, n = sum(ok), n_event = sum(sv[smp[ok], ep]),
      logFC_TvsN = if (nrow(d)) d$logFC[1] else NA,
      FDR_TvsN  = if (nrow(d)) d$adj.P.Val[1] else NA,
      HR_per_SD = s$coefficients[1,"exp(coef)"], CI_low = s$conf.int[1,"lower .95"],
      CI_high = s$conf.int[1,"upper .95"], p_value = s$coefficients[1,"Pr(>|z|)"])
  }
}
mres <- do.call(rbind, mrows)
write.csv(mres, paste0(V4, "portals_our_tcga_mirna_reference.csv"), row.names = FALSE)
cat("\n=== our TCGA miRNA reference ===\n"); print(mres, row.names = FALSE)
