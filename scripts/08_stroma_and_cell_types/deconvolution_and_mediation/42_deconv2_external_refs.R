#!/usr/bin/env Rscript
## 42_deconv2_external_refs.R
## Two EXTERNAL, independently generated estimates of tumour composition for the same
## TCGA-BRCA tumours, neither produced by me:
##  (1) Thorsson et al. 2018 (Immunity 48:812) PanImmune official CIBERSORT LM22
##      relative fractions, run by the TCGA PanCanAtlas group on Kallisto-quantified
##      reads. GDC file TCGA.Kallisto.fullIDs.cibersort.relative.tsv
##      -> an external check on my in-house nu-SVR CIBERSORT re-implementation.
##  (2) ABSOLUTE tumour purity (TCGA_mastercalls.abs_tables_JSedit.fixed.txt), a
##      DNA-based (copy-number) estimate. 1 - purity is a NON-TRANSCRIPTOMIC estimate
##      of non-malignant content and therefore cannot be circular with COL1A1 mRNA.
## Output: cache/v3/deconv2/external_estimates.rds
##         results/v3/deconv_external_reference_agreement.csv
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/42_deconv2_external.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
dir.create(file.path(BASE, "cache/v3/deconv2"), recursive = TRUE, showWarnings = FALSE)

SS <- readLines(file.path(BASE, "cache/v3/deconv/samples_primary_tumour.txt"))
logf("target samples: ", length(SS))

## ------------------------------------------------------------------ (1) Thorsson ----
TH <- fread(file.path(BASE, "data/deconv/thorsson_cibersort_relative.tsv"))
logf("Thorsson file: ", nrow(TH), " rows x ", ncol(TH), " cols; cancer types = ", uniqueN(TH$CancerType))
TH <- TH[CancerType == "BRCA"]
logf("BRCA rows: ", nrow(TH))
## SampleID like TCGA.3C.AAAU.01A.11R.A41B.07 -> TCGA-3C-AAAU-01
TH[, short := {
  p <- tstrsplit(SampleID, "\\.")
  paste(p[[1]], p[[2]], p[[3]], substr(p[[4]], 1, 2), sep = "-")
}]
logf("unique short IDs: ", uniqueN(TH$short), "; duplicated short IDs: ", sum(duplicated(TH$short)))
TH <- TH[!duplicated(short)]
setkey(TH, short)
ov <- intersect(SS, TH$short)
logf("overlap with our 1,097 primary tumours: ", length(ov))
cellcols <- setdiff(colnames(TH), c("SampleID", "CancerType", "short"))
logf("Thorsson columns: ", paste(cellcols, collapse = ", "))
THM <- as.matrix(TH[ov, ..cellcols]); rownames(THM) <- ov
mode(THM) <- "numeric"
## keep only the 22 LM22 fractions (the file also carries P.value/Correlation/RMSE)
qc <- intersect(c("P.value", "Correlation", "RMSE"), colnames(THM))
logf("QC columns present: ", paste(qc, collapse = ", "))
if (length(qc)) {
  logf("Thorsson median mixture correlation = ", round(median(THM[, "Correlation"], na.rm = TRUE), 3),
       "; samples with P.value < 0.05 = ", sum(THM[, "P.value"] < 0.05, na.rm = TRUE), "/", nrow(THM))
}
frac <- setdiff(colnames(THM), qc)
THF <- THM[, frac, drop = FALSE]
logf("Thorsson LM22 fractions: ", ncol(THF), " types; row sums median = ",
     round(median(rowSums(THF, na.rm = TRUE)), 4))

## agreement with my in-house CIBERSORT
CB <- readRDS(file.path(BASE, "cache/v3/deconv/out/CIBERSORT_LM22.rds"))
logf("in-house CIBERSORT_LM22: ", nrow(CB), " x ", ncol(CB), " -> ", paste(head(colnames(CB), 3), collapse = ", "))
norm_nm <- function(v) tolower(gsub("[^a-z0-9]", "", tolower(v)))
mapcol <- match(norm_nm(colnames(THF)), norm_nm(colnames(CB)))
logf("LM22 column name matching: ", sum(!is.na(mapcol)), "/", ncol(THF), " matched")
agr <- rbindlist(lapply(which(!is.na(mapcol)), function(k) {
  a <- THF[ov, k]; b <- CB[ov, mapcol[k]]
  data.table(cell_type = colnames(THF)[k],
             spearman = cor(a, b, method = "spearman", use = "complete.obs"),
             pearson  = cor(a, b, method = "pearson",  use = "complete.obs"),
             mean_thorsson = mean(a, na.rm = TRUE), mean_inhouse = mean(b, na.rm = TRUE))
}))
logf("\n=== in-house nu-SVR CIBERSORT vs official Thorsson CIBERSORT (n = ", length(ov), " tumours) ===")
for (i in seq_len(nrow(agr))) logf(sprintf("%-32s rho = %+0.3f  r = %+0.3f  mean %.4f vs %.4f",
  agr$cell_type[i], agr$spearman[i], agr$pearson[i], agr$mean_thorsson[i], agr$mean_inhouse[i]))
logf("median Spearman across the 22 LM22 types = ", round(median(agr$spearman, na.rm = TRUE), 3))
fwrite(agr, file.path(BASE, "results/v3/deconv_external_reference_agreement.csv"))

## ------------------------------------------------------------------ (2) ABSOLUTE ----
AB <- fread(file.path(BASE, "data/deconv/TCGA_ABSOLUTE_purity.txt"))
logf("\nABSOLUTE table: ", nrow(AB), " rows; cols = ", paste(colnames(AB), collapse = ", "))
AB[, short := substr(array, 1, 15)]
AB <- AB[!duplicated(short)]
ab <- AB[match(SS, short)]
pur <- ab$purity
names(pur) <- SS
logf("ABSOLUTE purity available for ", sum(!is.na(pur)), "/", length(SS), " tumours; ",
     "median purity = ", round(median(pur, na.rm = TRUE), 3))

EXT <- list(thorsson_cibersort = THF, absolute_purity = pur, samples = SS)
saveRDS(EXT, file.path(BASE, "cache/v3/deconv2/external_estimates.rds"))
logf("WROTE cache/v3/deconv2/external_estimates.rds")
logf("DONE 42")
