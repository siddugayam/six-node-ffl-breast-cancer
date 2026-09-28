# TIMER2.0 cross-check: correlate our TCGA-BRCA miRNA / gene values against TIMER2.0's own
# published deconvolution estimates (infiltration_estimation_for_tcga.csv, the exact file
# TIMER2.0 serves for download). This is an independent check on our CAF/fibroblast finding.
suppressPackageStartupMessages({library(data.table)})
RES <- "/path/to/revision/results/v4/"
DAT <- "/path/to/revision/data/"

ti <- fread(cmd = paste0("zcat ", DAT, "portals/timer2_infiltration_estimation_for_tcga.csv.gz"))
setnames(ti, 1, "sample")
cat("TIMER2 rows:", nrow(ti), "cols:", ncol(ti), "\n")

# TIMER barcodes are TCGA-XX-XXXX-01 ; ours are the same form
mir <- readRDS(paste0(DAT, "brca_mirna_expr.rds"))
gex <- readRDS(paste0(DAT, "brca_gene_expr.rds"))
stopifnot(!any(duplicated(rownames(gex))))

ti$sample <- toupper(ti$sample)
common_mir <- intersect(ti$sample, colnames(mir))
common_gex <- intersect(ti$sample, colnames(gex))
cat("overlap with our miRNA matrix:", length(common_mir),
    "| with our gene matrix:", length(common_gex), "\n")

infil_cols <- setdiff(names(ti), "sample")
tim <- as.data.frame(ti); rownames(tim) <- tim$sample

# focus columns: every fibroblast/endothelial/purity-ish estimate + the main immune ones
focus <- grep("fibroblast|Endothelial|uncharacterized|T cell CD8|Macrophage_TIMER|B cell_TIMER|Neutrophil_TIMER|T cell CD4\\+_TIMER|Myeloid dendritic cell_TIMER",
              infil_cols, value = TRUE, ignore.case = TRUE)
cat("focus cell types:", length(focus), "\n")

feats <- list(
  `hsa-miR-130a-3p` = mir["hsa-miR-130a-3p", ],
  `hsa-miR-130a-5p` = mir["hsa-miR-130a-5p", ],
  `hsa-miR-29a-3p`  = mir["hsa-miR-29a-3p", ]
)
for (g in c("COL1A1","COL3A1","PTEN","SMAD4","TGFBR2","DICER1","KLF4","VEGFA","EZH2","MYC",
            "DCN","LUM","FAP","POSTN","THY1","PDGFRB","EPCAM","PTPRC","ESR1","VIM"))
  if (g %in% rownames(gex)) feats[[g]] <- gex[g, ]

out <- list()
for (fn in names(feats)) {
  v <- feats[[fn]]
  smp <- intersect(names(v), rownames(tim))
  smp <- smp[grepl("-01$", smp)]                       # primary tumours only
  for (cc in focus) {
    a <- as.numeric(v[smp]); b <- as.numeric(tim[smp, cc])
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 30 || length(unique(b[ok])) < 5) next
    ct <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman"))
    out[[length(out) + 1]] <- data.frame(
      feature = fn, timer_cell_type = cc, method = sub(".*_", "", cc),
      n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
  }
}
res <- do.call(rbind, out)
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(res$feature, res$p), ]
write.csv(res, paste0(RES, "portals_timer2_infiltration_correlations.csv"), row.names = FALSE)
cat("\n=== miR-130a-3p vs fibroblast estimates ===\n")
print(res[res$feature == "hsa-miR-130a-3p" &
          grepl("fibroblast", res$timer_cell_type, ignore.case = TRUE), ], row.names = FALSE)
cat("\n=== miR-29a-3p vs fibroblast estimates ===\n")
print(res[res$feature == "hsa-miR-29a-3p" &
          grepl("fibroblast", res$timer_cell_type, ignore.case = TRUE), ], row.names = FALSE)
cat("\n=== COL1A1 / COL3A1 vs fibroblast estimates ===\n")
print(res[res$feature %in% c("COL1A1","COL3A1") &
          grepl("fibroblast", res$timer_cell_type, ignore.case = TRUE), ], row.names = FALSE)
cat("\nrows written:", nrow(res), "\n")
