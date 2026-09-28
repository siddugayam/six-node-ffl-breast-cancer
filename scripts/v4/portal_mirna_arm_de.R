# Arm-level (5p/3p) tumour-vs-normal limma DE for the cross-check miRNAs.
# results/BRCA_DEX_mirnas.csv is keyed on precursor names ("hsa-miR-130a"), so it cannot supply
# the arm-level logFC the manuscript quotes for hsa-miR-130a-3p. Recomputed here.
suppressPackageStartupMessages(library(limma))
DAT <- "/path/to/revision/data/"
V4  <- "/path/to/revision/results/v4/"
mir <- readRDS(paste0(DAT, "brca_mirna_expr.rds"))
ph  <- readRDS(paste0(DAT, "brca_pheno.rds"))
stopifnot(!any(duplicated(rownames(mir))))
grp <- ph$sample_type[match(colnames(mir), ph$sample)]
keep <- grp %in% c("Primary Tumor", "Solid Tissue Normal")
E <- mir[, keep]; g <- factor(grp[keep], levels = c("Solid Tissue Normal", "Primary Tumor"))
cat("tumours:", sum(g == "Primary Tumor"), " normals:", sum(g == "Solid Tissue Normal"), "\n")
des <- model.matrix(~ g)
fit <- eBayes(lmFit(E, des))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
stopifnot(!all(grepl("^[0-9]+$", rownames(tt))))            # rule 4 guard
stopifnot(identical(rownames(tt), rownames(E)))
tt$feature <- rownames(tt)
sel <- c("hsa-miR-130a-3p","hsa-miR-130a-5p","hsa-miR-29a-3p","hsa-miR-29b-3p",
         "hsa-miR-29c-3p","hsa-let-7b-5p","hsa-let-7e-5p","hsa-miR-130b-3p","hsa-miR-301a-3p")
out <- tt[tt$feature %in% sel, c("feature","logFC","AveExpr","t","P.Value","adj.P.Val")]
out$direction <- ifelse(out$logFC > 0, "UP in tumour", "DOWN in tumour")
out$n_tumour <- sum(g == "Primary Tumor"); out$n_normal <- sum(g == "Solid Tissue Normal")
write.csv(out, paste0(V4, "portals_our_tcga_mirna_arm_DE.csv"), row.names = FALSE)
print(out, row.names = FALSE)
