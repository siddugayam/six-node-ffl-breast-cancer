#!/usr/bin/env Rscript
# v4/05 -- (B) ssGSEA / GSVA per-sample scores for the top pathways (survival-arm input)
suppressPackageStartupMessages({
  library(GSVA); library(data.table); library(BiocParallel); library(survival)
})
ROOT  <- "/path/to/revision"
RES   <- file.path(ROOT, "results", "v4")
CACHE <- file.path(ROOT, "cache", "v4")
set.seed(20260908)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

sets <- readRDS(file.path(CACHE, "genesets.rds"))
A    <- fread(file.path(RES, "gsea_all_results.csv"))

## ---- choose "top pathways": per ranked list x interpretable collection ------
PRIMARY <- c("H", "C2_CP_KEGG", "C2_CP_REACTOME", "C6", "C3_MIR_LEGACY", "C3_TFT_LEGACY")
TOPN    <- 25L
sel <- A[padj < 0.05 & collection %in% PRIMARY]
sel[, absNES := abs(NES)]
setorder(sel, ranked_list, collection, -absNES)
top <- sel[, head(.SD, TOPN), by = .(ranked_list, collection)]
chosen <- unique(top[, .(collection, pathway)])
## always include every FFL class set + the two miR-29 positive-control sets
chosen <- rbind(chosen,
                data.table(collection = "FFL_CLASS", pathway = names(sets$FFL_CLASS)),
                data.table(collection = "C3_MIR_LEGACY", pathway = "TGGTGCT_MIR29A_MIR29B_MIR29C"),
                data.table(collection = "C3_MIR_MIRDB",  pathway = "MIR29A_3P"))
chosen <- unique(chosen)
gs <- list()
for (i in seq_len(nrow(chosen))) {
  cl <- chosen$collection[i]; pw <- chosen$pathway[i]
  if (!is.null(sets[[cl]][[pw]])) gs[[paste0(cl, "|", pw)]] <- sets[[cl]][[pw]]
}
msg("selected ", length(gs), " gene sets for per-sample scoring")
fwrite(data.table(set_id = names(gs), collection = sub("\\|.*", "", names(gs)),
                  pathway = sub("^[^|]*\\|", "", names(gs)), size = lengths(gs)),
       file.path(RES, "gsea_ssgsea_setlist.csv"))

## ---- expression ------------------------------------------------------------
expr <- readRDS(file.path(ROOT, "data", "brca_gene_expr.rds"))
stopifnot(!any(duplicated(rownames(expr))))
expr <- expr[matrixStats::rowVars(expr) > 0, , drop = FALSE]
msg("expression for scoring: ", nrow(expr), " x ", ncol(expr))

## ---- ssGSEA ----------------------------------------------------------------
t0 <- Sys.time()
p_ss <- ssgseaParam(exprData = expr, geneSets = gs, minSize = 5, maxSize = 5000,
                    alpha = 0.25, normalize = TRUE)
ss <- gsva(p_ss, verbose = FALSE, BPPARAM = MulticoreParam(workers = 8))
msg("ssGSEA done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min; ",
    nrow(ss), " sets x ", ncol(ss), " samples")
saveRDS(ss, file.path(RES, "gsea_ssgsea_scores_matrix.rds"))
fwrite(data.table(set_id = rownames(ss), as.data.table(ss)),
       file.path(RES, "gsea_ssgsea_scores_matrix.csv"))

## ---- GSVA (Gaussian kernel) as a second estimator ---------------------------
t0 <- Sys.time()
p_gv <- gsvaParam(exprData = expr, geneSets = gs, minSize = 5, maxSize = 5000, kcdf = "Gaussian")
gv <- gsva(p_gv, verbose = FALSE, BPPARAM = MulticoreParam(workers = 8))
msg("GSVA done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")
saveRDS(gv, file.path(RES, "gsea_gsva_scores_matrix.rds"))
fwrite(data.table(set_id = rownames(gv), as.data.table(gv)),
       file.path(RES, "gsea_gsva_scores_matrix.csv"))
cm <- sapply(intersect(rownames(ss), rownames(gv)),
             function(r) cor(ss[r, ], gv[r, colnames(ss)], method = "spearman"))
msg("ssGSEA vs GSVA per-set Spearman: median ", signif(median(cm), 4),
    " (min ", signif(min(cm), 4), ", max ", signif(max(cm), 4), ")")

## ---- sample annotation -----------------------------------------------------
ph  <- readRDS(file.path(ROOT, "data", "brca_pheno.rds"))
srv <- readRDS(file.path(ROOT, "data", "brca_survival.rds"))
mir <- readRDS(file.path(ROOT, "data", "brca_mirna_expr.rds"))
ann <- data.table(sample = colnames(ss))
ann[, sample_type := ph$sample_type[match(sample, ph$sample)]]
ann[, mir29a_3p := ifelse(sample %in% colnames(mir), mir["hsa-miR-29a-3p", match(sample, colnames(mir))], NA_real_)]
for (v in c("OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time",
            "ajcc_pathologic_tumor_stage", "age_at_initial_pathologic_diagnosis"))
  ann[[v]] <- srv[[v]][match(ann$sample, srv$sample)]
fwrite(ann, file.path(RES, "gsea_ssgsea_sample_annotation.csv"))
msg("annotation written: ", nrow(ann), " samples; tumours=",
    sum(ann$sample_type == "Primary Tumor", na.rm = TRUE))

## ---- quick per-set readouts for the survival arm ---------------------------
tum <- ann[sample_type == "Primary Tumor" & !is.na(OS.time) & OS.time > 0]
res <- rbindlist(lapply(rownames(ss), function(r) {
  s <- ss[r, tum$sample]
  z <- as.numeric(scale(s))
  cx <- tryCatch(summary(coxph(Surv(tum$OS.time, tum$OS) ~ z)), error = function(e) NULL)
  pf <- tum[!is.na(PFI.time) & PFI.time > 0]
  zp <- as.numeric(scale(ss[r, pf$sample]))
  cp <- tryCatch(summary(coxph(Surv(pf$PFI.time, pf$PFI) ~ zp)), error = function(e) NULL)
  tn <- ann[!is.na(sample_type)]
  wtn <- tryCatch(wilcox.test(ss[r, tn$sample] ~ tn$sample_type), error = function(e) NULL)
  data.table(set_id = r,
             OS_HR_perSD = if (is.null(cx)) NA_real_ else cx$coefficients[1, 2],
             OS_p = if (is.null(cx)) NA_real_ else cx$coefficients[1, 5],
             OS_n = nrow(tum), OS_events = sum(tum$OS, na.rm = TRUE),
             PFI_HR_perSD = if (is.null(cp)) NA_real_ else cp$coefficients[1, 2],
             PFI_p = if (is.null(cp)) NA_real_ else cp$coefficients[1, 5],
             PFI_n = nrow(pf), PFI_events = sum(pf$PFI, na.rm = TRUE),
             tumour_minus_normal = mean(ss[r, tn$sample[tn$sample_type == "Primary Tumor"]]) -
                                   mean(ss[r, tn$sample[tn$sample_type == "Solid Tissue Normal"]]),
             tumour_vs_normal_wilcox_p = if (is.null(wtn)) NA_real_ else wtn$p.value)
}))
res[, `:=`(OS_fdr = p.adjust(OS_p, "BH"), PFI_fdr = p.adjust(PFI_p, "BH"),
           tumour_vs_normal_fdr = p.adjust(tumour_vs_normal_wilcox_p, "BH"))]
setorder(res, OS_p)
fwrite(res, file.path(RES, "gsea_ssgsea_set_associations.csv"))
msg("set-level associations: OS FDR<0.05 -> ", sum(res$OS_fdr < 0.05, na.rm = TRUE),
    " ; PFI FDR<0.05 -> ", sum(res$PFI_fdr < 0.05, na.rm = TRUE))
print(head(res[, .(set_id, OS_HR_perSD, OS_p, OS_fdr, PFI_HR_perSD, PFI_p)], 15))
msg("DONE 05")
