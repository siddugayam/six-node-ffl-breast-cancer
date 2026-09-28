## 22_metabolic_tumour_vs_normal.R ----------------------------------------
## Part B: tumour vs normal for every metabolic score, with FDR.
suppressPackageStartupMessages({library(data.table); library(limma)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
ss <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
gv <- readRDS(file.path(RES,"metabolic_gsva_scores.rds"))
ph <- readRDS(file.path(ROOT,"data/brca_pheno.rds"))
inv <- fread(file.path(RES,"metabolic_geneset_inventory.csv"))

ph <- ph[match(colnames(ss), ph$sample), ]
stopifnot(identical(ph$sample, colnames(ss)), !any(is.na(ph$sample_type)))
grp <- factor(ifelse(ph$sample_type == "Primary Tumor", "Tumor", "Normal"), levels=c("Normal","Tumor"))
cat("n Tumor:", sum(grp=="Tumor"), " n Normal:", sum(grp=="Normal"), "\n")

run_tvn <- function(M, tag) {
  des <- model.matrix(~ grp)
  fit <- eBayes(lmFit(M, des))
  tt  <- topTable(fit, coef="grpTumor", number=Inf, sort.by="none")
  stopifnot(identical(rownames(tt), rownames(M)))   # guard vs limma integer-index failure
  stopifnot(!all(grepl("^[0-9]+$", rownames(tt))))
  w <- apply(M, 1, function(v) {
    x <- v[grp=="Tumor"]; y <- v[grp=="Normal"]
    wt <- suppressWarnings(wilcox.test(x, y))
    auc <- (wt$statistic / (length(x)*length(y)))
    c(p_wilcox=unname(wt$p.value), auc_tumour_higher=unname(auc),
      mean_T=mean(x), mean_N=mean(y),
      cohens_d=(mean(x)-mean(y))/sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y))/(length(x)+length(y)-2)))
  })
  data.table(set_id=rownames(M), score=tag,
             logFC_T_minus_N=tt$logFC, t=tt$t, p_limma=tt$P.Value,
             FDR_limma=p.adjust(tt$P.Value, "BH"),
             p_wilcox=w["p_wilcox",], FDR_wilcox=p.adjust(w["p_wilcox",], "BH"),
             auc_tumour_higher=w["auc_tumour_higher",],
             mean_tumour=w["mean_T",], mean_normal=w["mean_N",], cohens_d=w["cohens_d",])
}
r_ss <- run_tvn(ss, "ssGSEA"); r_gv <- run_tvn(gv, "GSVA")
out <- rbind(r_ss, r_gv)
out <- merge(out, inv[, .(set_id, collection, n_genes_in_universe, source, note)], by="set_id", all.x=TRUE)
setorder(out, score, p_limma)
fwrite(out, file.path(RES,"metabolic_tumour_vs_normal.csv"))

## direction agreement between methods
w <- dcast(out, set_id ~ score, value.var=c("cohens_d","FDR_limma"))
cat("\nssGSEA vs GSVA Cohen's d correlation:",
    round(cor(w$cohens_d_ssGSEA, w$cohens_d_GSVA, use="complete.obs"),3),
    " sign agreement:", round(mean(sign(w$cohens_d_ssGSEA)==sign(w$cohens_d_GSVA)),3), "\n")
cat("both FDR<0.05:", sum(w$FDR_limma_ssGSEA<0.05 & w$FDR_limma_GSVA<0.05), "of", nrow(w), "\n")
fwrite(w, file.path(RES,"metabolic_tumour_vs_normal_method_concordance.csv"))

s <- out[score=="ssGSEA"]
cat("\nsignificant (limma FDR<0.05):", sum(s$FDR_limma<0.05), "of", nrow(s), "\n")
cat("\n--- TOP 25 UP IN TUMOUR (ssGSEA, by Cohen's d) ---\n")
print(s[FDR_limma<0.05][order(-cohens_d)][1:25, .(set_id, collection, n=n_genes_in_universe,
      d=round(cohens_d,3), t=round(t,1), FDR=signif(FDR_limma,3))])
cat("\n--- TOP 25 DOWN IN TUMOUR (ssGSEA, by Cohen's d) ---\n")
print(s[FDR_limma<0.05][order(cohens_d)][1:25, .(set_id, collection, n=n_genes_in_universe,
      d=round(cohens_d,3), t=round(t,1), FDR=signif(FDR_limma,3))])
cat("\n--- by collection: n significant / n tested (ssGSEA) ---\n")
print(s[, .(tested=.N, sig=sum(FDR_limma<0.05), up=sum(FDR_limma<0.05 & cohens_d>0),
            down=sum(FDR_limma<0.05 & cohens_d<0)), by=collection])
cat("DONE 22\n")
