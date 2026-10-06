## ==========================================================================
## M8_analysisA_forest.R
##  A3) hub DE direction: TCGA tumour vs GTEx normal breast (cross-study),
##      re-derived here rather than taken from the v2 workflow
##  A4) consolidate every direction-concordance test into one table
##  Z ) one forest-plot-ready table covering C, D and E
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
S  <- readRDS(file.path(CA,"genesets.rds"))
CO <- readRDS(file.path(CA,"mrna_cohorts.rds"))

## ---------------- A3: TCGA tumour vs GTEx normal ---------------------------
cat("\n=========== A3  TCGA tumour vs GTEx normal (cross-study) ===========\n")
DE <- fread("results/BRCA_DEX_genes.csv")
stopifnot(mean(grepl("^[A-Za-z]", DE$feature)) > 0.95)     ## Rule 4 guard
fc <- setNames(DE$logFC, DE$feature); q <- setNames(DE$adj.P.Val, DE$feature)
Xt <- CO$TCGA_BRCA$X; Xg <- CO$GTEx_breast$X
sh <- intersect(rownames(Xt), rownames(Xg))
cat("genes shared TCGA/GTEx:", length(sh), "\n")
pct <- function(X){ apply(X, 2, function(v) rank(v, na.last="keep")/sum(is.finite(v))) }
Pt <- pct(Xt[sh,,drop=FALSE]); Pg <- pct(Xg[sh,,drop=FALSE])
dpc <- rowMeans(Pt, na.rm=TRUE) - rowMeans(Pg, na.rm=TRUE)
A3 <- list()
for(lab in c("network_hubs","all_network_nodes")){
  gs <- if(lab=="network_hubs") S$HUB_PROTEIN else S$ALL_NETWORK_PROTEIN
  g <- intersect(gs, sh); g <- g[g %in% names(fc)]
  g <- g[is.finite(q[g]) & q[g] < 0.05]
  n <- length(g); conc <- sum(sign(dpc[g]) == sign(fc[g]))
  bt <- binom.test(conc, n, 0.5)
  A3[[lab]] <- data.table(cohort="GTEx_breast", layer="mRNA", feature_set=lab,
    n_tested=n, n_concordant=conc, pct=round(100*conc/n,2), binom_p=bt$p.value,
    n_normal=ncol(Xg), n_tumour=ncol(Xt),
    spearman_logFC=cor(dpc[g], fc[g], method="spearman"))
  cat(sprintf("  %-18s %3d/%3d (%.1f%%) p=%.3g rho=%.3f\n", lab, conc, n,
      100*conc/n, bt$p.value, cor(dpc[g], fc[g], method="spearman")))
}
A3 <- rbindlist(A3)
A3[, `:=`(note="cross-study / cross-platform: within-sample rank percentiles, TCGA tumours vs GTEx normal breast")]

## ---------------- A4: consolidate ------------------------------------------
A <- fread(file.path(OUT,"multicohort_A_DE_direction_concordance.csv"))
A[, note := "tumour vs normal within the same GEO series / TCGA"]
ACC <- c(GSE42568="GSE42568", GSE45827="GSE45827", GSE10780="GSE10780",
  TCGA_BRCA="TCGA-BRCA", GSE26659="GSE26659", GSE40525="GSE40525",
  GSE45666="GSE45666", GTEx_breast="GTEx v8")
A <- rbind(A, A3, fill=TRUE)
A[, accession := ACC[cohort]]
## cohorts that CANNOT support analysis A, stated explicitly
noDE <- data.table(cohort=c("METABRIC","GSE96058","GSE20685","GSE21653","GSE22219",
                            "GSE58812","GSE19783","GSE22216","GSE37405"),
  accession=c("cBioPortal brca_metabric","GSE96058","GSE20685","GSE21653","GSE22219",
              "GSE58812","GSE19783","GSE22216","GSE37405"),
  layer=c(rep("mRNA",6),"mRNA+miRNA","miRNA","miRNA"),
  feature_set="NOT_POSSIBLE_no_normal_samples", n_tested=NA_integer_,
  n_concordant=NA_integer_, pct=NA_real_, binom_p=NA_real_,
  n_normal=0L, n_tumour=c(1980L,3273L,327L,266L,216L,107L,101L,210L,153L),
  spearman_logFC=NA_real_,
  note="tumour-only series: no matched normal, so no tumour-vs-normal DE direction test is possible; the v2 workflow's hub co-expression sign concordance (results/v2/newcohorts_summary.csv) is the substitute for these cohorts")
A <- rbind(A, noDE, fill=TRUE)
setcolorder(A, c("cohort","accession","layer","feature_set","n_tested","n_concordant",
                 "pct","binom_p","spearman_logFC","n_tumour","n_normal","note"))
fwrite(A, file.path(OUT,"multicohort_A_DE_direction_concordance.csv"))
cat("\n--- consolidated analysis A ---\n")
print(A[!is.na(pct), .(cohort,layer,feature_set,n_concordant,n_tested,pct,
                       binom_p=signif(binom_p,3))])

## ---------------- Z: forest-plot-ready table -------------------------------
cat("\n=========== Z  forest-plot-ready table ===========\n")
cols <- c("analysis","feature","stratum","row_type","cohort","accession","endpoint",
          "model","n","nevent","estimate","lo","hi","p","weight_pct","k","I2","tau2",
          "Q","p_Q","scale","note")
mk <- function(dt) { for(cc in setdiff(cols, names(dt))) dt[[cc]] <- NA; dt[, ..cols] }

CD <- fread(file.path(OUT,"multicohort_CD_module_meta.csv"))
CD[, `:=`(analysis=paste0(analysis,"__",pool), feature=module, stratum="all",
          estimate=HR, scale="hazard ratio per SD of the module score")]
E  <- fread(file.path(OUT,"multicohort_E_subtype_meta.csv"))
E[, `:=`(analysis="module_HR_subtype_stratified", feature=module, stratum=stratum,
  row_type="RE_pooled", cohort="POOLED (DerSimonian-Laird)", accession="",
  endpoint="primary per cohort", model="univariate", estimate=HR, weight_pct=100,
  tau2=NA_real_, scale="hazard ratio per SD of the module score",
  note=paste0("cohorts: ", cohorts))]
Ep <- fread(file.path(OUT,"multicohort_E_subtype_cox_percohort.csv"))
Ep[, `:=`(analysis="module_HR_subtype_stratified", feature=module, row_type="study",
  accession="", model="univariate", estimate=HR, weight_pct=NA_real_,
  scale="hazard ratio per SD of the module score", note=stratum_source)]
H <- fread(file.path(OUT,"multicohort_D_3node_vs_higherorder_headtohead.csv"))
H[, `:=`(analysis="FFL_3node_vs_higherorder_mutually_adjusted",
  feature=ifelse(term=="z3","FFL_3NODE_UNION","FFL_HIGHER_ONLY"), stratum="all",
  row_type=ifelse(is.na(cohort),"RE_pooled","study"),
  cohort=ifelse(is.na(cohort),"POOLED (DerSimonian-Laird)",cohort),
  accession="", model="mutually adjusted", estimate=HR,
  scale="hazard ratio per SD, both module scores in one Cox model",
  note="the two scores are correlated r=0.72-0.88 within cohort; this is a collinear contrast, not two independent effects")]
Z <- rbindlist(list(mk(CD), mk(E), mk(Ep), mk(H)), use.names=TRUE)
Z <- Z[is.finite(estimate)]
fwrite(Z, file.path(OUT,"multicohort_FOREST_TABLE.csv"))
cat("forest table rows:", nrow(Z), " analyses:", length(unique(Z$analysis)), "\n")

## ---------------- headline summary ----------------------------------------
P <- Z[row_type=="RE_pooled"]
head1 <- P[analysis %in% c("module_HR_primary_univariate__all_cohorts",
  "module_HR_primary_univariate__independent_of_TCGA",
  "module_HR_primary_adjusted__all_cohorts",
  "module_HR_subtype_stratified","FFL_3node_vs_higherorder_mutually_adjusted")]
setorder(head1, analysis, p)
fwrite(head1, file.path(OUT,"multicohort_HEADLINE_POOLED.csv"))
cat("\nWROTE", file.path(OUT,"multicohort_HEADLINE_POOLED.csv"), nrow(head1), "pooled rows\n")
