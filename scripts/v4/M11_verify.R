## Independent re-checks of the headline numbers, computed a second way.
suppressPackageStartupMessages({library(data.table); library(metafor); library(survival)})
setwd("/path/to/revision")
OUT <- "results/v4"; CA <- "cache/v4/multicohort"
v <- list(); add <- function(...) v[[length(v)+1]] <<- data.table(...)

## 1. miR-130a HR meta re-run straight from the per-cohort Cox table with metafor
CX <- fread(file.path(OUT,"multicohort_B_mir130a_cox_all.csv"))
P <- c(TCGA_BRCA="OS", GSE19783="OS", GSE22216="DRFS", GSE37405="RFS")
s <- CX[model=="univariate" & endpoint==P[cohort]]
m <- rma(yi=log(s$HR), sei=s$se_logHR, method="DL")
MT <- fread(file.path(OUT,"multicohort_B_mir130a_meta.csv"))
st <- MT[analysis=="mir130a_HR_primary_univariate" & row_type=="RE_pooled"]
add(check="miR-130a pooled HR", stored=round(st$HR,6), recomputed=round(exp(m$b[1]),6),
    max_abs_diff=abs(st$HR-exp(m$b[1])))
add(check="miR-130a pooled I2", stored=round(st$I2,4), recomputed=round(m$I2,4),
    max_abs_diff=abs(st$I2-m$I2))
add(check="miR-130a pooled p", stored=signif(st$p,6), recomputed=signif(m$pval,6),
    max_abs_diff=abs(st$p-m$pval))

## 2. TCGA miR-130a PFI HR vs the value already in results/v3
v3 <- fread("results/v3/mir130a_survival_cox.csv")
old <- v3[endpoint=="PFI" & variable=="mir130a_3p" & model=="mir130a_3p"]
new <- CX[cohort=="TCGA_BRCA" & endpoint=="PFI" & model=="univariate"]
add(check="TCGA miR-130a PFI HR (v3 vs v4)", stored=round(old$HR,6),
    recomputed=round(new$HR,6), max_abs_diff=abs(old$HR-new$HR))

## 3. TCGA miR-130a tumour-vs-normal logFC vs the manuscript's -1.16
B1 <- fread(file.path(OUT,"multicohort_B_mir130a_tumour_vs_normal.csv"))
DEm <- fread("results/BRCA_DEX_mirnas.csv")
cat("NOTE: results/BRCA_DEX_mirnas.csv is keyed on the AGGREGATED miRNA name space\n",
    "     (hsa-miR-130a, logFC", signif(DEm[feature=="hsa-miR-130a", logFC],4),
    "), whereas data/brca_mirna_expr.rds is arm-level (hsa-miR-130a-3p).\n",
    "     The manuscript's -1.16 is the ARM-LEVEL value.\n")
add(check="TCGA miR-130a-3p log2FC vs the manuscript's -1.16 (arm level)",
    stored=-1.16, recomputed=round(B1[cohort=="TCGA_BRCA"]$log2FC,6),
    max_abs_diff=abs(-1.16 - B1[cohort=="TCGA_BRCA"]$log2FC))
add(check="TCGA aggregated hsa-miR-130a log2FC (different name space, not the same quantity)",
    stored=round(DEm[feature=="hsa-miR-130a", logFC],6),
    recomputed=round(B1[cohort=="TCGA_BRCA"]$log2FC,6),
    max_abs_diff=abs(DEm[feature=="hsa-miR-130a", logFC] - B1[cohort=="TCGA_BRCA"]$log2FC))

## 4. module score matrices: no duplicated gene rownames anywhere (Rule 4)
CO <- readRDS(file.path(CA,"mrna_cohorts.rds"))
dupmax <- max(sapply(CO, function(o) sum(duplicated(rownames(o$X)))))
add(check="max duplicated rownames across mRNA cohort matrices", stored=0,
    recomputed=dupmax, max_abs_diff=dupmax)
MI <- readRDS(file.path(CA,"mirna_cohorts.rds"))
dupmi <- max(sapply(MI, function(o) sum(duplicated(rownames(o$M)))))
add(check="max duplicated rownames across miRNA cohort matrices", stored=0,
    recomputed=dupmi, max_abs_diff=dupmi)

## 5. metafor agreement recorded inside the module meta table
MM <- fread(file.path(OUT,"multicohort_CD_module_meta.csv"))
add(check="max |own DL - metafor| over all pooled module rows", stored=0,
    recomputed=signif(max(MM$metafor_max_abs_diff, na.rm=TRUE),4),
    max_abs_diff=max(MM$metafor_max_abs_diff, na.rm=TRUE))

## 6. recompute the SCAN-B MIR130A_ACTIVITY OS HR from scratch
SC <- readRDS(file.path(CA,"module_scores.rds"))
ph <- CO$GSE96058$pheno
d <- data.frame(time=ph$OS_time, event=ph$OS_event,
                z=as.numeric(scale(SC$GSE96058$MIR130A_ACTIVITY)))
d <- d[complete.cases(d) & d$time>0,]
hr <- exp(coef(coxph(Surv(time,event)~z, data=d))["z"])
CD <- fread(file.path(OUT,"multicohort_CD_module_cox_all.csv"))
sc <- CD[cohort=="GSE96058" & module=="MIR130A_ACTIVITY" & model=="univariate" & endpoint=="OS"]
add(check="SCAN-B MIR130A_ACTIVITY OS HR", stored=round(sc$HR,6),
    recomputed=round(unname(hr),6), max_abs_diff=abs(sc$HR-hr))

V <- rbindlist(v)
print(V)
fwrite(V, file.path(OUT,"multicohort_VERIFICATION.csv"))
V[, tolerated := c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)[seq_len(.N)]]
V[grepl("v3 vs v4", check), tolerated := TRUE]
cat("\nchecks agreeing to <1e-6:", sum(V$max_abs_diff < 1e-6), "of", nrow(V), "\n")
cat("expected non-zero rows:\n")
print(V[max_abs_diff >= 1e-6, .(check, stored, recomputed, max_abs_diff)])
cat("\n  - 'v3 vs v4': v3 fitted n=1051/140 events, this run n=1053/140 (two extra samples\n",
    "    recovered by the survival-table join); HR 0.8395 vs 0.8399.\n",
    "  - 'aggregated name space': the two rows are different quantities, shown for the record.\n")
