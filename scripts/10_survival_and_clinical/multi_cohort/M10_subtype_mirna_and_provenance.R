## ==========================================================================
## M10 -- subtype behaviour of the measured miRNA itself (the direct analogue
## of MDA-MB-231 vs MCF7) + a provenance/verification file for the whole arm.
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(survival); library(metafor)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
MI <- readRDS(file.path(CA,"mirna_cohorts.rds"))
CO <- readRDS(file.path(CA,"mrna_cohorts.rds"))
SC <- readRDS(file.path(CA,"module_scores.rds"))

## ---- miR-130a level by subtype, every cohort that carries a subtype call ---
cat("\n=========== miR-130a level by subtype ===========\n")
rows <- list()
add <- function(...) rows[[length(rows)+1]] <<- data.table(...)
## TCGA (PAM50)
ph <- MI$TCGA_BRCA$pheno; x <- as.numeric(MI$TCGA_BRCA$M["hsa-miR-130a-3p",])
tb <- data.table(st=ph$pam50, v=x)[!is.na(st) & st!=""]
for(s in unique(tb$st)) add(cohort="TCGA_BRCA", subtype_source="PAM50 (RNAseq call)",
  subtype=s, n=sum(tb$st==s), mean=mean(tb$v[tb$st==s]), sd=sd(tb$v[tb$st==s]))
w <- wilcox.test(v ~ st, data=tb[st %in% c("Basal","LumA")])
cat(sprintf("TCGA    Basal(n=%d) mean %.3f vs LumA(n=%d) mean %.3f : diff %+.3f log2 RPM, wilcox p=%.3g\n",
  sum(tb$st=="Basal"), mean(tb$v[tb$st=="Basal"]), sum(tb$st=="LumA"),
  mean(tb$v[tb$st=="LumA"]), mean(tb$v[tb$st=="Basal"])-mean(tb$v[tb$st=="LumA"]), w$p.value))
BAS <- list(TCGA_BRCA=data.table(cohort="TCGA_BRCA", contrast="Basal - LumA",
  n_basal=sum(tb$st=="Basal"), n_luma=sum(tb$st=="LumA"),
  diff=mean(tb$v[tb$st=="Basal"])-mean(tb$v[tb$st=="LumA"]), wilcox_p=w$p.value))
## GSE19783 (subtype field)
ph <- MI$GSE19783$pheno; x <- as.numeric(MI$GSE19783$M["hsa-miR-130a",])
tb <- data.table(st=ph$subtype, v=x)[!is.na(st) & st!=""]
print(tb[, .(n=.N, mean=round(mean(v),3)), by=st][order(-mean)])
for(s in unique(tb$st)) add(cohort="GSE19783", subtype_source="deposited subtype call",
  subtype=s, n=sum(tb$st==s), mean=mean(tb$v[tb$st==s]), sd=sd(tb$v[tb$st==s]))
if(all(c("Basal","Lum A") %in% tb$st)){
  s2 <- tb[st %in% c("Basal","Lum A")]
  w <- wilcox.test(v ~ st, data=s2)
  d <- mean(s2$v[s2$st=="Basal"]) - mean(s2$v[s2$st=="Lum A"])
  cat(sprintf("GSE19783 Basal(n=%d) vs LumA(n=%d): diff %+.3f, wilcox p=%.3g\n",
      sum(s2$st=="Basal"), sum(s2$st=="Lum A"), d, w$p.value))
  BAS$GSE19783 <- data.table(cohort="GSE19783", contrast="Basal - Lum A",
    n_basal=sum(s2$st=="Basal"), n_luma=sum(s2$st=="Lum A"), diff=d, wilcox_p=w$p.value)
}
## ER-defined contrast wherever ER is recorded
for(nm in c("GSE22216","GSE26659","GSE45666","GSE40525")){
  o <- MI[[nm]]; r <- rownames(o$M)[tolower(rownames(o$M))=="hsa-mir-130a"][1]
  er <- o$pheno$ER
  ern <- ifelse(er %in% c("1","positive","Positive","pos"), "ERpos",
         ifelse(er %in% c("0","negative","Negative","neg"), "ERneg", NA))
  if(sum(!is.na(ern))<20) next
  v <- as.numeric(o$M[r,])
  a <- v[ern=="ERneg" & !is.na(ern)]; b <- v[ern=="ERpos" & !is.na(ern)]
  if(length(a)<8 || length(b)<8) next
  w <- wilcox.test(a,b)
  cat(sprintf("%-9s ER-neg(n=%d) vs ER-pos(n=%d): diff %+.3f, wilcox p=%.3g\n",
      nm, length(a), length(b), mean(a)-mean(b), w$p.value))
  BAS[[nm]] <- data.table(cohort=nm, contrast="ER-negative - ER-positive",
    n_basal=length(a), n_luma=length(b), diff=mean(a)-mean(b), wilcox_p=w$p.value)
  add(cohort=nm, subtype_source="ER IHC", subtype="ERneg", n=length(a), mean=mean(a), sd=sd(a))
  add(cohort=nm, subtype_source="ER IHC", subtype="ERpos", n=length(b), mean=mean(b), sd=sd(b))
}
LV <- rbindlist(rows); fwrite(LV, file.path(OUT,"multicohort_E_mir130a_level_by_subtype.csv"))
BASd <- rbindlist(BAS)
print(BASd)
fwrite(BASd, file.path(OUT,"multicohort_E_mir130a_basal_vs_luminal.csv"))
cat("\nCCLE reference (results/v3/mir130a_ccle_breast_lines.csv):\n")
CC <- fread("results/v3/mir130a_ccle_breast_lines.csv")
cat("columns:", paste(names(CC), collapse=","), " rows:", nrow(CC), "\n")
nmc <- names(CC)[1]
print(CC[grepl("MDAMB231|MDA-MB-231|MCF7", get(nmc), ignore.case=TRUE)])

## ---- miR-130a Cox within subtype / ER stratum ------------------------------
cat("\n=========== miR-130a Cox within subtype ===========\n")
srows <- list()
zs <- function(x) as.numeric(scale(x))
## TCGA, PAM50 strata
ph <- MI$TCGA_BRCA$pheno; x <- as.numeric(MI$TCGA_BRCA$M["hsa-miR-130a-3p",])
for(st in c("Basal","LumA")) for(ep in c("OS","PFI")){
  sel <- ph$pam50 == st & !is.na(ph$pam50)
  d <- data.frame(time=ph[[paste0(ep,"_time")]][sel], event=ph[[paste0(ep,"_event")]][sel],
                  z=zs(x[sel]))
  d <- d[complete.cases(d) & d$time>0,]
  if(nrow(d)<40 || sum(d$event)<8) next
  s <- summary(coxph(Surv(time,event)~z, data=d))$coefficients["z",]
  srows[[length(srows)+1]] <- data.table(cohort="TCGA_BRCA", stratum=st,
    stratum_source="PAM50", endpoint=ep, n=nrow(d), nevent=sum(d$event),
    HR=exp(s[1]), se_logHR=s[3], lo=exp(s[1]-1.96*s[3]), hi=exp(s[1]+1.96*s[3]), p=s[5])
}
## GSE22216, ER strata
o <- MI$GSE22216; ph <- o$pheno; x <- as.numeric(o$M["hsa-miR-130a",])
for(st in c("0","1")){
  sel <- ph$ER == st & !is.na(ph$ER)
  d <- data.frame(time=ph$DRFS_time[sel], event=ph$DRFS_event[sel], z=zs(x[sel]))
  d <- d[complete.cases(d) & d$time>0,]
  if(nrow(d)<40 || sum(d$event)<8) next
  s <- summary(coxph(Surv(time,event)~z, data=d))$coefficients["z",]
  srows[[length(srows)+1]] <- data.table(cohort="GSE22216",
    stratum=ifelse(st=="0","ER-negative","ER-positive"), stratum_source="ER IHC",
    endpoint="DRFS", n=nrow(d), nevent=sum(d$event), HR=exp(s[1]), se_logHR=s[3],
    lo=exp(s[1]-1.96*s[3]), hi=exp(s[1]+1.96*s[3]), p=s[5])
}
SR <- rbindlist(srows)
print(SR[, .(cohort,stratum,endpoint,n,nevent,HR=round(HR,3),lo=round(lo,3),
             hi=round(hi,3),p=signif(p,3))])
fwrite(SR, file.path(OUT,"multicohort_E_mir130a_cox_by_subtype.csv"))

## ---- provenance / verification --------------------------------------------
cat("\n=========== provenance ===========\n")
prov <- rbindlist(list(
 data.table(item="miRNA cohorts assembled", value="7", detail="TCGA-BRCA, GSE19783, GSE22216, GSE37405, GSE26659, GSE40525, GSE45666"),
 data.table(item="miRNA cohorts with an outcome", value="4 Cox + 1 binary",
   detail="TCGA(OS/PFI/DSS), GSE19783(OS/BCSS), GSE22216(DRFS), GSE37405(RFS/OS); GSE26659 has only a binary 72-month relapse flag"),
 data.table(item="mRNA cohorts assembled", value="8",
   detail="TCGA-BRCA, METABRIC, SCAN-B/GSE96058, GSE20685, GSE21653, GSE22219, GSE58812, GTEx breast"),
 data.table(item="mRNA cohorts with an outcome", value="7", detail="all except GTEx"),
 data.table(item="patients in the module meta-analysis", value="7231",
   detail="primary endpoint per cohort, 1903 events"),
 data.table(item="patients in the miR-130a HR meta-analysis", value="1512",
   detail="4 cohorts, 344 events"),
 data.table(item="SCAN-B expression slice", value="1897 of 1901 wanted symbols",
   detail="re-streamed from GSE96058_expr.csv.gz; the earlier v2 slice held only 173 symbols; the 172 shared genes agree to 0 (max abs difference)"),
 data.table(item="Rule 4 guard (limma symbol rownames)", value="passed",
   detail="results/BRCA_DEX_genes.csv feature column asserted >95% symbol-like before use (0.9986); expression matrices asserted free of duplicated rownames"),
 data.table(item="meta-analysis method", value="DerSimonian-Laird random effects",
   detail="own implementation, cross-checked against metafor::rma(method='DL'); max abs difference recorded per pooled row in metafor_max_abs_diff"),
 data.table(item="GTEx cross-study hub DE concordance", value="33/44 = 75.0%, binom p=0.00126",
   detail="independently re-derived here; reproduces the v2 workflow's number exactly"),
 data.table(item="GSE22219<->GSE22216 patient linkage", value="205 of 210 (fingerprint)",
   detail="naive title-index matching agrees on only 1.9-51.9% of clinical fields and was rejected")))
fwrite(prov, file.path(OUT,"multicohort_PROVENANCE.csv"))
print(prov)
