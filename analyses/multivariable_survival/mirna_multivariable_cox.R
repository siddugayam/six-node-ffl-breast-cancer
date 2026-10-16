## mirna_multivariable_cox.R -- unadjusted and covariate-adjusted Cox models for the 11 miRNAs of the survival
## meta-analysis (Table S8), Schoenfeld tests, fixed- and random-effects pooling, BH over the 11 miRNAs, the
## harmonised overall-survival sensitivity and the pool without the stromal score (see SETTINGS.md).
## Cohort loading, pick() and the within-cohort z-scoring are copied from
## scripts/10_survival_and_clinical/mirna_survival_meta_analysis/c11_mirna_cox_meta.R.
## Root of the repository: environment variable FFL_ROOT. Outputs go to the folder of this script; no sample-level
## values are written.
## usage: FFL_ROOT=<repository checkout> Rscript mirna_multivariable_cox.R
suppressPackageStartupMessages({library(survival); library(metafor); library(data.table)})
ROOT <- Sys.getenv("FFL_ROOT"); stopifnot(nzchar(ROOT))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
setwd(ROOT)
B <- "cache/tissue_and_cohorts/cohorts/built"

TARGETS <- c("miR-21","miR-195","miR-204","miR-383","miR-124","miR-155","miR-429",
             "miR-141","miR-34a","miR-101","miR-29a")
RX <- c("miR-21"="^hsa-mir-21(-[0-9])?(-[35]p)?$",
        "miR-195"="^hsa-mir-195(-[35]p)?$", "miR-204"="^hsa-mir-204(-[35]p)?$",
        "miR-383"="^hsa-mir-383(-[35]p)?$", "miR-124"="^hsa-mir-124a?(-[123])?(-[35]p)?$",
        "miR-155"="^hsa-mir-155(-[35]p)?$", "miR-429"="^hsa-mir-429(-[35]p)?$",
        "miR-141"="^hsa-mir-141(-[35]p)?$", "miR-34a"="^hsa-mir-34a(-[35]p)?$",
        "miR-101"="^hsa-mir-101(-[12])?(-[35]p)?$",
        "miR-29a"="^hsa-mir-29a(-[35]p)?$")
pick <- function(M, tgt){
  rn <- rownames(M)
  nm <- tolower(gsub("miR","mir",rn))
  star <- grepl("\\*", nm)
  nm2 <- gsub("\\*","",nm)
  hit <- which(grepl(RX[[tgt]], nm2) & !star)
  if(!length(hit)) hit <- which(grepl(RX[[tgt]], nm2))
  if(!length(hit)) return(NULL)
  mu <- rowMeans(M[hit,,drop=FALSE], na.rm=TRUE)
  mu[!is.finite(mu)] <- -Inf
  if(all(!is.finite(mu))) return(NULL)
  i <- hit[which.max(mu)]
  list(row=as.numeric(M[i,]), name=rn[i], n_matched=length(hit),
       all=paste(rn[hit], collapse=";"))
}

## stromal score (CAF-A): mean row-z over the signature genes present (score_of of N2_cohort_analysis.R)
library(matrixStats)
zrow <- function(M){ s <- rowSds(M, na.rm=TRUE); s[s==0] <- NA; (M - rowMeans(M, na.rm=TRUE))/s }
score_of <- function(X, genes){
  g <- intersect(genes, rownames(X)); if(length(g) < 2) return(list(v=NULL, n=length(g)))
  Z <- zrow(X[g,,drop=FALSE]); Z <- Z[rowSums(is.na(Z))==0,,drop=FALSE]
  list(v=colMeans(Z), n=nrow(Z))
}
cafA_genes <- readLines("cache/external_cohorts/cafA_signature_genes.txt")
stopifnot(length(cafA_genes) == 123)

## ------------------------------------------------------------------ cohorts --
CO <- readRDS("cache/metabolism/multicohort/mirna_cohorts.rds")
COH <- list()
add <- function(name, M, ph, endpoint, os=NULL, strata=NULL, covars, covar_note){
  COH[[name]] <<- list(name=name, M=M, ph=ph, endpoint=endpoint, os=os, strata=strata, covars=covars,
                       covar_note=covar_note)
}
fac <- function(x, lev) factor(ifelse(x %in% lev, x, NA), levels=lev)

## TCGA-BRCA
p <- CO$TCGA_BRCA$pheno
cl <- fread("data/brca_clinicalMatrix.tsv", sep="\t", quote="")
cl <- cl[match(p$sample, cl$sampleID)]; stopifnot(identical(cl$sampleID, p$sample))
st <- sub("^Stage ", "", p$stage); st <- sub("[ABC]$", "", st)
p$stage4 <- fac(st, c("I","II","III","IV"))
p$ER_ihc <- fac(cl$breast_carcinoma_estrogen_receptor_status, c("Positive","Negative"))
p$PAM50 <- fac(cl$PAM50Call_RNAseq, c("LumA","LumB","Her2","Basal","Normal"))
fib <- fread("results/systems_and_dynamics/deconv_fibroblast_estimates_extended.csv")
p$stromal <- fib$meanZ_CAF_A_estimate[match(substr(p$sample,1,15), fib$sample)]
p$age <- as.numeric(p$age)
add("TCGA_BRCA", CO$TCGA_BRCA$M, p, c("OS_time","OS_event"), os=c("OS_time","OS_event"),
    covars=c("age","stage4","ER_ihc","PAM50","stromal"), covar_note="age; stage I-IV; ER (IHC); PAM50; stromal score")

## GSE19783
p <- CO$GSE19783$pheno
g <- readRDS("cache/external2/gse19783.rds")
stopifnot(identical(colnames(g$GSE19783_mRNA), colnames(CO$GSE19783$M)), identical(p$sample, colnames(CO$GSE19783$M)))
sc <- score_of(g$GSE19783_mRNA, cafA_genes); n_caf_19783 <- sc$n
p$stromal <- as.numeric(sc$v)
p$ER2 <- fac(p$ER, c("Positive","Negative"))
p$subtype5 <- fac(p$subtype, c("Lum A","Lum B","ERBB2","Basal","Normal-like"))
add("GSE19783", CO$GSE19783$M, p, c("OS_time","OS_event"), os=c("OS_time","OS_event"),
    covars=c("ER2","subtype5","stromal"), covar_note="ER; PAM50 subtype; stromal score")

## GSE22216 (matched mRNA: GSE22219)
p <- CO$GSE22216$pheno
cb <- readRDS("cache/external_cohorts/cohorts_base.rds")
lk <- read.csv("cache/metabolism/multicohort/GSE22219_GSE22216_link.csv", stringsAsFactors=FALSE)
sc <- score_of(cb$GSE22219$X, cafA_genes); n_caf_22219 <- sc$n
p$stromal <- as.numeric(sc$v[lk$gsm_mrna[match(p$gsm, lk$gsm_mirna)]])
p$grade3 <- factor(p$grade, levels=c(1,2,3))
p$age <- as.numeric(p$age); p$size <- as.numeric(p$size); p$node <- as.numeric(p$node); p$ER <- as.numeric(p$ER)
add("GSE22216", CO$GSE22216$M, p, c("DRFS_time","DRFS_event"),
    covars=c("age","grade3","size","node","ER","stromal"), covar_note="age; grade; size; nodes; ER; stromal score")

## GSE37405
p <- CO$GSE37405$pheno
for(v in c("age","size","node")) p[[v]] <- as.numeric(p[[v]])
add("GSE37405", CO$GSE37405$M, p, c("RFS_time","RFS_event"), os=c("OS_time","OS_event"), strata="platform",
    covars=c("age","size","node"), covar_note="age; size; nodes (ER-positive cohort; platform stratum)")

## GSE59829
E <- as.matrix(fread(file.path(B,"GSE59829_expr.csv")), rownames=1)
P <- fread(file.path(B,"GSE59829_pheno.csv")); stopifnot(identical(colnames(E), P$gsm))
ph <- data.frame(sample=P$gsm,
                 DMFS_time=suppressWarnings(as.numeric(P[["dfs (months)"]])),
                 DMFS_event=suppressWarnings(as.numeric(P[["distant_metastasis"]])),
                 age=suppressWarnings(as.numeric(P[["age"]])),
                 size=suppressWarnings(as.numeric(P[["tumor_size"]])),
                 ER=P[["esr1_status"]], stringsAsFactors=FALSE)
ph$ER2 <- factor(ph$ER)
add("GSE59829", E, ph, c("DMFS_time","DMFS_event"), covars=c("age","size","ER2"), covar_note="age; size; ER")

## GSE78870
E <- as.matrix(fread(file.path(B,"GSE78870_expr.csv")), rownames=1)
P <- fread(file.path(B,"GSE78870_pheno.csv")); stopifnot(identical(colnames(E), P$gsm))
ph <- data.frame(sample=P$gsm,
                 TTP_time=suppressWarnings(as.numeric(P[["time-to-progression (in months)"]])),
                 TTP_event=suppressWarnings(as.numeric(P[["disease progression (event)"]])),
                 age=suppressWarnings(as.numeric(P[["age at start therapy"]])),
                 grade3=factor(P[["grade"]], levels=c(1,2,3)),
                 ER2=fac(P[["estrogen receptor status"]], c("positive","negative")), stringsAsFactors=FALSE)
add("GSE78870", E, ph, c("TTP_time","TTP_event"), covars=c("age","grade3","ER2"), covar_note="age; grade; ER")

## covariate availability (counts only)
cov_rows <- list()
for(cn in names(COH)){
  co <- COH[[cn]]
  for(v in c(co$covars)) cov_rows[[length(cov_rows)+1]] <- data.frame(cohort=cn, covariate=v,
      n_samples=nrow(co$ph), n_available=sum(!is.na(co$ph[[v]])),
      levels=if(is.factor(co$ph[[v]])) paste(sprintf("%s %d", levels(co$ph[[v]]), as.integer(table(co$ph[[v]]))), collapse="; ") else "")
}
COVTAB <- rbindlist(cov_rows)
COVTAB <- rbind(COVTAB, data.table(cohort=c("GSE19783","GSE22216"), covariate="stromal score genes present",
                                   n_samples=NA, n_available=c(n_caf_19783, n_caf_22219), levels="of 123"))
fwrite(COVTAB, file.path(HERE, "cohort_covariates.csv"))

## ---------------------------------------------------------------- Cox fits --
fit1 <- function(co, z, endp, covars){
  ph <- co$ph
  d <- data.frame(z=z, time=suppressWarnings(as.numeric(ph[[endp[1]]])), ev=suppressWarnings(as.numeric(ph[[endp[2]]])))
  if(!is.null(co$strata)) d$stratum <- ph[[co$strata]]
  for(cv in covars) d[[cv]] <- ph[[cv]]
  d <- d[complete.cases(d) & is.finite(d$time) & d$time > 0, , drop=FALSE]
  for(cv in covars) if(is.factor(d[[cv]])) d[[cv]] <- droplevels(d[[cv]])
  if(nrow(d) < 25 || sum(d$ev) < 8) return(NULL)
  f <- paste0("Surv(time, ev) ~ ", paste(c("z", covars), collapse=" + "),
              if(!is.null(co$strata)) " + strata(stratum)" else "")
  fit <- try(coxph(as.formula(f), data=d), silent=TRUE)
  if(inherits(fit, "try-error")) return(NULL)
  s <- summary(fit)$coefficients
  zp <- try(cox.zph(fit, transform="km"), silent=TRUE)
  zg <- if(inherits(zp, "try-error")) NA else zp$table["GLOBAL","p"]
  zz <- if(inherits(zp, "try-error")) NA else zp$table["z","p"]
  data.frame(n=nrow(d), nevent=sum(d$ev), logHR=s["z","coef"], se=s["z","se(coef)"], HR=exp(s["z","coef"]),
             lo=exp(s["z","coef"]-1.96*s["z","se(coef)"]), hi=exp(s["z","coef"]+1.96*s["z","se(coef)"]),
             p=s["z","Pr(>|z|)"], schoenfeld_global_p=zg, schoenfeld_mirna_p=zz,
             covariates=if(length(covars)) paste(covars, collapse="+") else "none", stringsAsFactors=FALSE)
}

rows <- list()
for(cn in names(COH)){
  co <- COH[[cn]]; M <- co$M; stopifnot(ncol(M)==nrow(co$ph))
  for(tg in TARGETS){
    pk <- pick(M, tg); if(is.null(pk)) next
    x <- pk$row; z <- rep(NA_real_, length(x))
    if(!is.null(co$strata)){ s <- co$ph[[co$strata]]; for(lv in unique(s)){ i <- which(s==lv); z[i] <- as.numeric(scale(x[i])) } }
    else z <- as.numeric(scale(x))
    eps <- list(primary=co$endpoint); if(!is.null(co$os)) eps$OS <- co$os
    for(en in names(eps)){
      cvn <- setdiff(co$covars, "stromal")
      for(model in c("unadjusted","adjusted","adjusted_without_stromal")){
        cv <- switch(model, unadjusted=character(0), adjusted=co$covars, adjusted_without_stromal=cvn)
        if(model=="adjusted_without_stromal" && !("stromal" %in% co$covars)) next
        r <- fit1(co, z, eps[[en]], cv)
        if(is.null(r)) next
        rows[[length(rows)+1]] <- cbind(data.frame(cohort=cn, miRNA=tg, probe=pk$name,
          endpoint=if(en=="primary") sub("_time$","",co$endpoint[1]) else "OS", endpoint_set=en, model=model,
          covariate_set=if(model=="unadjusted") "none" else if(model=="adjusted") co$covar_note else sub("; stromal score","",co$covar_note),
          stringsAsFactors=FALSE), r)
      }
    }
  }
}
COX <- rbindlist(rows)
fwrite(COX, file.path(HERE, "percohort_cox.csv"))

## ---------------------------------------------------------------- pooling ---
pool <- function(df, label){
  out <- list()
  for(tg in TARGETS){
    s <- df[miRNA==tg]; if(nrow(s) < 2) next
    re <- rma(yi=s$logHR, sei=s$se, method="DL"); fe <- rma(yi=s$logHR, sei=s$se, method="FE")
    out[[length(out)+1]] <- data.frame(analysis=label, miRNA=tg, k=nrow(s), cohorts=paste(s$cohort, collapse=";"),
      n=sum(s$n), nevent=sum(s$nevent),
      HR_FE=exp(fe$b[1]), lo_FE=exp(fe$ci.lb), hi_FE=exp(fe$ci.ub), p_FE=fe$pval,
      HR_RE=exp(re$b[1]), lo_RE=exp(re$ci.lb), hi_RE=exp(re$ci.ub), p_RE=re$pval,
      I2=re$I2, tau2=re$tau2, Q=re$QE, p_Q=re$QEp, stringsAsFactors=FALSE)
  }
  o <- rbindlist(out); o[, q_BH_FE := p.adjust(p_FE, "BH")]; o[, q_BH_RE := p.adjust(p_RE, "BH")]; o
}
P_UNADJ <- pool(COX[endpoint_set=="primary" & model=="unadjusted"], "primary endpoint, unadjusted")
P_ADJ   <- pool(COX[endpoint_set=="primary" & model=="adjusted"], "primary endpoint, adjusted (maximal available)")
P_NOSTR <- pool(COX[endpoint_set=="primary" & (model=="adjusted_without_stromal" |
                     (model=="adjusted" & !grepl("stromal", covariate_set)))], "primary endpoint, adjusted without the stromal score")
P_OSU   <- pool(COX[endpoint_set=="OS" & model=="unadjusted"], "OS only, unadjusted")
P_OSA   <- pool(COX[endpoint_set=="OS" & model=="adjusted"], "OS only, adjusted (maximal available)")
fwrite(rbind(P_UNADJ, P_ADJ), file.path(HERE, "pooled_adjusted.csv"))
fwrite(P_NOSTR, file.path(HERE, "pooled_adjusted_without_stromal.csv"))
fwrite(rbind(P_OSU, P_OSA), file.path(HERE, "pooled_os_harmonised.csv"))

## ---------------------------------------------- Table S8 reproduction (control) --
S8 <- fread("supplementary_tables/TableS8_mirna_survival_meta.csv")[row_type=="RE_pooled"]
R8 <- merge(S8[, .(miRNA, k_S8=k, HR_S8=HR, lo_S8=lo, hi_S8=hi, p_S8=p, I2_S8=I2, q_S8=q_BH)],
            P_UNADJ[, .(miRNA, k, HR=HR_RE, lo=lo_RE, hi=hi_RE, p=p_RE, I2, q=q_BH_RE)], by="miRNA", all=TRUE)
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-300)
R8[, identical_within_1e9 := k==k_S8 & rel(HR,HR_S8)<1e-9 & rel(lo,lo_S8)<1e-9 & rel(hi,hi_S8)<1e-9 & rel(p,p_S8)<1e-9 &
                             abs(I2-I2_S8)<1e-6 & rel(q,q_S8)<1e-9]
fwrite(R8, file.path(HERE, "table_s8_reproduction.csv"))
cat("Table S8 reproduced (11 miRNAs, relative tolerance 1e-9):", all(R8$identical_within_1e9), "\n")
stopifnot(nrow(R8) == 11, all(R8$identical_within_1e9))
print(P_UNADJ[, .(miRNA, k, HR_RE=signif(HR_RE,3), lo=signif(lo_RE,3), hi=signif(hi_RE,3), p=signif(p_RE,3), I2=signif(I2,3))])
print(P_ADJ[, .(miRNA, k, HR_RE=signif(HR_RE,3), lo=signif(lo_RE,3), hi=signif(hi_RE,3), p=signif(p_RE,3), q=signif(q_BH_RE,3), HR_FE=signif(HR_FE,3), p_FE=signif(p_FE,3), I2=signif(I2,3))])
