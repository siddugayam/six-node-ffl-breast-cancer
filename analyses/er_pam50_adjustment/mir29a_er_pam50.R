## mir29a_er_pam50.R -- hsa-miR-29a and the other ten miRNAs of Table S8: Cox models for each cohort's primary endpoint,
## unadjusted and adjusted for oestrogen-receptor status and/or PAM50 subtype, in the five miRNA cohorts that have ER or
## subtype annotation (TCGA-BRCA, GSE19783, GSE22216, GSE59829, GSE78870), both on the complete cases of the adjusted
## model; fixed- and random-effects pools with I2 and BH over the 11 miRNAs; Schoenfeld tests (see SETTINGS.md).
## Cohort loading, pick() and z-scoring are those of mirna_multivariable_cox.R (multivariable survival analysis), which
## copies them from c11_mirna_cox_meta.R.
## Root of the repository: environment variable FFL_ROOT. No sample-level values are written.
## usage: FFL_ROOT=<repository checkout> Rscript mir29a_er_pam50.R
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
add("TCGA_BRCA", CO$TCGA_BRCA$M, p, c("OS_time","OS_event"), covars=c("ER_ihc","PAM50"), covar_note="ER (IHC); PAM50")

## GSE19783
p <- CO$GSE19783$pheno
g <- readRDS("cache/external2/gse19783.rds")
stopifnot(identical(colnames(g$GSE19783_mRNA), colnames(CO$GSE19783$M)), identical(p$sample, colnames(CO$GSE19783$M)))
sc <- score_of(g$GSE19783_mRNA, cafA_genes); n_caf_19783 <- sc$n
p$stromal <- as.numeric(sc$v)
p$ER2 <- fac(p$ER, c("Positive","Negative"))
p$subtype5 <- fac(p$subtype, c("Lum A","Lum B","ERBB2","Basal","Normal-like"))
add("GSE19783", CO$GSE19783$M, p, c("OS_time","OS_event"), covars=c("ER2","subtype5"), covar_note="ER; PAM50 subtype")

## GSE22216 (matched mRNA: GSE22219)
p <- CO$GSE22216$pheno
cb <- readRDS("cache/external_cohorts/cohorts_base.rds")
lk <- read.csv("cache/metabolism/multicohort/GSE22219_GSE22216_link.csv", stringsAsFactors=FALSE)
sc <- score_of(cb$GSE22219$X, cafA_genes); n_caf_22219 <- sc$n
p$stromal <- as.numeric(sc$v[lk$gsm_mrna[match(p$gsm, lk$gsm_mirna)]])
p$grade3 <- factor(p$grade, levels=c(1,2,3))
p$age <- as.numeric(p$age); p$size <- as.numeric(p$size); p$node <- as.numeric(p$node); p$ER <- as.numeric(p$ER)
add("GSE22216", CO$GSE22216$M, p, c("DRFS_time","DRFS_event"), covars=c("ER"), covar_note="ER")

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
add("GSE59829", E, ph, c("DMFS_time","DMFS_event"), covars=c("ER2"), covar_note="ER")

## GSE78870
E <- as.matrix(fread(file.path(B,"GSE78870_expr.csv")), rownames=1)
P <- fread(file.path(B,"GSE78870_pheno.csv")); stopifnot(identical(colnames(E), P$gsm))
ph <- data.frame(sample=P$gsm,
                 TTP_time=suppressWarnings(as.numeric(P[["time-to-progression (in months)"]])),
                 TTP_event=suppressWarnings(as.numeric(P[["disease progression (event)"]])),
                 age=suppressWarnings(as.numeric(P[["age at start therapy"]])),
                 grade3=factor(P[["grade"]], levels=c(1,2,3)),
                 ER2=fac(P[["estrogen receptor status"]], c("positive","negative")), stringsAsFactors=FALSE)
add("GSE78870", E, ph, c("TTP_time","TTP_event"), covars=c("ER2"), covar_note="ER")


## ---------------------------------------------------------------- Cox fits --
fit1 <- function(co, z, covars, model){
  ph <- co$ph; endp <- co$endpoint
  d <- data.frame(z=z, time=suppressWarnings(as.numeric(ph[[endp[1]]])), ev=suppressWarnings(as.numeric(ph[[endp[2]]])))
  for(cv in covars) d[[cv]] <- ph[[cv]]
  d <- d[complete.cases(d) & is.finite(d$time) & d$time > 0, , drop=FALSE]
  for(cv in covars) if(is.factor(d[[cv]])) d[[cv]] <- droplevels(d[[cv]])
  if(nrow(d) < 25 || sum(d$ev) < 8) return(NULL)
  vars <- if(model == "adjusted") c("z", covars) else "z"
  fit <- try(coxph(as.formula(paste("Surv(time, ev) ~", paste(vars, collapse=" + "))), data=d), silent=TRUE)
  if(inherits(fit, "try-error")) return(NULL)
  s <- summary(fit)$coefficients
  zp <- try(cox.zph(fit, transform="km"), silent=TRUE)
  data.frame(n=nrow(d), nevent=sum(d$ev), logHR=s["z","coef"], se=s["z","se(coef)"], HR=exp(s["z","coef"]),
             lo=exp(s["z","coef"]-1.96*s["z","se(coef)"]), hi=exp(s["z","coef"]+1.96*s["z","se(coef)"]), p=s["z","Pr(>|z|)"],
             schoenfeld_global_p=if(inherits(zp,"try-error")) NA else zp$table["GLOBAL","p"],
             schoenfeld_mirna_p=if(inherits(zp,"try-error")) NA else zp$table["z","p"])
}
rows <- list()
for(cn in names(COH)){
  co <- COH[[cn]]; M <- co$M; stopifnot(ncol(M)==nrow(co$ph))
  for(tg in TARGETS){
    pk <- pick(M, tg); if(is.null(pk)) next
    z <- as.numeric(scale(pk$row))
    for(model in c("unadjusted_same_cases","adjusted")){
      r <- fit1(co, z, co$covars, if(model=="adjusted") "adjusted" else "unadjusted"); if(is.null(r)) next
      rows[[length(rows)+1]] <- cbind(data.frame(cohort=cn, miRNA=tg, probe=pk$name, endpoint=sub("_time$","",co$endpoint[1]),
        model=model, covariates=if(model=="adjusted") co$covar_note else "none (same cases as the adjusted model)",
        stringsAsFactors=FALSE), r)
    }
  }
}
COX <- rbindlist(rows)
fwrite(COX, file.path(HERE, "mir29a_er_pam50_percohort.csv"))
pool <- function(df, label){
  out <- list()
  for(tg in TARGETS){
    s <- df[miRNA==tg]; if(nrow(s) < 2) next
    re <- rma(yi=s$logHR, sei=s$se, method="DL"); fe <- rma(yi=s$logHR, sei=s$se, method="FE")
    out[[length(out)+1]] <- data.frame(analysis=label, miRNA=tg, k=nrow(s), cohorts=paste(s$cohort, collapse=";"),
      n=sum(s$n), nevent=sum(s$nevent), HR_FE=exp(fe$b[1]), lo_FE=exp(fe$ci.lb), hi_FE=exp(fe$ci.ub), p_FE=fe$pval,
      HR_RE=exp(re$b[1]), lo_RE=exp(re$ci.lb), hi_RE=exp(re$ci.ub), p_RE=re$pval, I2=re$I2, tau2=re$tau2, Q=re$QE, p_Q=re$QEp,
      stringsAsFactors=FALSE)
  }
  o <- rbindlist(out); o[, q_BH_FE := p.adjust(p_FE, "BH")]; o[, q_BH_RE := p.adjust(p_RE, "BH")]; o
}
P <- rbind(pool(COX[model=="unadjusted_same_cases"], "unadjusted, same cases"), pool(COX[model=="adjusted"], "adjusted for ER and/or PAM50"))
fwrite(P, file.path(HERE, "mir29a_er_pam50_pooled.csv"))
print(P[miRNA=="miR-29a"]); print(COX[miRNA=="miR-29a", .(cohort, model, covariates, n, nevent, HR, lo, hi, p, schoenfeld_mirna_p)])
