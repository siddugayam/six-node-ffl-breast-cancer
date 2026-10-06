## ==========================================================================
## c11_mirna_cox_meta.R  -- Part A
## Cox models for the 10 prioritised miRNAs (+ miR-29a control)
## in every breast-cancer miRNA cohort with an outcome, then DerSimonian-Laird
## random-effects meta-analysis per miRNA.
## ==========================================================================
suppressPackageStartupMessages({library(survival); library(metafor); library(data.table)})
setwd("/path/to/revision")
OUT <- "results/v6"; B <- "cache/v6/cohorts/built"
set.seed(1)

TARGETS <- c("miR-21","miR-195","miR-204","miR-383","miR-124","miR-155","miR-429",
             "miR-141","miR-34a","miR-101","miR-29a")
RX <- c("miR-21"="^hsa-mir-21(-[0-9])?(-[35]p)?$",
        "miR-195"="^hsa-mir-195(-[35]p)?$", "miR-204"="^hsa-mir-204(-[35]p)?$",
        "miR-383"="^hsa-mir-383(-[35]p)?$", "miR-124"="^hsa-mir-124a?(-[123])?(-[35]p)?$",
        "miR-155"="^hsa-mir-155(-[35]p)?$", "miR-429"="^hsa-mir-429(-[35]p)?$",
        "miR-141"="^hsa-mir-141(-[35]p)?$", "miR-34a"="^hsa-mir-34a(-[35]p)?$",
        "miR-101"="^hsa-mir-101(-[12])?(-[35]p)?$",
        "miR-29a"="^hsa-mir-29a(-[35]p)?$")

## pick the dominant (highest mean) non-star probe/row for a target miRNA
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

## ---------------------------------------------------------------- cohorts --
CO <- readRDS("cache/v4/multicohort/mirna_cohorts.rds")
COH <- list()
add <- function(name, accession, platform, M, ph, endpoints, strata=NULL, covars=NULL, note=""){
  COH[[name]] <<- list(name=name, accession=accession, platform=platform, M=M, ph=ph,
                       endpoints=endpoints, strata=strata, covars=covars, note=note)
}
add("TCGA_BRCA","TCGA-BRCA","Illumina HiSeq miRNA-seq (log2 RPM)",
    CO$TCGA_BRCA$M, CO$TCGA_BRCA$pheno,
    list(OS=c("OS_time","OS_event"), PFI=c("PFI_time","PFI_event"), DSS=c("DSS_time","DSS_event")),
    covars="age", note="reference cohort (discovery)")
add("GSE19783","GSE19783","Agilent-019118 human miRNA 2.0 (GPL8227)",
    CO$GSE19783$M, CO$GSE19783$pheno,
    list(OS=c("OS_time","OS_event"), BCSS=c("BCSS_time","BCSS_event")),
    note="Oslo MicMa; time variable is the deposited disease-free-survival time")
add("GSE22216","GSE22216 (= GSE22220 miRNA arm)","Illumina Human v1 MicroRNA expression beadchip (GPL8178)",
    CO$GSE22216$M, CO$GSE22216$pheno, list(DRFS=c("DRFS_time","DRFS_event")),
    covars=c("age","grade","size","node"), note="Buffa et al., 210 early primary BC")
p <- CO$GSE37405$pheno
add("GSE37405","GSE37405","Exiqon miRCURY LNA (GPL13703/14149/15462)",
    CO$GSE37405$M, p, list(RFS=c("RFS_time","RFS_event"), OS=c("OS_time","OS_event")),
    strata="platform", covars=c("age","size","node"),
    note="high-risk ER+ on adjuvant tamoxifen (DBCG); 3 array generations, platform as stratum")

## ---- NEW cohort 1: GSE59829 ----------------------------------------------
E <- as.matrix(fread(file.path(B,"GSE59829_expr.csv")), rownames=1)
P <- fread(file.path(B,"GSE59829_pheno.csv"))
stopifnot(identical(colnames(E), P$gsm))
ph <- data.frame(sample=P$gsm,
                 DMFS_time=suppressWarnings(as.numeric(P[["dfs (months)"]])),
                 DMFS_event=suppressWarnings(as.numeric(P[["distant_metastasis"]])),
                 age=suppressWarnings(as.numeric(P[["age"]])),
                 size=suppressWarnings(as.numeric(P[["tumor_size"]])),
                 ER=P[["esr1_status"]], stringsAsFactors=FALSE)
add("GSE59829","GSE59829","Illumina Human v2 MicroRNA expression beadchip (GPL8179)", E, ph,
    list(DMFS=c("DMFS_time","DMFS_event")), covars=c("age","size"),
    note="NEW in v6: 123 primary BC, distant-metastasis-free survival (Ivan/D'Assoro et al.)")

## ---- NEW cohort 2: GSE78870 ----------------------------------------------
E <- as.matrix(fread(file.path(B,"GSE78870_expr.csv")), rownames=1)
P <- fread(file.path(B,"GSE78870_pheno.csv"))
stopifnot(identical(colnames(E), P$gsm))
ph <- data.frame(sample=P$gsm,
                 TTP_time=suppressWarnings(as.numeric(P[["time-to-progression (in months)"]])),
                 TTP_event=suppressWarnings(as.numeric(P[["disease progression (event)"]])),
                 DFI_time=suppressWarnings(as.numeric(P[["disease-free interval (in months)"]])),
                 DFI_event=suppressWarnings(as.numeric(P[["disease relapse (event)"]])),
                 age=suppressWarnings(as.numeric(P[["age at start therapy"]])), stringsAsFactors=FALSE)
add("GSE78870","GSE78870","TaqMan microRNA Low-Density Array pools A and B v2.0 (GPL20662)", E, ph,
    list(TTP=c("TTP_time","TTP_event"), DFI=c("DFI_time","DFI_event")), covars="age",
    note="NEW in v6: 106 ER+ BC on first-line aromatase inhibitor; time to progression")

PRIMARY <- c(TCGA_BRCA="OS", GSE19783="OS", GSE22216="DRFS", GSE37405="RFS",
             GSE59829="DMFS", GSE78870="TTP")

## ---------------------------------------------------------------- Cox ------
rows <- list(); probes <- list()
for(cn in names(COH)){
  co <- COH[[cn]]; M <- co$M; ph <- co$ph
  stopifnot(ncol(M)==nrow(ph))
  for(tg in TARGETS){
    pk <- pick(M, tg)
    if(is.null(pk)) next
    x <- pk$row
    ## z-score within cohort (within platform stratum if declared)
    z <- rep(NA_real_, length(x))
    if(!is.null(co$strata)){
      s <- ph[[co$strata]]
      for(lv in unique(s)){ i <- which(s==lv); z[i] <- as.numeric(scale(x[i])) }
    } else z <- as.numeric(scale(x))
    probes[[length(probes)+1]] <- data.frame(cohort=cn, accession=co$accession,
        platform=co$platform, miRNA=tg, probe=pk$name, n_probes_matched=pk$n_matched,
        all_matched=pk$all, stringsAsFactors=FALSE)
    for(ep in names(co$endpoints)){
      tv <- co$endpoints[[ep]][1]; ev <- co$endpoints[[ep]][2]
      d <- data.frame(z=z, time=suppressWarnings(as.numeric(ph[[tv]])),
                      ev=suppressWarnings(as.numeric(ph[[ev]])))
      if(!is.null(co$strata)) d$stratum <- ph[[co$strata]]
      for(cv in co$covars) d[[cv]] <- suppressWarnings(as.numeric(ph[[cv]]))
      for(model in c("univariate","adjusted")){
        if(model=="adjusted" && is.null(co$covars)) next
        vars <- if(model=="univariate") "z" else c("z", co$covars)
        dd <- d[, c("z","time","ev", if(!is.null(co$strata)) "stratum",
                    if(model=="adjusted") co$covars), drop=FALSE]
        dd <- dd[complete.cases(dd) & is.finite(dd$time) & dd$time>0, , drop=FALSE]
        if(nrow(dd)<25 || sum(dd$ev)<8) next
        f <- paste0("Surv(time, ev) ~ ", paste(vars, collapse=" + "),
                    if(!is.null(co$strata)) " + strata(stratum)" else "")
        fit <- try(coxph(as.formula(f), data=dd), silent=TRUE)
        if(inherits(fit,"try-error")) next
        s <- summary(fit)$coefficients
        rows[[length(rows)+1]] <- data.frame(cohort=cn, accession=co$accession,
          platform=co$platform, miRNA=tg, probe=pk$name, endpoint=ep, model=model,
          n=nrow(dd), nevent=sum(dd$ev), logHR=s["z","coef"], se=s["z","se(coef)"],
          HR=exp(s["z","coef"]), lo=exp(s["z","coef"]-1.96*s["z","se(coef)"]),
          hi=exp(s["z","coef"]+1.96*s["z","se(coef)"]), p=s["z","Pr(>|z|)"],
          primary=as.integer(ep==PRIMARY[[cn]]), stringsAsFactors=FALSE)
      }
    }
  }
}
COX <- rbindlist(rows); PROBE <- rbindlist(probes)
fwrite(COX, file.path(OUT,"mirna_meta_percohort_cox.csv"))
fwrite(PROBE, file.path(OUT,"cohorts_mirna_probe_mapping.csv"))
cat("Cox fits:", nrow(COX), " cohorts:", length(unique(COX$cohort)), "\n")

## ------------------------------------------------------- meta-analysis -----
pool <- function(df, label){
  out <- list()
  for(tg in unique(df$miRNA)){
    s <- df[miRNA==tg]
    if(nrow(s)<2) next
    m <- try(rma(yi=s$logHR, sei=s$se, method="DL"), silent=TRUE)
    if(inherits(m,"try-error")) next
    w <- weights(m)
    for(i in seq_len(nrow(s)))
      out[[length(out)+1]] <- data.frame(analysis=label, miRNA=tg, row_type="study",
        cohort=s$cohort[i], accession=s$accession[i], endpoint=s$endpoint[i],
        n=s$n[i], nevent=s$nevent[i], HR=s$HR[i], lo=s$lo[i], hi=s$hi[i], p=s$p[i],
        weight_pct=w[i], k=nrow(s), I2=NA, tau2=NA, Q=NA, p_Q=NA, stringsAsFactors=FALSE)
    out[[length(out)+1]] <- data.frame(analysis=label, miRNA=tg, row_type="RE_pooled",
      cohort="POOLED (DerSimonian-Laird)", accession="", endpoint="primary per cohort",
      n=sum(s$n), nevent=sum(s$nevent), HR=exp(m$b[1]), lo=exp(m$ci.lb), hi=exp(m$ci.ub),
      p=m$pval, weight_pct=100, k=m$k, I2=m$I2, tau2=m$tau2, Q=m$QE, p_Q=m$QEp,
      stringsAsFactors=FALSE)
  }
  rbindlist(out)
}
PRIM <- COX[model=="univariate" & primary==1]
M1 <- pool(PRIM, "primary_endpoint_univariate")
M2 <- pool(PRIM[cohort!="TCGA_BRCA"], "primary_endpoint_univariate_noTCGA")
M3 <- pool(COX[model=="adjusted" & primary==1], "primary_endpoint_adjusted")
META <- rbindlist(list(M1,M2,M3))
fwrite(META, file.path(OUT,"mirna_meta_pooled.csv"))

## forest-plot-ready table (primary univariate)
FOR <- META[analysis=="primary_endpoint_univariate"]
FOR[, `:=`(label=ifelse(row_type=="RE_pooled","POOLED (RE, DL)",paste0(cohort," (",endpoint,")")),
           ci=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))]
fwrite(FOR, file.path(OUT,"mirna_meta_forest_table.csv"))

cat("\n================ POOLED (primary endpoint, univariate, per SD) ==========\n")
pp <- META[analysis=="primary_endpoint_univariate" & row_type=="RE_pooled"][order(p)]
print(as.data.frame(pp[,.(miRNA,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
                          p=signif(p,3),I2=round(I2,1),p_Q=signif(p_Q,3))]))
cat("\n============ POOLED without TCGA (independent replication) ==============\n")
pp2 <- META[analysis=="primary_endpoint_univariate_noTCGA" & row_type=="RE_pooled"][order(p)]
print(as.data.frame(pp2[,.(miRNA,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
                          p=signif(p,3),I2=round(I2,1))]))
cat("\n============ POOLED adjusted ============\n")
pp3 <- META[analysis=="primary_endpoint_adjusted" & row_type=="RE_pooled"][order(p)]
print(as.data.frame(pp3[,.(miRNA,k,n,nevent,HR=round(HR,3),lo=round(lo,3),hi=round(hi,3),
                          p=signif(p,3),I2=round(I2,1))]))
saveRDS(COH, "cache/v6/cohorts/built/mirna_cohorts_v6.rds")
