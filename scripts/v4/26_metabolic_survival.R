## 26_metabolic_survival.R -------------------------------------------------
## Part E: Cox models for every metabolic score (OS, PFI, DSS), per SD,
## adjusted for age at diagnosis and AJCC pathologic stage. BH-FDR per endpoint.
suppressPackageStartupMessages({library(data.table); library(survival)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
ss <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
gv <- readRDS(file.path(RES,"metabolic_gsva_scores.rds"))
ph <- readRDS(file.path(ROOT,"data/brca_pheno.rds"))
sv <- as.data.table(readRDS(file.path(ROOT,"data/brca_survival.rds")))
inv <- fread(file.path(RES,"metabolic_geneset_inventory.csv"))

tum <- ph$sample[ph$sample_type=="Primary Tumor"]
sv <- sv[sample %in% intersect(tum, colnames(ss))]
sv[, stage_raw := `ajcc_pathologic_tumor_stage`]
sv[, stage := fcase(grepl("Stage IV", stage_raw), "IV",
                    grepl("Stage III", stage_raw), "III",
                    grepl("Stage II", stage_raw), "II",
                    grepl("Stage I($| |A|B|C)", stage_raw), "I",
                    default = NA_character_)]
sv[, stage := factor(stage, levels=c("I","II","III","IV"))]
sv[, age := age_at_initial_pathologic_diagnosis]
cat("tumours with score+survival:", nrow(sv), "\n")
print(table(sv$stage, useNA="ifany"))
cat("missing age:", sum(is.na(sv$age)), "\n")
keep <- sv[!is.na(stage) & !is.na(age)]
cat("analysis set (stage + age complete):", nrow(keep), "\n")
for (ep in c("OS","PFI","DSS"))
  cat(ep, "events:", sum(keep[[ep]]==1, na.rm=TRUE), " n with time:", sum(!is.na(keep[[paste0(ep,".time")]])), "\n")

run_cox <- function(M, tag) {
  s <- keep$sample
  Z <- t(scale(t(M[, s, drop=FALSE])))      # per-SD across the analysed tumours
  rbindlist(lapply(c("OS","PFI","DSS"), function(ep) {
    ev <- keep[[ep]]; ti <- keep[[paste0(ep,".time")]]
    ok <- !is.na(ev) & !is.na(ti) & ti > 0
    y <- Surv(ti[ok], ev[ok]); age <- keep$age[ok]; stg <- keep$stage[ok]
    rbindlist(lapply(rownames(Z), function(p) {
      x <- Z[p, ok]
      f <- tryCatch(coxph(y ~ x + age + stg), error=function(e) NULL)
      if (is.null(f)) return(NULL)
      cf <- summary(f)$coefficients
      f0 <- tryCatch(coxph(y ~ x), error=function(e) NULL)
      cf0 <- if (is.null(f0)) rep(NA_real_,5) else summary(f0)$coefficients["x",]
      zph <- tryCatch(cox.zph(f)$table["x","p"], error=function(e) NA_real_)
      data.table(set_id=p, endpoint=ep, score=tag, n=sum(ok), events=sum(ev[ok]==1),
                 HR_per_SD=exp(cf["x","coef"]), lo=exp(cf["x","coef"]-1.96*cf["x","se(coef)"]),
                 hi=exp(cf["x","coef"]+1.96*cf["x","se(coef)"]), z=cf["x","z"], p=cf["x","Pr(>|z|)"],
                 HR_unadjusted=exp(cf0[1]), p_unadjusted=cf0[5], p_ph_assumption=zph)
    }))
  }))
}
out <- rbind(run_cox(ss,"ssGSEA"), run_cox(gv,"GSVA"))
out[, FDR := p.adjust(p,"BH"), by=.(endpoint, score)]
out[, FDR_unadjusted := p.adjust(p_unadjusted,"BH"), by=.(endpoint, score)]
out <- merge(out, inv[,.(set_id, collection, n_genes_in_universe)], by="set_id", all.x=TRUE)
setorder(out, score, endpoint, p)
fwrite(out, file.path(RES,"metabolic_survival_cox.csv"))

cat("\n== significant at FDR<0.05, age+stage adjusted ==\n")
print(out[, .(tested=.N, sig_FDR05=sum(FDR<0.05), sig_raw_p05=sum(p<0.05)), by=.(score, endpoint)])
for (ep in c("OS","PFI","DSS")) {
  cat("\n---", ep, "(ssGSEA, age+stage adjusted) top 12 by p ---\n")
  print(out[score=="ssGSEA" & endpoint==ep][order(p)][1:12,
    .(set_id=substr(set_id,1,58), collection, HR=round(HR_per_SD,3),
      CI=paste0(round(lo,2),"-",round(hi,2)), p=signif(p,3), FDR=signif(FDR,3))])
}
## concordance between methods on the prognostic direction
w <- dcast(out[endpoint=="OS"], set_id ~ score, value.var="HR_per_SD")
cat("\nOS HR ssGSEA vs GSVA Pearson:", round(cor(log(w$ssGSEA), log(w$GSVA), use="complete.obs"),3), "\n")
## does adjustment matter?
cat("median |log HR| unadjusted vs adjusted (OS, ssGSEA):",
    round(median(abs(log(out[score=="ssGSEA"&endpoint=="OS", HR_unadjusted]))),4), "vs",
    round(median(abs(log(out[score=="ssGSEA"&endpoint=="OS", HR_per_SD]))),4), "\n")
cat("DONE 26\n")
