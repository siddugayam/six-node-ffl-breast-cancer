#!/usr/bin/env Rscript
## 45_farmer_survival.R
## Task E(i): do the FFL module score and the Farmer stromal score predict survival in
## TCGA-BRCA and METABRIC, and does the FFL module carry any prognostic information that
## the published stromal programme does not?  In METABRIC we also test the *interaction*
## with chemotherapy, because Farmer's claim is that the signature is predictive of
## chemoresistance rather than prognostic.
suppressPackageStartupMessages({
  library(data.table); library(survival); library(GSVA); library(matrixStats)
})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v5/45_farmer_survival.log")
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); set.seed(11)

sets <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
KEY <- c("FARMER_STROMAL","FFL_3node","FFL_higher_order_only","FFL_6node_all",
         "MS_exemplar_module","ESTIMATE_STROMAL_noCOL","CAF_scRNA_50","WEST_DTF_FIBROMATOSIS",
         "WINSLOW_STROMAL_SIG1","FINAK_SDPP","HALLMARK_EMT","HALLMARK_TGF_BETA",
         "CHANG_WOUND_UP_VANTVEER","NABA_CORE_MATRISOME","TRIULZI_ECM")

cox_one <- function(tt, ee, x, covars = NULL, label = "") {
  ok <- is.finite(tt) & tt > 0 & is.finite(ee) & is.finite(x)
  d <- data.frame(time = tt[ok], ev = ee[ok], sc = as.numeric(scale(x[ok])))
  if (!is.null(covars)) d <- cbind(d, covars[ok, , drop = FALSE])
  d <- d[complete.cases(d), ]
  if (nrow(d) < 50 || sum(d$ev == 1) < 10) return(NULL)
  f <- if (is.null(covars)) coxph(Surv(time, ev) ~ sc, data = d)
       else coxph(as.formula(paste("Surv(time, ev) ~ sc +",
                                   paste(colnames(covars), collapse = " + "))), data = d)
  s <- summary(f)
  data.table(model = label, n = nrow(d), n_event = sum(d$ev == 1),
             HR_per_SD = s$coefficients["sc","exp(coef)"],
             CI_low = s$conf.int["sc","lower .95"], CI_high = s$conf.int["sc","upper .95"],
             p = s$coefficients["sc","Pr(>|z|)"], C_index = s$concordance[1])
}

## ==================================================================== TCGA-BRCA =========
sc   <- readRDS(file.path(BASE, "cache/v5/farmer_tcga_scores.rds"))
surv <- as.data.table(readRDS(file.path(BASE, "data/brca_survival_clean.rds")))
S <- sc$samples
cl <- surv[match(S, sample_id)]
logf("TCGA samples with a score: ", length(S), " | with survival row: ", sum(!is.na(cl$sample_id)))
SS <- sc$ssgsea
rows <- list(); k <- 0L
for (nm in intersect(KEY, rownames(SS))) {
  x <- as.numeric(SS[nm, ])
  for (ep in c("OS","PFI","DSS")) {
    tt <- cl[[paste0(ep, ".time")]]; ee <- cl[[ep]]
    r1 <- cox_one(tt, ee, x, NULL, "univariable")
    cv <- data.frame(age = cl$age, stage = factor(cl$stage_group))
    r2 <- cox_one(tt, ee, x, cv, "adjusted_age_stage")
    for (r in list(r1, r2)) if (!is.null(r)) { k <- k + 1L
      rows[[k]] <- cbind(data.table(cohort = "TCGA-BRCA", endpoint = ep, score = nm), r) }
  }
}
## head to head: does the FFL module add over the Farmer stromal score?
h2h <- list(); h <- 0L
for (mod in c("FFL_3node","FFL_higher_order_only","FFL_6node_all","MS_exemplar_module")) {
  for (ep in c("OS","PFI","DSS")) {
    tt <- cl[[paste0(ep,".time")]]; ee <- cl[[ep]]
    ok <- is.finite(tt) & tt > 0 & is.finite(ee)
    d <- data.frame(time = tt[ok], ev = ee[ok],
                    m = as.numeric(scale(SS[mod, ok])), f = as.numeric(scale(SS["FARMER_STROMAL", ok])))
    d <- d[complete.cases(d), ]
    if (nrow(d) < 50 || sum(d$ev == 1) < 10) next
    Sv <- Surv(d$time, d$ev)
    mm <- coxph(Sv ~ m, d); mf <- coxph(Sv ~ f, d); mb <- coxph(Sv ~ m + f, d)
    a1 <- anova(mf, mb); a2 <- anova(mm, mb)
    h <- h + 1L
    h2h[[h]] <- data.table(cohort = "TCGA-BRCA", endpoint = ep, module = mod,
      n = nrow(d), n_event = sum(d$ev == 1),
      HR_module_alone = summary(mm)$coefficients[1,"exp(coef)"], p_module_alone = summary(mm)$coefficients[1,"Pr(>|z|)"],
      HR_farmer_alone = summary(mf)$coefficients[1,"exp(coef)"], p_farmer_alone = summary(mf)$coefficients[1,"Pr(>|z|)"],
      HR_module_adj_farmer = summary(mb)$coefficients["m","exp(coef)"], p_module_adj_farmer = summary(mb)$coefficients["m","Pr(>|z|)"],
      HR_farmer_adj_module = summary(mb)$coefficients["f","exp(coef)"], p_farmer_adj_module = summary(mb)$coefficients["f","Pr(>|z|)"],
      LRT_p_module_adds_over_farmer = a1[2,"Pr(>|Chi|)"],
      LRT_p_farmer_adds_over_module = a2[2,"Pr(>|Chi|)"],
      C_module = concordance(mm)$concordance, C_farmer = concordance(mf)$concordance,
      C_both = concordance(mb)$concordance)
  }
}

## ==================================================================== METABRIC ==========
mb <- readRDS(file.path(BASE, "data/metabric.rds"))
M <- mb$M; CL <- as.data.table(mb$cl)
logf("METABRIC expression: ", paste(dim(M), collapse = " x "))
M <- M[rowSums(is.na(M)) == 0 & rowSds(M) > 0, , drop = FALSE]
gs_mb <- lapply(sets, function(s) intersect(s, rownames(M)))
gs_mb <- gs_mb[sapply(gs_mb, length) >= 2]
logf("METABRIC: signatures scorable = ", length(gs_mb),
     " | FARMER genes present = ", length(gs_mb$FARMER_STROMAL), "/", length(sets$FARMER_STROMAL))
SM <- gsva(ssgseaParam(M, gs_mb, minSize = 2, normalize = TRUE), verbose = FALSE)
saveRDS(SM, file.path(BASE, "cache/v5/farmer_metabric_scores.rds"))
CL <- CL[match(colnames(M), SAMPLE_ID)]
stopifnot(identical(CL$SAMPLE_ID, colnames(M)))
CL[, os_ev := as.integer(grepl("^1", OS_STATUS))]
CL[, rfs_ev := as.integer(grepl("^1", RFS_STATUS))]
CL[, chemo := ifelse(CHEMOTHERAPY == "YES", 1L, ifelse(CHEMOTHERAPY == "NO", 0L, NA_integer_))]
logf("METABRIC: OS events ", sum(CL$os_ev, na.rm=TRUE), " | RFS events ", sum(CL$rfs_ev, na.rm=TRUE),
     " | chemotherapy YES ", sum(CL$chemo == 1, na.rm=TRUE), " NO ", sum(CL$chemo == 0, na.rm=TRUE))
for (nm in intersect(KEY, rownames(SM))) {
  x <- as.numeric(SM[nm, ])
  for (ep in c("OS","RFS")) {
    tt <- if (ep == "OS") CL$OS_MONTHS else CL$RFS_MONTHS
    ee <- if (ep == "OS") CL$os_ev else CL$rfs_ev
    r1 <- cox_one(tt, ee, x, NULL, "univariable")
    cv <- data.frame(age = CL$AGE_AT_DIAGNOSIS, grade = CL$GRADE,
                     er = factor(CL$ER_STATUS), nodes = CL$LYMPH_NODES_EXAMINED_POSITIVE)
    r2 <- cox_one(tt, ee, x, cv, "adjusted_age_grade_ER_nodes")
    for (r in list(r1, r2)) if (!is.null(r)) { k <- k + 1L
      rows[[k]] <- cbind(data.table(cohort = "METABRIC", endpoint = ep, score = nm), r) }
  }
}
allc <- rbindlist(rows)
allc[, q_BH := p.adjust(p, "BH"), by = .(cohort, endpoint, model)]
fwrite(allc, file.path(BASE, "results/v5/farmer_survival_cox.csv"))
logf("WROTE results/v5/farmer_survival_cox.csv rows=", nrow(allc))
logf("\n=== univariable Cox, HR per SD ===")
for (co in unique(allc$cohort)) for (ep in unique(allc[cohort==co, endpoint])) {
  s <- allc[cohort==co & endpoint==ep & model=="univariable"][order(p)]
  for (i in seq_len(nrow(s))) logf(sprintf("  %-10s %-4s %-26s HR=%.3f (%.3f-%.3f) p=%.3g q=%.3g C=%.3f",
    co, ep, s$score[i], s$HR_per_SD[i], s$CI_low[i], s$CI_high[i], s$p[i], s$q_BH[i], s$C_index[i]))
}

## ---------------- METABRIC: chemotherapy interaction (predictive vs prognostic) ---------
logf("\n=== METABRIC: score x chemotherapy interaction (Farmer's claim is *predictive*) ===")
int_rows <- list(); ii <- 0L
for (nm in intersect(KEY, rownames(SM))) {
  x <- as.numeric(scale(SM[nm, ]))
  for (ep in c("OS","RFS")) {
    tt <- if (ep=="OS") CL$OS_MONTHS else CL$RFS_MONTHS
    ee <- if (ep=="OS") CL$os_ev else CL$rfs_ev
    d <- data.frame(time = tt, ev = ee, sc = x, chemo = CL$chemo,
                    age = CL$AGE_AT_DIAGNOSIS, grade = CL$GRADE, er = factor(CL$ER_STATUS),
                    nodes = CL$LYMPH_NODES_EXAMINED_POSITIVE)
    d <- d[complete.cases(d) & is.finite(d$time) & d$time > 0, ]
    if (nrow(d) < 100) next
    f0 <- coxph(Surv(time, ev) ~ sc + chemo + age + grade + er + nodes, data = d)
    f1 <- coxph(Surv(time, ev) ~ sc * chemo + age + grade + er + nodes, data = d)
    lr <- anova(f0, f1)
    dch <- d[d$chemo == 1, ]; dno <- d[d$chemo == 0, ]
    hc <- coxph(Surv(time, ev) ~ sc + age + grade + er + nodes, data = dch)
    hn <- coxph(Surv(time, ev) ~ sc + age + grade + er + nodes, data = dno)
    ii <- ii + 1L
    int_rows[[ii]] <- data.table(score = nm, endpoint = ep, n = nrow(d),
      n_chemo = nrow(dch), n_nochemo = nrow(dno),
      HR_in_chemo = summary(hc)$coefficients["sc","exp(coef)"],
      p_in_chemo  = summary(hc)$coefficients["sc","Pr(>|z|)"],
      HR_in_nochemo = summary(hn)$coefficients["sc","exp(coef)"],
      p_in_nochemo  = summary(hn)$coefficients["sc","Pr(>|z|)"],
      interaction_coef = summary(f1)$coefficients["sc:chemo","coef"],
      interaction_p = summary(f1)$coefficients["sc:chemo","Pr(>|z|)"],
      LRT_p_interaction = lr[2,"Pr(>|Chi|)"])
  }
}
itab <- rbindlist(int_rows)
itab[, q_BH := p.adjust(interaction_p, "BH"), by = endpoint]
fwrite(itab, file.path(BASE, "results/v5/farmer_metabric_chemo_interaction.csv"))
for (i in seq_len(nrow(itab))) { r <- itab[i]
  logf(sprintf("  %-26s %-4s chemo n=%d HR=%.3f p=%.3g | no-chemo n=%d HR=%.3f p=%.3g | interaction p=%.3g",
    r$score, r$endpoint, r$n_chemo, r$HR_in_chemo, r$p_in_chemo,
    r$n_nochemo, r$HR_in_nochemo, r$p_in_nochemo, r$interaction_p)) }

h2t <- rbindlist(h2h)
fwrite(h2t, file.path(BASE, "results/v5/farmer_vs_ffl_module_headtohead_survival.csv"))
logf("\n=== TCGA head-to-head: does the FFL module add anything over the Farmer stromal score? ===")
for (i in seq_len(nrow(h2t))) { r <- h2t[i]
  logf(sprintf("  %-4s %-22s module alone HR=%.3f p=%.3g | Farmer alone HR=%.3f p=%.3g | module|Farmer HR=%.3f p=%.3g | LRT module adds p=%.3g",
    r$endpoint, r$module, r$HR_module_alone, r$p_module_alone, r$HR_farmer_alone, r$p_farmer_alone,
    r$HR_module_adj_farmer, r$p_module_adj_farmer, r$LRT_p_module_adds_over_farmer)) }
logf("DONE 45")
