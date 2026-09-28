#!/usr/bin/env Rscript
# 15_external_metabric.R
# EXTERNAL COHORT 2: METABRIC (cBioPortal datahub study brca_metabric),
# Illumina HT-12 v3 microarray, n=1980 tumours with expression, OS and RFS.
# Replicates (a) the univariable hub Cox results and (b) the module-score survival
# association found (or not found) in TCGA-BRCA.
# NOTE METABRIC has NO miRNA assay, so every module is additionally recomputed in
# TCGA restricted to its protein-coding members, giving a like-for-like comparison.
# Output: results/external_METABRIC_cox.csv, results/external_METABRIC_modules.csv

suppressPackageStartupMessages({
  library(data.table); library(survival); library(GSVA)
})
BASE <- "/path/to/revision"
EXT  <- file.path(BASE, "cache/external")
LOG  <- file.path(BASE, "logs", "15_external_metabric.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")
set.seed(1234)

## ------------------------------------------------------------- expression ----
E <- fread(file.path(EXT, "metabric_data_mrna_illumina_microarray.txt"))
logf("METABRIC raw expression: ", nrow(E), " rows x ", ncol(E) - 2, " samples")
E <- E[Hugo_Symbol != "" & !is.na(Hugo_Symbol)]
scols <- setdiff(names(E), c("Hugo_Symbol", "Entrez_Gene_Id"))
M <- as.matrix(E[, ..scols]); rownames(M) <- E$Hugo_Symbol
ord <- order(apply(M, 1, function(z) stats::var(z, na.rm = TRUE)), decreasing = TRUE)
M <- M[ord, , drop = FALSE]; M <- M[!duplicated(rownames(M)), , drop = FALSE]
v <- apply(M, 1, function(z) stats::sd(z, na.rm = TRUE))
M <- M[is.finite(v) & v > 0, , drop = FALSE]
logf("METABRIC symbol-level matrix (dupes collapsed by max variance, zero-var dropped): ",
     paste(dim(M), collapse = " x "),
     " ; value range ", paste(round(range(M, na.rm = TRUE), 2), collapse = " to "))

## --------------------------------------------------------------- clinical ----
cp <- fread(file.path(EXT, "metabric_data_clinical_patient.txt"), skip = "PATIENT_ID")
cs <- fread(file.path(EXT, "metabric_data_clinical_sample.txt"), skip = "PATIENT_ID")
cl <- merge(cp, cs, by = "PATIENT_ID")
logf("clinical merged rows=", nrow(cl))
cl <- cl[SAMPLE_ID %in% colnames(M)]
logf("clinical rows with expression=", nrow(cl))

cl[, OS_time  := suppressWarnings(as.numeric(OS_MONTHS))]
cl[, OS_event := ifelse(grepl("^1", OS_STATUS), 1L, ifelse(grepl("^0", OS_STATUS), 0L, NA_integer_))]
cl[, RFS_time := suppressWarnings(as.numeric(RFS_MONTHS))]
cl[, RFS_event:= ifelse(grepl("^1", RFS_STATUS), 1L, ifelse(grepl("^0", RFS_STATUS), 0L, NA_integer_))]
cl[, age := suppressWarnings(as.numeric(AGE_AT_DIAGNOSIS))]
cl[, TUMOR_STAGE := as.character(TUMOR_STAGE)]
cl[, stage_group := factor(fcase(TUMOR_STAGE == "0", "0", TUMOR_STAGE == "1", "I",
                                 TUMOR_STAGE == "2", "II", TUMOR_STAGE == "3", "III",
                                 TUMOR_STAGE == "4", "IV", default = NA_character_),
                           levels = c("0", "I", "II", "III", "IV"))]
cl[, grade := factor(fifelse(as.character(GRADE) %in% c("1", "2", "3"), as.character(GRADE), NA_character_))]
logf("OS: n non-missing time=", sum(!is.na(cl$OS_time)), " events=", sum(cl$OS_event == 1, na.rm = TRUE))
logf("RFS: n non-missing time=", sum(!is.na(cl$RFS_time)), " events=", sum(cl$RFS_event == 1, na.rm = TRUE))
logf("stage: ", paste(sprintf("%s=%d", names(table(cl$stage_group, useNA = "ifany")),
                              table(cl$stage_group, useNA = "ifany")), collapse = " "))
logf("grade: ", paste(sprintf("%s=%d", names(table(cl$grade, useNA = "ifany")),
                              table(cl$grade, useNA = "ifany")), collapse = " "))

M <- M[, cl$SAMPLE_ID, drop = FALSE]
stopifnot(identical(colnames(M), cl$SAMPLE_ID))
saveRDS(list(M = M, cl = cl), file.path(BASE, "data/metabric.rds"))

## ==================== A. univariable Cox on network nodes ===================
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
prot  <- nodes[type != "miRNA", name]
tested <- intersect(prot, rownames(M))
logf("network protein-coding nodes present on the METABRIC array: ", length(tested),
     " / ", length(prot))
missing_nodes <- setdiff(prot, rownames(M))
logf("absent from METABRIC array (", length(missing_nodes), "): ",
     paste(head(missing_nodes, 40), collapse = ", "))

run_cox <- function(x, tt, ee) {
  ok <- is.finite(x) & is.finite(tt) & tt > 0 & is.finite(ee)
  if (sum(ok) < 50 || sum(ee[ok] == 1) < 10) return(NULL)
  xs <- as.numeric(scale(x[ok]))
  f <- try(coxph(Surv(tt[ok], ee[ok]) ~ xs), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  s <- summary(f)
  list(n = sum(ok), ev = sum(ee[ok] == 1), HR = s$coefficients[1, "exp(coef)"],
       lo = s$conf.int[1, "lower .95"], hi = s$conf.int[1, "upper .95"],
       p = s$coefficients[1, "Pr(>|z|)"], C = s$concordance[1])
}

rows <- list(); k <- 0L
for (gn in tested) {
  x <- as.numeric(M[gn, ])
  for (ep in c("OS", "RFS")) {
    tt <- cl[[paste0(ep, "_time")]]; ee <- cl[[paste0(ep, "_event")]]
    r <- run_cox(x, tt, ee); if (is.null(r)) next
    k <- k + 1L
    rows[[k]] <- data.table(feature = gn, endpoint = ep, n = r$n, n_event = r$ev,
                            HR_per_SD = r$HR, CI_low = r$lo, CI_high = r$hi,
                            p_value = r$p, C_index = r$C)
  }
}
mb <- rbindlist(rows)
mb[, q_value := p.adjust(p_value, method = "BH"), by = endpoint]
fwrite(mb, file.path(BASE, "results/external_METABRIC_cox.csv"))
logf("WROTE results/external_METABRIC_cox.csv rows=", nrow(mb))
for (ep in c("OS", "RFS"))
  logf(sprintf("METABRIC %s: tested=%d nominal p<0.05=%d BH q<0.05=%d",
               ep, nrow(mb[endpoint == ep]), sum(mb[endpoint == ep]$p_value < 0.05),
               sum(mb[endpoint == ep]$q_value < 0.05)))

## ==================== B. module scores in BOTH cohorts ======================
sets <- readRDS(file.path(BASE, "data/ffl_module_sets.rds"))
mods_all <- list(
  MS_exemplar_module    = sets$MS_MODULE,
  FFL_3node             = sets$N3,
  FFL_higher_order_only = sets$HIGHER_ONLY,
  FFL_6node_all         = sets$N6)
# protein-coding restriction (METABRIC has no miRNA assay)
mods_pc <- lapply(mods_all, function(s) intersect(s, prot))
for (nmm in names(mods_pc))
  logf(sprintf("module %-22s all=%d protein-coding=%d in-METABRIC=%d",
               nmm, length(mods_all[[nmm]]), length(mods_pc[[nmm]]),
               length(intersect(mods_pc[[nmm]], rownames(M)))))
mods_mb <- lapply(mods_pc, function(s) intersect(s, rownames(M)))
mods_mb <- mods_mb[sapply(mods_mb, length) >= 2]

logf("running GSVA on METABRIC ...")
gs_mb <- gsva(gsvaParam(M, mods_mb, kcdf = "Gaussian", minSize = 2), verbose = FALSE)
Zmb <- t(scale(t(M)))
mz_mb <- t(sapply(mods_mb, function(s) colMeans(Zmb[s, , drop = FALSE], na.rm = TRUE)))
logf("METABRIC module scores: GSVA ", paste(dim(gs_mb), collapse = "x"),
     " meanZ ", paste(dim(mz_mb), collapse = "x"))

analyse <- function(S, meth, cohort, clin, eps) {
  out <- list(); i <- 0L
  for (mn in rownames(S)) {
    scz <- as.numeric(scale(as.numeric(S[mn, ])))
    for (ep in names(eps)) {
      tt <- clin[[eps[[ep]][1]]]; ee <- clin[[eps[[ep]][2]]]
      ok <- is.finite(tt) & tt > 0 & is.finite(ee) & is.finite(scz)
      if (sum(ok) < 50) next
      Sv <- Surv(tt[ok], ee[ok])
      f <- coxph(Sv ~ scz[ok]); s <- summary(f)
      grp <- factor(ifelse(scz[ok] > stats::median(scz[ok]), "high", "low"),
                    levels = c("low", "high"))
      sd_ <- survdiff(Sv ~ grp)
      lr <- stats::pchisq(sd_$chisq, df = length(sd_$n) - 1, lower.tail = FALSE)
      # multivariable: age + stage (+ grade where available)
      dd <- data.frame(time = tt[ok], ev = ee[ok], sc = scz[ok],
                       age = clin$age[ok], stage = clin$stage_group[ok])
      dd <- dd[stats::complete.cases(dd), ]
      mvHR <- mvP <- NA_real_; mvN <- NA_integer_
      if (nrow(dd) > 50 && sum(dd$ev == 1) > 10 && nlevels(droplevels(dd$stage)) > 1) {
        dd$stage <- droplevels(dd$stage)
        f2 <- try(coxph(Surv(time, ev) ~ sc + age + stage, data = dd), silent = TRUE)
        if (!inherits(f2, "try-error")) {
          s2 <- summary(f2); mvHR <- s2$coefficients["sc", "exp(coef)"]
          mvP <- s2$coefficients["sc", "Pr(>|z|)"]; mvN <- nrow(dd)
        }
      }
      i <- i + 1L
      out[[i]] <- data.table(cohort = cohort, score_method = meth, module = mn,
                             endpoint = ep, n = sum(ok), n_event = sum(ee[ok] == 1),
                             uni_HR_per_SD = s$coefficients[1, "exp(coef)"],
                             uni_CI_low = s$conf.int[1, "lower .95"],
                             uni_CI_high = s$conf.int[1, "upper .95"],
                             uni_p = s$coefficients[1, "Pr(>|z|)"],
                             uni_C_index = s$concordance[1], logrank_p = lr,
                             mv_HR_per_SD = mvHR, mv_p = mvP, mv_n = mvN)
    }
  }
  rbindlist(out)
}

mb_mod <- rbind(
  analyse(gs_mb, "GSVA",  "METABRIC", cl, list(OS = c("OS_time", "OS_event"),
                                               RFS = c("RFS_time", "RFS_event"))),
  analyse(mz_mb, "meanZ", "METABRIC", cl, list(OS = c("OS_time", "OS_event"),
                                               RFS = c("RFS_time", "RFS_event"))))

## same protein-coding modules recomputed in TCGA, for a like-for-like comparison
gex <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
pheno <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
tsurv <- as.data.table(readRDS(file.path(BASE, "data/brca_survival_clean.rds")))
tum <- Reduce(intersect, list(colnames(gex),
                              pheno$sample[pheno$sample_type == "Primary Tumor"],
                              tsurv$sample_id))
G <- gex[, tum, drop = FALSE]
vg <- apply(G, 1, function(z) stats::sd(z, na.rm = TRUE)); G <- G[is.finite(vg) & vg > 0, ]
tcl <- tsurv[sample_id %in% tum]; setkey(tcl, sample_id); tcl <- tcl[tum]
stopifnot(identical(tcl$sample_id, colnames(G)))
mods_tc <- lapply(mods_pc, function(s) intersect(s, rownames(G)))
mods_tc <- mods_tc[sapply(mods_tc, length) >= 2]
logf("TCGA protein-coding-restricted scoring on ", ncol(G), " primary tumours")
gs_tc <- gsva(gsvaParam(G, mods_tc, kcdf = "Gaussian", minSize = 2), verbose = FALSE)
Ztc <- t(scale(t(G)))
mz_tc <- t(sapply(mods_tc, function(s) colMeans(Ztc[s, , drop = FALSE], na.rm = TRUE)))
tcl[, OS_time := OS.time]; tcl[, OS_event := OS]
tcl[, RFS_time := PFI.time]; tcl[, RFS_event := PFI]   # PFI is TCGA's RFS analogue
tc_mod <- rbind(
  analyse(gs_tc, "GSVA",  "TCGA_proteincoding_only", tcl,
          list(OS = c("OS_time", "OS_event"), RFS = c("RFS_time", "RFS_event"))),
  analyse(mz_tc, "meanZ", "TCGA_proteincoding_only", tcl,
          list(OS = c("OS_time", "OS_event"), RFS = c("RFS_time", "RFS_event"))))

allmod <- rbind(mb_mod, tc_mod)
fwrite(allmod, file.path(BASE, "results/external_METABRIC_modules.csv"))
logf("WROTE results/external_METABRIC_modules.csv rows=", nrow(allmod))
for (i in seq_len(nrow(allmod))) {
  x <- allmod[i]
  logf(sprintf("%-24s %-5s %-22s %-4s n=%d ev=%d HR=%.3f (%.3f-%.3f) p=%.4g C=%.3f logrank=%.4g | MV HR=%.3f p=%.4g",
               x$cohort, x$score_method, x$module, x$endpoint, x$n, x$n_event,
               x$uni_HR_per_SD, x$uni_CI_low, x$uni_CI_high, x$uni_p, x$uni_C_index,
               x$logrank_p, x$mv_HR_per_SD, x$mv_p))
}

## ---- head-to-head in METABRIC: does higher-order add over 3-node? ----------
h2h <- list(); j <- 0L
for (meth in c("GSVA", "meanZ")) {
  S <- if (meth == "GSVA") gs_mb else mz_mb
  s3  <- as.numeric(scale(as.numeric(S["FFL_3node", ])))
  shi <- as.numeric(scale(as.numeric(S["FFL_higher_order_only", ])))
  for (ep in c("OS", "RFS")) {
    tt <- cl[[paste0(ep, "_time")]]; ee <- cl[[paste0(ep, "_event")]]
    ok <- is.finite(tt) & tt > 0 & is.finite(ee)
    d <- data.frame(time = tt[ok], ev = ee[ok], s3 = s3[ok], shi = shi[ok])
    Sv <- Surv(d$time, d$ev)
    m3 <- coxph(Sv ~ s3, data = d); mh <- coxph(Sv ~ shi, data = d)
    mb2 <- coxph(Sv ~ s3 + shi, data = d)
    a1 <- anova(m3, mb2); a2 <- anova(mh, mb2)
    cc <- concordance(m3, mh)
    dC <- diff(cc$concordance); seD <- sqrt(cc$var[1] + cc$var[2] - 2 * cc$var[3])
    j <- j + 1L
    h2h[[j]] <- data.table(cohort = "METABRIC", score_method = meth, endpoint = ep,
                           n = nrow(d), n_event = sum(d$ev == 1),
                           r_3node_vs_higher = stats::cor(d$s3, d$shi),
                           HR_3node = summary(m3)$coefficients[1, "exp(coef)"],
                           p_3node = summary(m3)$coefficients[1, "Pr(>|z|)"],
                           C_3node = cc$concordance[1],
                           HR_higher = summary(mh)$coefficients[1, "exp(coef)"],
                           p_higher = summary(mh)$coefficients[1, "Pr(>|z|)"],
                           C_higher = cc$concordance[2],
                           deltaC_higher_minus_3node = dC, se_deltaC = seD,
                           p_deltaC = 2 * stats::pnorm(-abs(dC / seD)),
                           LRT_p_higher_adds = a1[2, "Pr(>|Chi|)"],
                           LRT_p_3node_adds  = a2[2, "Pr(>|Chi|)"])
  }
}
h2h <- rbindlist(h2h)
fwrite(h2h, file.path(BASE, "results/external_METABRIC_3node_vs_higher.csv"))
logf("WROTE results/external_METABRIC_3node_vs_higher.csv")
for (i in seq_len(nrow(h2h))) {
  x <- h2h[i]
  logf(sprintf("H2H METABRIC %-5s %-4s | 3node HR=%.3f p=%.3g C=%.4f | higher HR=%.3f p=%.3g C=%.4f | dC=%+.4f p=%.3g | LRT higher-adds p=%.4g, 3node-adds p=%.4g",
               x$score_method, x$endpoint, x$HR_3node, x$p_3node, x$C_3node,
               x$HR_higher, x$p_higher, x$C_higher, x$deltaC_higher_minus_3node,
               x$p_deltaC, x$LRT_p_higher_adds, x$LRT_p_3node_adds))
}
logf("DONE")
