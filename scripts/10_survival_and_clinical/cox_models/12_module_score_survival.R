#!/usr/bin/env Rscript
# 12_module_score_survival.R
# PART 1B  module-level scores (GSVA + mean-z) -> Cox, KM median split,
#          multivariable Cox adjusted for age and AJCC stage.
# PART 1C  head-to-head: 3-node-FFL score vs higher-order-only (>=4-node) score.
#          C-index, nested-model likelihood-ratio test, added information.
# Output: results/module_score_survival.csv (+ KM figures)

suppressPackageStartupMessages({
  library(data.table); library(survival); library(GSVA); library(Biobase)
})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs", "12_module_score_survival.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")
set.seed(1234)

sets  <- readRDS(file.path(BASE, "data/ffl_module_sets.rds"))
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
ntype <- setNames(nodes$type, nodes$name)

gex <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mex <- readRDS(file.path(BASE, "data/brca_mirna_expr_canonical.rds"))
pheno <- readRDS(file.path(BASE, "data/brca_pheno.rds"))
surv  <- as.data.table(readRDS(file.path(BASE, "data/brca_survival_clean.rds")))

tumour <- pheno$sample[pheno$sample_type == "Primary Tumor"]
common <- Reduce(intersect, list(colnames(gex), colnames(mex), tumour, surv$sample_id))
logf("samples with gene + miRNA + survival, primary tumour: ", length(common))

X <- rbind(gex[, common, drop = FALSE], mex[, common, drop = FALSE])
# drop zero-variance rows (GSVA cannot rank them)
v <- apply(X, 1, function(z) stats::sd(z, na.rm = TRUE))
X <- X[is.finite(v) & v > 0, , drop = FALSE]
logf("combined expression matrix for scoring: ", paste(dim(X), collapse = " x "),
     " (", nrow(gex), " gene rows + ", nrow(mex), " miRNA rows before variance filter)")

## ------------------------------------------------------------- module defs --
mods <- list(
  MS_exemplar_module      = sets$MS_MODULE,
  MS_exemplar_TFs         = c("NFKB1", "RELA", "SP1"),
  MS_exemplar_collagens   = c("COL1A1", "COL3A1"),
  MS_exemplar_miRNAs      = c("hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c",
                              "hsa-let-7b", "hsa-let-7e"),
  FFL_3node               = sets$N3,
  FFL_higher_order_only   = sets$HIGHER_ONLY,
  FFL_6node_all           = sets$N6,
  FFL_4node_added         = setdiff(sets$N4, sets$N3),
  FFL_5node_added         = setdiff(sets$N5, sets$N4),
  FFL_6node_added         = setdiff(sets$N6, sets$N5)
)
mods <- lapply(mods, function(s) intersect(s, rownames(X)))
for (nmm in names(mods)) logf(sprintf("module %-24s n_members_with_expression=%d",
                                      nmm, length(mods[[nmm]])))
mods <- mods[sapply(mods, length) >= 2]

## ------------------------------------------------------------------- GSVA ---
logf("running GSVA (kcdf=Gaussian) ...")
par <- gsvaParam(X, mods, kcdf = "Gaussian", minSize = 2, maxSize = Inf)
gs  <- gsva(par, verbose = FALSE)
logf("GSVA scores: ", paste(dim(gs), collapse = " x "))

## ------------------------------------------------ mean-z sensitivity score --
Z <- t(scale(t(X)))                      # per-feature z across samples
meanz <- t(sapply(mods, function(s) colMeans(Z[s, , drop = FALSE], na.rm = TRUE)))
logf("mean-z scores: ", paste(dim(meanz), collapse = " x "))

score_tabs <- list(GSVA = gs, meanZ = meanz)

## ------------------------------------------------------------- clinical -----
cl <- surv[sample_id %in% common]
setkey(cl, sample_id)
cl <- cl[common]
stopifnot(identical(cl$sample_id, common))
logf("clinical rows aligned: ", nrow(cl),
     " | age non-missing=", sum(!is.na(cl$age)),
     " | stage non-missing=", sum(!is.na(cl$stage_group)))
logf("events among these samples: OS=", sum(cl$OS == 1, na.rm = TRUE),
     " PFI=", sum(cl$PFI == 1, na.rm = TRUE), " DSS=", sum(cl$DSS == 1, na.rm = TRUE))

## --------------------------------------------------------------- analysis ---
rows <- list(); r <- 0L
km_store <- list()

for (meth in names(score_tabs)) {
  S <- score_tabs[[meth]]
  for (mname in rownames(S)) {
    sc <- as.numeric(S[mname, ])
    scz <- as.numeric(scale(sc))
    for (ep in c("OS", "PFI", "DSS")) {
      tt <- cl[[paste0(ep, ".time")]]; ee <- cl[[ep]]
      ok <- is.finite(tt) & tt > 0 & is.finite(ee) & is.finite(scz)
      if (sum(ok) < 50) next
      Sv <- Surv(tt[ok], ee[ok])

      ## univariable
      f1 <- coxph(Sv ~ scz[ok]); s1 <- summary(f1)

      ## median split + log-rank
      grp <- factor(ifelse(scz[ok] > stats::median(scz[ok]), "high", "low"),
                    levels = c("low", "high"))
      sd_ <- survdiff(Sv ~ grp)
      lr_p <- stats::pchisq(sd_$chisq, df = length(sd_$n) - 1, lower.tail = FALSE)
      fkm <- coxph(Sv ~ grp); skm <- summary(fkm)

      ## multivariable: + age + stage
      dfm <- data.frame(time = tt[ok], ev = ee[ok], sc = scz[ok],
                        age = cl$age[ok], stage = cl$stage_group[ok])
      dfm2 <- dfm[stats::complete.cases(dfm), ]
      mv_hr <- mv_lo <- mv_hi <- mv_p <- NA_real_; mv_n <- NA_integer_; mv_ev <- NA_integer_
      mv_c <- NA_real_
      if (nrow(dfm2) > 50 && sum(dfm2$ev == 1) > 10 &&
          nlevels(droplevels(dfm2$stage)) > 1) {
        f2 <- try(coxph(Surv(time, ev) ~ sc + age + stage, data = dfm2), silent = TRUE)
        if (!inherits(f2, "try-error")) {
          s2 <- summary(f2)
          mv_hr <- s2$coefficients["sc", "exp(coef)"]
          mv_lo <- s2$conf.int["sc", "lower .95"]
          mv_hi <- s2$conf.int["sc", "upper .95"]
          mv_p  <- s2$coefficients["sc", "Pr(>|z|)"]
          mv_n  <- nrow(dfm2); mv_ev <- sum(dfm2$ev == 1)
          mv_c  <- s2$concordance[1]
        }
      }

      r <- r + 1L
      rows[[r]] <- data.table(
        score_method = meth, module = mname, n_members = length(mods[[mname]]),
        endpoint = ep, n = sum(ok), n_event = sum(ee[ok] == 1),
        uni_HR_per_SD = s1$coefficients[1, "exp(coef)"],
        uni_CI_low  = s1$conf.int[1, "lower .95"],
        uni_CI_high = s1$conf.int[1, "upper .95"],
        uni_p = s1$coefficients[1, "Pr(>|z|)"],
        uni_C_index = s1$concordance[1],
        median_split_HR = skm$coefficients[1, "exp(coef)"],
        median_split_CI_low  = skm$conf.int[1, "lower .95"],
        median_split_CI_high = skm$conf.int[1, "upper .95"],
        logrank_p = lr_p,
        mv_HR_per_SD = mv_hr, mv_CI_low = mv_lo, mv_CI_high = mv_hi,
        mv_p = mv_p, mv_n = mv_n, mv_n_event = mv_ev, mv_C_index = mv_c)

      if (meth == "GSVA" && ep %in% c("OS", "PFI"))
        km_store[[paste(mname, ep, sep = "|")]] <- list(Sv = Sv, grp = grp, p = lr_p)
    }
  }
}
res <- rbindlist(rows)
res[, uni_q_BH := p.adjust(uni_p, method = "BH"), by = .(score_method, endpoint)]
setorder(res, score_method, endpoint, uni_p)
outf <- file.path(BASE, "results/module_score_survival.csv")
fwrite(res, outf)
logf("WROTE ", outf, " rows=", nrow(res))

logf("======== PART 1B: module Cox (GSVA) ========")
for (ep in c("OS", "PFI", "DSS")) {
  sub <- res[score_method == "GSVA" & endpoint == ep][order(uni_p)]
  for (i in seq_len(nrow(sub)))
    logf(sprintf("%-4s %-24s n=%d ev=%d  uniHR=%.3f (%.3f-%.3f) p=%.4g q=%.3g C=%.3f | logrank p=%.4g | MV(age+stage) HR=%.3f p=%.4g",
                 ep, sub$module[i], sub$n[i], sub$n_event[i], sub$uni_HR_per_SD[i],
                 sub$uni_CI_low[i], sub$uni_CI_high[i], sub$uni_p[i], sub$uni_q_BH[i],
                 sub$uni_C_index[i], sub$logrank_p[i], sub$mv_HR_per_SD[i], sub$mv_p[i]))
}
logf("======== sanity: mean-z gives the same direction? ========")
chk <- merge(res[score_method == "GSVA",  .(module, endpoint, HR_gsva = uni_HR_per_SD, p_gsva = uni_p)],
             res[score_method == "meanZ", .(module, endpoint, HR_meanz = uni_HR_per_SD, p_meanz = uni_p)],
             by = c("module", "endpoint"))
chk[, same_direction := sign(log(HR_gsva)) == sign(log(HR_meanz))]
logf("same direction in ", sum(chk$same_direction), "/", nrow(chk), " module-endpoint pairs")
for (i in seq_len(nrow(chk)))
  logf(sprintf("   %-24s %-4s GSVA HR=%.3f p=%.3g | meanZ HR=%.3f p=%.3g | same_dir=%s",
               chk$module[i], chk$endpoint[i], chk$HR_gsva[i], chk$p_gsva[i],
               chk$HR_meanz[i], chk$p_meanz[i], chk$same_direction[i]))

## =================== PART 1C: 3-node vs higher-order head-to-head ===========
logf("======== PART 1C: 3-node vs higher-order (>=4-node) head-to-head ========")
h2h <- list(); h <- 0L
for (meth in names(score_tabs)) {
  S <- score_tabs[[meth]]
  s3  <- as.numeric(scale(as.numeric(S["FFL_3node", ])))
  shi <- as.numeric(scale(as.numeric(S["FFL_higher_order_only", ])))
  cor_s <- stats::cor(s3, shi, method = "pearson")
  logf(sprintf("[%s] Pearson r(3-node score, higher-order-only score) = %.3f", meth, cor_s))
  for (ep in c("OS", "PFI", "DSS")) {
    tt <- cl[[paste0(ep, ".time")]]; ee <- cl[[ep]]
    ok <- is.finite(tt) & tt > 0 & is.finite(ee)
    d <- data.frame(time = tt[ok], ev = ee[ok], s3 = s3[ok], shi = shi[ok],
                    age = cl$age[ok], stage = cl$stage_group[ok])
    d <- d[stats::complete.cases(d[, c("time", "ev", "s3", "shi")]), ]
    Sv <- Surv(d$time, d$ev)
    m3   <- coxph(Sv ~ s3,        data = d)
    mhi  <- coxph(Sv ~ shi,       data = d)
    mboth<- coxph(Sv ~ s3 + shi,  data = d)
    # nested LRTs
    lrt_add_hi <- anova(m3,  mboth)      # does higher-order add over 3-node?
    lrt_add_3  <- anova(mhi, mboth)      # does 3-node add over higher-order?
    # C-index with SE, and a paired comparison of the two univariable models
    c3  <- concordance(m3);  chi <- concordance(mhi); cboth <- concordance(mboth)
    cc  <- concordance(m3, mhi)          # joint -> covariance for a paired test
    dC  <- diff(cc$concordance)
    seD <- sqrt(cc$var[1] + cc$var[2] - 2 * cc$var[3])
    z_c <- dC / seD
    p_c <- 2 * stats::pnorm(-abs(z_c))
    # clinical baseline and incremental value over clinical
    dcl <- d[stats::complete.cases(d), ]
    p_cl3 <- p_clhi <- NA_real_; c_cl <- c_cl3 <- c_clhi <- NA_real_
    if (nrow(dcl) > 50 && nlevels(droplevels(dcl$stage)) > 1) {
      Sc <- Surv(dcl$time, dcl$ev)
      mcl   <- coxph(Sc ~ age + stage,       data = dcl)
      mcl3  <- coxph(Sc ~ age + stage + s3,  data = dcl)
      mclhi <- coxph(Sc ~ age + stage + shi, data = dcl)
      p_cl3  <- anova(mcl, mcl3)[2, "Pr(>|Chi|)"]
      p_clhi <- anova(mcl, mclhi)[2, "Pr(>|Chi|)"]
      c_cl   <- concordance(mcl)$concordance
      c_cl3  <- concordance(mcl3)$concordance
      c_clhi <- concordance(mclhi)$concordance
    }
    h <- h + 1L
    h2h[[h]] <- data.table(
      score_method = meth, endpoint = ep, n = nrow(d), n_event = sum(d$ev == 1),
      r_3node_vs_higher = cor_s,
      HR_3node = summary(m3)$coefficients[1, "exp(coef)"],
      p_3node  = summary(m3)$coefficients[1, "Pr(>|z|)"],
      C_3node  = c3$concordance,
      HR_higher = summary(mhi)$coefficients[1, "exp(coef)"],
      p_higher  = summary(mhi)$coefficients[1, "Pr(>|z|)"],
      C_higher  = chi$concordance,
      C_both    = cboth$concordance,
      deltaC_higher_minus_3node = dC, se_deltaC = seD, p_deltaC = p_c,
      LRT_chisq_higher_adds = lrt_add_hi[2, "Chisq"],
      LRT_p_higher_adds     = lrt_add_hi[2, "Pr(>|Chi|)"],
      LRT_chisq_3node_adds  = lrt_add_3[2, "Chisq"],
      LRT_p_3node_adds      = lrt_add_3[2, "Pr(>|Chi|)"],
      C_clinical = c_cl, C_clinical_plus_3node = c_cl3, C_clinical_plus_higher = c_clhi,
      LRT_p_3node_over_clinical = p_cl3, LRT_p_higher_over_clinical = p_clhi)
  }
}
h2h <- rbindlist(h2h)
outf2 <- file.path(BASE, "results/module_3node_vs_higher_order.csv")
fwrite(h2h, outf2)
logf("WROTE ", outf2, " rows=", nrow(h2h))
for (i in seq_len(nrow(h2h))) {
  x <- h2h[i]
  logf(sprintf("[%s] %-4s n=%d ev=%d | 3node HR=%.3f p=%.3g C=%.4f | higher HR=%.3f p=%.3g C=%.4f | dC=%+.4f (SE %.4f) p=%.3g",
               x$score_method, x$endpoint, x$n, x$n_event, x$HR_3node, x$p_3node, x$C_3node,
               x$HR_higher, x$p_higher, x$C_higher, x$deltaC_higher_minus_3node,
               x$se_deltaC, x$p_deltaC))
  logf(sprintf("        LRT higher-adds-over-3node chisq=%.3f p=%.4g | 3node-adds-over-higher chisq=%.3f p=%.4g",
               x$LRT_chisq_higher_adds, x$LRT_p_higher_adds,
               x$LRT_chisq_3node_adds, x$LRT_p_3node_adds))
  logf(sprintf("        clinical C=%.4f -> +3node C=%.4f (LRT p=%.4g) | +higher C=%.4f (LRT p=%.4g)",
               x$C_clinical, x$C_clinical_plus_3node, x$LRT_p_3node_over_clinical,
               x$C_clinical_plus_higher, x$LRT_p_higher_over_clinical))
}

## ------------------------------------------------------------ KM figures ----
dir.create(file.path(BASE, "figures"), showWarnings = FALSE)
pdf(file.path(BASE, "figures/KM_module_scores.pdf"), width = 11, height = 8)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (k in names(km_store)) {
  z <- km_store[[k]]
  fit <- survfit(z$Sv ~ z$grp)
  plot(fit, col = c("#2c7fb8", "#d95f02"), lwd = 2, xlab = "days", ylab = "survival prob",
       main = sprintf("%s\nlog-rank p = %.3g", k, z$p))
  legend("bottomleft", legend = c("low", "high"), col = c("#2c7fb8", "#d95f02"),
         lwd = 2, bty = "n", cex = 0.8)
}
par(op); dev.off()
logf("WROTE figures/KM_module_scores.pdf")
saveRDS(list(gsva = gs, meanz = meanz, samples = common, mods = mods),
        file.path(BASE, "data/tcga_module_scores.rds"))
logf("DONE")
