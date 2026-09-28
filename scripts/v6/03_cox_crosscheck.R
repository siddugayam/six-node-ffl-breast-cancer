# Cross-check HPA pathology-atlas prognostic calls against our own TCGA-BRCA and
# METABRIC Cox / best-cutoff results for the same genes.
suppressPackageStartupMessages({library(survival)})
setwd("/path/to/revision")
gs <- c("E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3",
        "CCND2","COL1A1","STAT5A","MYBL2","FN1","PDGFRB","MET","CXCL12","MMP14","PLAU",
        "COL3A1","POSTN","NFKB1","RELA","SP1","ETS1")

fit_one <- function(x, time, event){
  ok <- is.finite(x) & is.finite(time) & is.finite(event) & time > 0
  x <- x[ok]; time <- time[ok]; event <- event[ok]
  if (length(unique(x)) < 5 || sum(event) < 10) return(NULL)
  z <- as.numeric(scale(x))
  cx <- coxph(Surv(time, event) ~ z)
  s  <- summary(cx)
  # HPA-style best cut-off scan: percentile 20..80, min log-rank p
  qs <- quantile(x, probs = seq(0.20, 0.80, by = 0.01), na.rm = TRUE)
  best_p <- NA_real_; best_cut <- NA_real_; best_dir <- NA_character_
  for (cut in unique(qs)) {
    grp <- x > cut
    if (min(table(grp)) < 10) next
    sd_ <- survdiff(Surv(time, event) ~ grp)
    p <- pchisq(sd_$chisq, df = 1, lower.tail = FALSE)
    if (is.na(best_p) || p < best_p) {
      best_p <- p; best_cut <- cut
      # observed vs expected: dir = "unfavourable" if high group has more deaths than expected
      oe_high <- sd_$obs[2] / sd_$exp[2]
      best_dir <- if (oe_high > 1) "unfavourable" else "favourable"
    }
  }
  data.frame(n = length(x), n_event = sum(event),
             HR_per_SD = unname(s$coefficients[1, "exp(coef)"]),
             CI_low  = unname(s$conf.int[1, "lower .95"]),
             CI_high = unname(s$conf.int[1, "upper .95"]),
             p_cox = unname(s$coefficients[1, "Pr(>|z|)"]),
             C_index = unname(s$concordance[1]),
             bestcut_p = best_p, bestcut_quantile = mean(x <= best_cut),
             bestcut_direction = best_dir,
             cox_direction = ifelse(unname(s$coefficients[1,"coef"]) > 0, "unfavourable", "favourable"),
             stringsAsFactors = FALSE)
}

res <- list()

## ---- TCGA-BRCA ----
ge <- readRDS("data/brca_gene_expr_symbol.rds")
sv <- readRDS("data/brca_survival.rds")
rownames(sv) <- sv$sample
common <- intersect(colnames(ge), rownames(sv))
cat("TCGA overlap samples:", length(common), "\n")
sv2 <- sv[common, ]
for (g in gs) {
  if (!g %in% rownames(ge)) next
  x <- as.numeric(ge[g, common])
  for (ep in c("OS","DSS")) {
    tt <- as.numeric(sv2[[paste0(ep, ".time")]]); ee <- as.numeric(sv2[[ep]])
    r <- fit_one(x, tt, ee)
    if (!is.null(r)) res[[length(res)+1]] <- cbind(cohort = "TCGA-BRCA", gene = g, endpoint = ep, r)
  }
}

## ---- METABRIC ----
mb <- readRDS("data/metabric.rds")
M <- mb$M; cl <- as.data.frame(mb$cl)
sid <- if ("SAMPLE_ID" %in% names(cl)) cl$SAMPLE_ID else cl$PATIENT_ID
rownames(cl) <- sid
common2 <- intersect(colnames(M), rownames(cl))
cat("METABRIC overlap samples:", length(common2), "\n")
cl2 <- cl[common2, ]
for (g in gs) {
  if (!g %in% rownames(M)) next
  x <- as.numeric(M[g, common2])
  for (ep in c("OS","RFS")) {
    tt <- as.numeric(cl2[[paste0(ep, "_time")]]); ee <- as.numeric(cl2[[paste0(ep, "_event")]])
    r <- fit_one(x, tt, ee)
    if (!is.null(r)) res[[length(res)+1]] <- cbind(cohort = "METABRIC", gene = g, endpoint = ep, r)
  }
}

out <- do.call(rbind, res)
out$q_BH <- NA_real_
for (k in unique(paste(out$cohort, out$endpoint))) {
  i <- paste(out$cohort, out$endpoint) == k
  out$q_BH[i] <- p.adjust(out$p_cox[i], "BH")
}
dir.create("results/v6", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "results/v6/hpa_crosscheck_our_cox.csv", row.names = FALSE)
cat("rows:", nrow(out), "\n")
print(head(out, 8))
