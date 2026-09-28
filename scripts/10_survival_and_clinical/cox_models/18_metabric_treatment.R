#!/usr/bin/env Rscript
# 18_metabric_treatment.R
# METABRIC records which patients actually received chemotherapy, hormone therapy and
# radiotherapy, so the closest legitimate computational test is:
#   (i)  is the hub / module score prognostic WITHIN the treated subgroup, and
#   (ii) is there a score x treatment INTERACTION on overall survival
#        (i.e. does the score identify patients who do worse despite treatment)?
# A significant interaction is the minimum evidence needed before the word
# "resistance" can be used. Output: results/metabric_treatment_interaction.csv

suppressPackageStartupMessages({ library(data.table); library(survival) })
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs", "18_metabric_treatment.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")

d <- readRDS(file.path(BASE, "data/metabric.rds"))
M <- d$M; cl <- copy(d$cl)
stopifnot(identical(colnames(M), cl$SAMPLE_ID))

# module scores: mean z-score over the protein-coding members (METABRIC has no miRNA
# assay). mean-z is used here rather than GSVA because script 12/15 showed the two
# give the same direction, and mean-z is cheap and fully deterministic.
sets  <- readRDS(file.path(BASE, "data/ffl_module_sets.rds"))
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
prot  <- nodes[type != "miRNA", name]
mods <- lapply(list(MS_exemplar_module    = sets$MS_MODULE,
                    FFL_3node             = sets$N3,
                    FFL_higher_order_only = sets$HIGHER_ONLY,
                    FFL_6node_all         = sets$N6),
               function(s) intersect(intersect(s, prot), rownames(M)))
mods <- mods[sapply(mods, length) >= 2]
for (nm in names(mods)) logf("module ", nm, " members in METABRIC = ", length(mods[[nm]]))
Zmb <- t(scale(t(M)))
gs <- t(sapply(mods, function(s) colMeans(Zmb[s, , drop = FALSE], na.rm = TRUE)))
logf("module score matrix (mean-z): ", paste(dim(gs), collapse = " x "))

for (v in c("CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY")) {
  cl[[v]] <- as.character(cl[[v]])
  logf(v, ": ", paste(sprintf("%s=%d", names(table(cl[[v]], useNA = "ifany")),
                              table(cl[[v]], useNA = "ifany")), collapse = " "))
}

named_hubs <- c("NFKB1", "RELA", "SP1", "ETS1", "COL1A1", "COL3A1", "VEGFA",
                "CCND2", "MYC", "E2F1", "TP53")
feats <- list()
for (h in intersect(named_hubs, rownames(M))) feats[[h]] <- as.numeric(M[h, ])
for (m in rownames(gs)) feats[[paste0("MODULE_meanZ:", m)]] <- as.numeric(gs[m, ])
logf("features tested: ", length(feats))

rows <- list(); k <- 0L
for (tx in c("CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY")) {
  trt <- factor(ifelse(cl[[tx]] == "YES", "treated", ifelse(cl[[tx]] == "NO", "untreated", NA)),
                levels = c("untreated", "treated"))
  for (ep in c("OS", "RFS")) {
    tt <- cl[[paste0(ep, "_time")]]; ee <- cl[[paste0(ep, "_event")]]
    for (fn in names(feats)) {
      x <- feats[[fn]]
      ok <- is.finite(tt) & tt > 0 & is.finite(ee) & is.finite(x) & !is.na(trt)
      if (sum(ok) < 100) next
      xs <- as.numeric(scale(x[ok]))
      dd <- data.frame(time = tt[ok], ev = ee[ok], x = xs, trt = droplevels(trt[ok]),
                       age = cl$age[ok])
      if (nlevels(dd$trt) < 2) next
      fi <- try(coxph(Surv(time, ev) ~ x * trt + age, data = dd), silent = TRUE)
      if (inherits(fi, "try-error")) next
      si <- summary(fi)
      irow <- grep(":", rownames(si$coefficients), value = TRUE)[1]
      # stratum-specific effects
      eff <- lapply(levels(dd$trt), function(lv) {
        s <- dd[dd$trt == lv, ]
        if (nrow(s) < 50 || sum(s$ev == 1) < 10) return(c(NA, NA, NA))
        f <- try(coxph(Surv(time, ev) ~ x, data = s), silent = TRUE)
        if (inherits(f, "try-error")) return(c(NA, NA, NA))
        ss <- summary(f)
        c(ss$coefficients[1, "exp(coef)"], ss$coefficients[1, "Pr(>|z|)"], nrow(s))
      })
      names(eff) <- levels(dd$trt)
      k <- k + 1L
      rows[[k]] <- data.table(
        treatment = tx, endpoint = ep, feature = fn, n = nrow(dd),
        n_event = sum(dd$ev == 1),
        HR_untreated = eff[["untreated"]][1], p_untreated = eff[["untreated"]][2],
        n_untreated = eff[["untreated"]][3],
        HR_treated = eff[["treated"]][1], p_treated = eff[["treated"]][2],
        n_treated = eff[["treated"]][3],
        interaction_HR = si$coefficients[irow, "exp(coef)"],
        interaction_p  = si$coefficients[irow, "Pr(>|z|)"])
    }
  }
}
res <- rbindlist(rows)
res[, interaction_q_BH := p.adjust(interaction_p, method = "BH")]
setorder(res, interaction_p)
fwrite(res, file.path(BASE, "results/metabric_treatment_interaction.csv"))
logf("WROTE results/metabric_treatment_interaction.csv rows=", nrow(res))
logf("interaction tests: ", nrow(res), " ; nominal p<0.05 = ", sum(res$interaction_p < 0.05),
     " ; BH q<0.05 = ", sum(res$interaction_q_BH < 0.05))
for (i in seq_len(min(20, nrow(res)))) {
  x <- res[i]
  logf(sprintf("%-16s %-4s %-28s untreated HR=%.3f (p=%.3g, n=%d) | treated HR=%.3f (p=%.3g, n=%d) | interaction HR=%.3f p=%.4g q=%.3g",
               x$treatment, x$endpoint, x$feature, x$HR_untreated, x$p_untreated,
               as.integer(x$n_untreated), x$HR_treated, x$p_treated,
               as.integer(x$n_treated), x$interaction_HR, x$interaction_p,
               x$interaction_q_BH))
}
logf("DONE")
