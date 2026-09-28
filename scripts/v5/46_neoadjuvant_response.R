#!/usr/bin/env Rscript
## 46_neoadjuvant_response.R
## Task E(ii): does the Farmer stromal score -- and does the FFL module score -- predict
## pathological complete response to neoadjuvant chemotherapy?  Seven retrievable GEO
## cohorts with pCR outcome; GSE25066 additionally has distant relapse-free survival.
suppressPackageStartupMessages({
  library(data.table); library(GSVA); library(matrixStats)
})
BASE  <- "/path/to/revision"
GEO   <- file.path(BASE, "cache/v5/geo")
LOG   <- file.path(BASE, "logs/v5/46_neoadjuvant_response.log")
logf  <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); set.seed(5)

COHORTS <- list(
  GSE25066 = list(platform="GPL96",  pcr_key="pathologic_response_pcr_rd", pcr_yes="pCR", pcr_no="RD",
                  er_key="er_status_ihc", er_pos="P", er_neg="N",
                  regimen="taxane-anthracycline (T/FAC-type), Hatzis 2011 JAMA", pmid="21558518"),
  GSE20194 = list(platform="GPL96",  pcr_key="pcr_vs_rd", pcr_yes="pCR", pcr_no="RD",
                  er_key="er_status", er_pos="P", er_neg="N",
                  regimen="T/FAC-type, MAQC-II breast set (Shi 2010 Nat Biotechnol)", pmid="20676074"),
  GSE41998 = list(platform="GPL571", pcr_key="pcr", pcr_yes="Yes", pcr_no="No",
                  er_key="er", er_pos="positive", er_neg="negative",
                  regimen="AC then ixabepilone or paclitaxel (Horak 2013 Clin Cancer Res)", pmid="23340299"),
  GSE22093 = list(platform="GPL96",  pcr_key="pcr.v.rd", pcr_yes="pCR", pcr_no="RD",
                  er_key="er positive vs negative by immunohistochemistry", er_pos="ERpos", er_neg="ERneg",
                  regimen="T/FAC or FAC, MDACC+Lyon FNA series (Iwamoto 2011 JNCI)", pmid="21191116"),
  GSE23988 = list(platform="GPL96",  pcr_key="pcr.v.rd", pcr_yes="pCR", pcr_no="RD",
                  er_key="er positive vs negative", er_pos="ERpos", er_neg="ERneg",
                  regimen="T/FAC, US FNA series (Iwamoto 2011 JNCI; Tabchy 2010)", pmid="21191116"),
  GSE32646 = list(platform="GPL570", pcr_key="pathologic response pcr ncr", pcr_yes="pCR", pcr_no="nCR",
                  er_key="er status ihc", er_pos="positive", er_neg="negative",
                  regimen="paclitaxel then FEC (Miyake 2012 Cancer Sci)", pmid="22320227"),
  GSE42822 = list(platform="GPL96",  pcr_key="pcr (1) vs rd (0)", pcr_yes="1", pcr_no="0",
                  er_key=NA, er_pos=NA, er_neg=NA,
                  regimen="FEC/TX +/- trastuzumab, USO series (Shen 2012 BMC Med Genomics)", pmid="23158478")
)

## ------------------------------------------------------------------ series parsing -----
read_series <- function(gse) {
  f <- file.path(GEO, paste0(gse, ".txt.gz"))
  L <- readLines(f, warn = FALSE)
  acc <- strsplit(grep("^!Sample_geo_accession", L, value = TRUE)[1], "\t")[[1]][-1]
  acc <- gsub('"', "", acc)
  ch <- grep("^!Sample_characteristics_ch1", L, value = TRUE)
  meta <- vector("list", length(acc)); names(meta) <- acc
  for (i in seq_along(acc)) meta[[i]] <- list()
  for (line in ch) {
    v <- gsub('"', "", strsplit(line, "\t")[[1]][-1])
    for (i in seq_along(v)) {
      if (!grepl(":", v[i], fixed = TRUE)) next
      kk <- trimws(sub(":.*$", "", v[i])); vv <- trimws(sub("^[^:]*:", "", v[i]))
      meta[[i]][[kk]] <- vv
    }
  }
  st <- grep("^!series_matrix_table_begin", L); en <- grep("^!series_matrix_table_end", L)
  tab <- fread(text = L[(st + 1):(en - 1)], sep = "\t", header = TRUE, data.table = FALSE)
  rn <- gsub('"', "", tab[[1]]); tab <- tab[, -1, drop = FALSE]
  colnames(tab) <- gsub('"', "", colnames(tab))
  X <- as.matrix(tab); rownames(X) <- rn
  X <- X[, acc[acc %in% colnames(X)], drop = FALSE]
  list(X = X, meta = meta[colnames(X)])
}
getf <- function(meta, key) sapply(meta, function(m) if (is.null(m[[key]])) NA_character_ else m[[key]])

## ------------------------------------------------------------------ probe -> symbol ----
annot_map <- function(platform) {
  f <- file.path(GEO, paste0(platform, ".annot.gz"))
  L <- readLines(f, warn = FALSE)
  st <- grep("^!platform_table_begin", L); en <- grep("^!platform_table_end", L)
  a <- fread(text = L[(st + 1):(en - 1)], sep = "\t", header = TRUE, data.table = TRUE,
             quote = "", fill = TRUE)
  a <- a[, .(ID, sym = `Gene symbol`)]
  a <- a[sym != "" & !is.na(sym)]
  a[, sym := sub("///.*$", "", sym)]        # first symbol of multi-mapping probes
  a
}
collapse_to_symbol <- function(X, amap) {
  m <- amap[match(rownames(X), ID)]
  keep <- !is.na(m$sym)
  X <- X[keep, , drop = FALSE]; sym <- m$sym[keep]
  mu <- rowMeans(X, na.rm = TRUE)
  ord <- order(sym, -mu)
  X <- X[ord, , drop = FALSE]; sym <- sym[ord]
  first <- !duplicated(sym)
  Y <- X[first, , drop = FALSE]; rownames(Y) <- sym[first]
  Y
}

sets <- readRDS(file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
KEY  <- c("FARMER_STROMAL","FFL_3node","FFL_higher_order_only","FFL_6node_all","MS_exemplar_module",
          "ESTIMATE_STROMAL_noCOL","CAF_scRNA_50","WEST_DTF_FIBROMATOSIS","WINSLOW_STROMAL_SIG1",
          "FINAK_SDPP","HALLMARK_EMT","HALLMARK_TGF_BETA","CHANG_WOUND_UP_VANTVEER",
          "NABA_CORE_MATRISOME","TRIULZI_ECM","COLLAGEN_PAIR")

auc_fun <- function(score, y) {           # Mann-Whitney AUC
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

store <- list(); rows <- list(); k <- 0L
for (gse in names(COHORTS)) {
  cf <- COHORTS[[gse]]
  logf("\n################ ", gse, "  (", cf$platform, ") — ", cf$regimen)
  s <- read_series(gse)
  X <- s$X
  logf("  probes x samples: ", paste(dim(X), collapse = " x "))
  rng <- range(X, na.rm = TRUE); logf("  value range: ", paste(round(rng, 2), collapse = " .. "))
  if (rng[2] > 60) { X[X < 1] <- 1; X <- log2(X); logf("  -> log2 transformed") }
  X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
  Y <- collapse_to_symbol(X, annot_map(cf$platform))
  Y <- Y[rowSds(Y) > 0, , drop = FALSE]
  logf("  genes after probe collapse: ", nrow(Y),
       " | Farmer genes present: ", length(intersect(sets$FARMER_STROMAL, rownames(Y))), "/50")
  gs <- lapply(sets, function(z) intersect(z, rownames(Y)))
  gs <- gs[sapply(gs, length) >= 2]
  SC <- gsva(ssgseaParam(Y, gs, minSize = 2, normalize = TRUE), verbose = FALSE)
  pcr_raw <- getf(s$meta, cf$pcr_key)
  y <- ifelse(pcr_raw == cf$pcr_yes, 1L, ifelse(pcr_raw == cf$pcr_no, 0L, NA_integer_))
  er <- if (is.na(cf$er_key)) rep(NA_character_, ncol(Y)) else {
    e <- getf(s$meta, cf$er_key)
    ifelse(e == cf$er_pos, "pos", ifelse(e == cf$er_neg, "neg", NA_character_)) }
  logf("  pCR = ", sum(y == 1, na.rm = TRUE), " / evaluable ", sum(!is.na(y)),
       " (", round(100 * mean(y, na.rm = TRUE), 1), "%) | ER+ ", sum(er == "pos", na.rm = TRUE),
       " ER- ", sum(er == "neg", na.rm = TRUE))
  store[[gse]] <- list(SC = SC, y = y, er = er, meta = s$meta, genes = rownames(Y))
  for (nm in intersect(KEY, rownames(SC))) {
    x <- as.numeric(scale(SC[nm, ]))
    ok <- !is.na(y) & is.finite(x)
    if (sum(ok) < 40 || sum(y[ok] == 1) < 5) next
    f <- glm(y[ok] ~ x[ok], family = binomial)
    cs <- summary(f)$coefficients
    ## ER-adjusted
    or_adj <- p_adj <- NA_real_; n_adj <- NA_integer_
    ok2 <- ok & !is.na(er)
    if (sum(ok2) >= 40 && length(unique(er[ok2])) == 2 && sum(y[ok2] == 1) >= 5) {
      f2 <- glm(y[ok2] ~ x[ok2] + factor(er[ok2]), family = binomial)
      c2 <- summary(f2)$coefficients
      or_adj <- exp(c2[2, 1]); p_adj <- c2[2, 4]; n_adj <- sum(ok2)
    }
    k <- k + 1L
    rows[[k]] <- data.table(cohort = gse, platform = cf$platform, regimen = cf$regimen,
      pmid = cf$pmid, score = nm, n = sum(ok), n_pCR = sum(y[ok] == 1),
      logOR_per_SD = cs[2, 1], se = cs[2, 2], OR_per_SD = exp(cs[2, 1]),
      OR_low = exp(cs[2, 1] - 1.96 * cs[2, 2]), OR_high = exp(cs[2, 1] + 1.96 * cs[2, 2]),
      p = cs[2, 4], AUC = auc_fun(x[ok], y[ok]),
      OR_per_SD_ER_adjusted = or_adj, p_ER_adjusted = p_adj, n_ER_adjusted = n_adj)
  }
}
res <- rbindlist(rows)
res[, q_BH := p.adjust(p, "BH"), by = cohort]
fwrite(res, file.path(BASE, "results/v5/farmer_neoadjuvant_pcr.csv"))
logf("\nWROTE results/v5/farmer_neoadjuvant_pcr.csv rows=", nrow(res))
saveRDS(store, file.path(BASE, "cache/v5/farmer_neoadjuvant_store.rds"))

## -------------------------------------------------------------------- meta-analysis ----
meta_one <- function(d) {
  w <- 1 / d$se^2; mu <- sum(w * d$logOR_per_SD) / sum(w); se <- sqrt(1 / sum(w))
  Q <- sum(w * (d$logOR_per_SD - mu)^2); df <- nrow(d) - 1
  tau2 <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (d$se^2 + tau2); mur <- sum(wr * d$logOR_per_SD) / sum(wr); ser <- sqrt(1 / sum(wr))
  data.table(k = nrow(d), n_total = sum(d$n), n_pCR_total = sum(d$n_pCR),
    OR_fixed = exp(mu), fixed_lo = exp(mu - 1.96 * se), fixed_hi = exp(mu + 1.96 * se),
    p_fixed = 2 * pnorm(-abs(mu / se)),
    OR_random = exp(mur), random_lo = exp(mur - 1.96 * ser), random_hi = exp(mur + 1.96 * ser),
    p_random = 2 * pnorm(-abs(mur / ser)),
    Q = Q, df = df, p_hetero = pchisq(Q, df, lower.tail = FALSE),
    I2 = max(0, 100 * (Q - df) / Q), tau2 = tau2)
}
mt <- res[, meta_one(.SD), by = score]
setorder(mt, OR_fixed)
fwrite(mt, file.path(BASE, "results/v5/farmer_neoadjuvant_meta.csv"))
logf("\n=== meta-analysis across the ", length(unique(res$cohort)), " neoadjuvant cohorts ",
     "(OR of pCR per +1 SD of score; OR<1 = resistance) ===")
for (i in seq_len(nrow(mt))) { r <- mt[i]
  logf(sprintf("  %-26s k=%d n=%4d pCR=%3d | fixed OR=%.3f (%.3f-%.3f) p=%.3g | random OR=%.3f (%.3f-%.3f) p=%.3g | I2=%.0f%%",
    r$score, r$k, r$n_total, r$n_pCR_total, r$OR_fixed, r$fixed_lo, r$fixed_hi, r$p_fixed,
    r$OR_random, r$random_lo, r$random_hi, r$p_random, r$I2)) }
logf("\n=== per-cohort detail for the Farmer stromal score and the FFL modules ===")
for (nm in c("FARMER_STROMAL","FFL_3node","FFL_higher_order_only","FFL_6node_all","COLLAGEN_PAIR")) {
  for (i in which(res$score == nm)) { r <- res[i]
    logf(sprintf("  %-26s %-9s n=%3d pCR=%2d OR=%.3f (%.3f-%.3f) p=%.3g AUC=%.3f | ER-adj OR=%.3f p=%.3g",
      nm, r$cohort, r$n, r$n_pCR, r$OR_per_SD, r$OR_low, r$OR_high, r$p, r$AUC,
      r$OR_per_SD_ER_adjusted, r$p_ER_adjusted)) }
}
logf("DONE 46 (part 1: pCR)")
