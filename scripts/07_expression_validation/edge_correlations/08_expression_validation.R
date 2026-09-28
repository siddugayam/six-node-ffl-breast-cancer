#!/usr/bin/env Rscript
# 08_expression_validation.R
# Expression-based validation of PREDICTED REGULATORY RELATIONSHIPS in TCGA-BRCA.
# ("differential expression cannot demonstrate that the proposed FFLs regulate anything").
#
# A) edge-level Spearman correlation of every predicted edge -> results/edge_correlation.csv
# B) degree- + expression-decile-matched random-pair null -> results/null_comparison_summary.csv
# C) sign-concordance test (replaces the inappropriate ANOVA) -> results/sign_concordance_summary.csv

REV <- "/path/to/revision"
LOG <- file.path(REV, "logs", "expr_validation.log")
dir.create(dirname(LOG), showWarnings = FALSE, recursive = TRUE)
con <- file(LOG, open = "wt")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n"); cat(msg, "\n", file = con); flush(con)
}
set.seed(1234)
say("=== 08_expression_validation.R START ===")

## ---------------------------------------------------------------- 1. LOAD
G0 <- readRDS(file.path(REV, "data", "brca_gene_expr.rds"))
M0 <- readRDS(file.path(REV, "data", "brca_mirna_expr_canonical.rds"))
ph <- readRDS(file.path(REV, "data", "brca_pheno.rds"))
edges <- read.delim(file.path(REV, "data", "canonical_edges.tsv"), stringsAsFactors = FALSE)
nodes <- read.delim(file.path(REV, "data", "canonical_nodes.tsv"), stringsAsFactors = FALSE)
tier  <- read.delim(file.path(REV, "data", "edge_evidence_tier.tsv"), stringsAsFactors = FALSE)
cov_tf  <- read.delim(file.path(REV, "results", "TF_target_TRRUST_coverage.tsv"), stringsAsFactors = FALSE)
cov_tfm <- read.delim(file.path(REV, "results", "TF_miRNA_TransmiR_coverage.tsv"), stringsAsFactors = FALSE)
say("gene matrix ", paste(dim(G0), collapse = "x"), "  miRNA matrix ", paste(dim(M0), collapse = "x"))
say("canonical edges ", nrow(edges), "  nodes ", nrow(nodes), "  tier rows ", nrow(tier))

tumour <- ph$sample[ph$sample_type == "Primary Tumor"]
sm <- sort(intersect(intersect(colnames(G0), tumour), intersect(colnames(M0), tumour)))
N_PAIRED <- length(sm)
say("Primary Tumor: gene assay ", sum(colnames(G0) %in% tumour), ", miRNA assay ",
    sum(colnames(M0) %in% tumour), ", PAIRED (both assays) n = ", N_PAIRED)
G <- G0[, sm, drop = FALSE]; M <- M0[, sm, drop = FALSE]
rm(G0, M0); invisible(gc())

sdG <- apply(G, 1, sd); sdM <- apply(M, 1, sd)
node_type   <- setNames(nodes$type, nodes$name)
gene_nodes  <- nodes$name[nodes$type %in% c("TF", "Gene")]
mirna_nodes <- nodes$name[nodes$type == "miRNA"]
gene_ok  <- intersect(gene_nodes,  rownames(G)[sdG > 0])
mirna_ok <- intersect(mirna_nodes, rownames(M)[sdM > 0])
say("usable network nodes: gene/TF ", length(gene_ok), "/", length(gene_nodes),
    " ; miRNA ", length(mirna_ok), "/", length(mirna_nodes))
say("unusable miRNA nodes: ", paste(setdiff(mirna_nodes, mirna_ok), collapse = ", "))
say("unusable gene nodes : ", paste(setdiff(gene_nodes, gene_ok), collapse = ", "))

## ------------------------------------------------- 2. RANKS -> SPEARMAN
zrow <- function(X) { mu <- rowMeans(X); s <- sqrt(rowSums((X - mu)^2)); (X - mu) / s }
Z <- rbind(zrow(t(apply(G[gene_ok, , drop = FALSE], 1, rank))),
           zrow(t(apply(M[mirna_ok, , drop = FALSE], 1, rank))))
say("z-scored rank matrix: ", nrow(Z), " features x ", ncol(Z), " samples")
rho_of <- function(a, b) {
  ok <- a %in% rownames(Z) & b %in% rownames(Z)
  out <- rep(NA_real_, length(a))
  if (any(ok)) out[ok] <- rowSums(Z[a[ok], , drop = FALSE] * Z[b[ok], , drop = FALSE])
  pmax(pmin(out, 1), -1)
}
p_of <- function(rho, n) {
  tt <- rho * sqrt((n - 2) / pmax(1 - rho^2, .Machine$double.xmin))
  2 * pt(-abs(tt), df = n - 2)
}
ct <- suppressWarnings(cor.test(G["COL1A1", ], G["COL3A1", ], method = "spearman", exact = FALSE))
say("SANITY cor.test rho(COL1A1,COL3A1)=", signif(unname(ct$estimate), 8),
    " vectorised=", signif(rho_of("COL1A1", "COL3A1"), 8),
    " ; p cor.test=", signif(ct$p.value, 4),
    " p vectorised=", signif(p_of(rho_of("COL1A1", "COL3A1"), N_PAIRED), 4))

## ------------------------------------------------- 3. PREDICTED SIGN PER EDGE
key <- function(s, t) paste(s, t, sep = "\r")
E <- edges
E$predicted_sign <- NA_integer_; E$sign_source <- NA_character_; E$annotation <- NA_character_

i <- E$edge_type == "miRNA_target"
E$predicted_sign[i] <- -1L; E$sign_source[i] <- "mechanism_miRNA_repression"; E$annotation[i] <- "miRNA_target"

tf_mode <- setNames(cov_tf$mode, key(cov_tf$source, cov_tf$target))
i <- E$edge_type == "TF_target"
md <- unname(tf_mode[key(E$source, E$target)[i]]); md[is.na(md) | md == ""] <- "no_TRRUST_record"
E$annotation[i] <- md
E$predicted_sign[i] <- unname(c(Activation = 1L, Repression = -1L)[md])
E$sign_source[i] <- ifelse(md %in% c("Activation", "Repression"), "TRRUST", "unannotated")

tfm_mode <- setNames(cov_tfm$mode, key(cov_tfm$source, cov_tfm$target))
i <- E$edge_type == "TF_miRNA"
md <- unname(tfm_mode[key(E$source, E$target)[i]]); md[is.na(md) | md == ""] <- "no_TransmiR_record"
E$annotation[i] <- md
E$predicted_sign[i] <- unname(c(Activation = 1L, Repression = -1L)[md])
E$sign_source[i] <- ifelse(md %in% c("Activation", "Repression"), "TransmiR", "unannotated")

i <- E$edge_type %in% c("gene_gene", "miRNA_miRNA")
E$annotation[i] <- paste0(E$edge_type[i], "_association"); E$sign_source[i] <- "no_directional_prediction"

tb <- table(E$edge_type, ifelse(is.na(E$predicted_sign), "unannotated", as.character(E$predicted_sign)))
say("PREDICTED SIGN ASSIGNMENT:"); cat(capture.output(print(tb)), sep = "\n")
cat(capture.output(print(tb)), sep = "\n", file = con, append = TRUE)

tk <- setNames(tier$tier, key(tier$source, tier$target))
E$evidence_tier <- unname(tk[key(E$source, E$target)])
E$evidence_tier[E$edge_type != "miRNA_target"] <- NA
say("evidence tier joined for ", sum(!is.na(E$evidence_tier)), " of ",
    sum(E$edge_type == "miRNA_target"), " miRNA_target edges")

## ------------------------------------------------- 4. EDGE CORRELATIONS
E$rho <- rho_of(E$source, E$target)
E$n_samples <- ifelse(is.na(E$rho), NA_integer_, N_PAIRED)
E$p <- ifelse(is.na(E$rho), NA_real_, p_of(E$rho, N_PAIRED))
E$fdr <- NA_real_; ok <- !is.na(E$p); E$fdr[ok] <- p.adjust(E$p[ok], method = "BH")
E$concordant <- ifelse(is.na(E$rho) | is.na(E$predicted_sign), NA, sign(E$rho) == E$predicted_sign)
say("edges with computable rho: ", sum(ok), " / ", nrow(E),
    "  (", nrow(E) - sum(ok), " dropped because a partner is unmeasured / zero-variance)")

write.csv(E[, c("source","target","edge_type","predicted_sign","rho","p","fdr","n_samples",
                "evidence_tier","concordant","annotation","sign_source","source_type","target_type")],
          file.path(REV, "results", "edge_correlation.csv"), row.names = FALSE)
say("WROTE results/edge_correlation.csv rows=", nrow(E))

## ------------------------------------------------- 5. MATCHED RANDOM-PAIR NULL
deg_out <- table(factor(edges$source, levels = nodes$name))
deg_in  <- table(factor(edges$target, levels = nodes$name))
deg_out <- setNames(as.integer(deg_out), names(deg_out))
deg_in  <- setNames(as.integer(deg_in),  names(deg_in))

pool <- c(gene_ok, mirna_ok)
mean_expr <- c(rowMeans(G[gene_ok, , drop = FALSE]), rowMeans(M[mirna_ok, , drop = FALSE]))
platform <- setNames(c(rep("gene", length(gene_ok)), rep("mirna", length(mirna_ok))), pool)
dec <- setNames(rep(NA_integer_, length(pool)), pool)
for (pl in c("gene", "mirna")) {
  nn <- pool[platform == pl]; v <- mean_expr[nn]
  br <- unique(quantile(v, probs = seq(0, 1, 0.1)))
  dec[nn] <- as.integer(cut(v, breaks = br, include.lowest = TRUE))
}
ntype <- node_type[pool]; names(ntype) <- pool
say("null pool ", length(pool), " nodes (", sum(platform == "gene"), " gene/TF, ",
    sum(platform == "mirna"), " miRNA); expression deciles within platform pool")

Eset <- new.env(hash = TRUE, parent = emptyenv())
for (r in seq_len(nrow(edges))) assign(key(edges$source[r], edges$target[r]), TRUE, envir = Eset)

cand_cache <- new.env(hash = TRUE, parent = emptyenv())
get_cand <- function(nd, role) {
  ck <- paste(nd, role, sep = "\r")
  if (exists(ck, envir = cand_cache, inherits = FALSE)) return(get(ck, envir = cand_cache))
  ty <- ntype[[nd]]; d0 <- dec[[nd]]
  dvec <- if (role == "out") deg_out else deg_in
  w <- 0L; cands <- character(0)
  repeat {
    cands <- setdiff(pool[ntype == ty & !is.na(dec) & abs(dec - d0) <= w], nd)
    if (length(cands) >= 8L || w >= 9L) break
    w <- w + 1L
  }
  if (length(cands) > 8L) {
    dd <- abs(log2(dvec[cands] + 1) - log2(dvec[[nd]] + 1))
    cands <- cands[order(dd)][seq_len(max(8L, ceiling(length(cands) / 2)))]
  }
  res <- list(cands = cands, w = w,
              degdiff = median(abs(log2(dvec[cands] + 1) - log2(dvec[[nd]] + 1))))
  assign(ck, res, envir = cand_cache); res
}

NPER <- 10L
er <- which(!is.na(E$rho))
say("generating ", NPER, " matched random pairs for each of ", length(er), " measurable edges")
capS <- vector("list", length(er)); capT <- vector("list", length(er)); capI <- vector("list", length(er))
winS <- integer(length(er)); winT <- integer(length(er))
ddS <- numeric(length(er)); ddT <- numeric(length(er))
for (j in seq_along(er)) {
  r <- er[j]
  cs <- get_cand(E$source[r], "out"); ct2 <- get_cand(E$target[r], "in")
  winS[j] <- cs$w; winT[j] <- ct2$w; ddS[j] <- cs$degdiff; ddT[j] <- ct2$degdiff
  got <- 0L; tries <- 0L; seen <- character(0); ss <- character(0); tt2 <- character(0)
  while (got < NPER && tries < 400L) {
    tries <- tries + 1L
    a <- cs$cands[sample.int(length(cs$cands), 1L)]
    b <- ct2$cands[sample.int(length(ct2$cands), 1L)]
    if (a == b) next
    kk <- key(a, b)
    if (kk %in% seen) next
    if (exists(kk, envir = Eset, inherits = FALSE)) next
    seen <- c(seen, kk); ss <- c(ss, a); tt2 <- c(tt2, b); got <- got + 1L
  }
  capS[[j]] <- ss; capT[[j]] <- tt2; capI[[j]] <- rep(r, length(ss))
  if (j %% 1000 == 0) say("  ... ", j, " / ", length(er))
}
NL <- data.frame(edge_row = unlist(capI), rnd_source = unlist(capS), rnd_target = unlist(capT),
                 stringsAsFactors = FALSE)
NL$real_source <- E$source[NL$edge_row]; NL$real_target <- E$target[NL$edge_row]
NL$edge_type <- E$edge_type[NL$edge_row]; NL$predicted_sign <- E$predicted_sign[NL$edge_row]
NL$evidence_tier <- E$evidence_tier[NL$edge_row]
NL$rho <- rho_of(NL$rnd_source, NL$rnd_target)
NL$p <- p_of(NL$rho, N_PAIRED)
NL$concordant <- ifelse(is.na(NL$predicted_sign), NA, sign(NL$rho) == NL$predicted_sign)
say("random pairs generated: ", nrow(NL), " (mean ", round(nrow(NL) / length(er), 3), " per edge)")
say("expression-decile match EXACT: source ", sum(winS == 0), "/", length(winS),
    " ; target ", sum(winT == 0), "/", length(winT))
say("median |log2(deg+1) difference| real node vs its candidate pool: source ",
    signif(median(ddS, na.rm = TRUE), 3), " ; target ", signif(median(ddT, na.rm = TRUE), 3))
write.csv(NL, file.path(REV, "results", "edge_correlation_null.csv"), row.names = FALSE)
say("WROTE results/edge_correlation_null.csv rows=", nrow(NL))

## ---- stratified permutation machinery -------------------------------------
null_mat_for <- function(rows, valcol) {
  K <- length(rows)
  mm <- matrix(NA_real_, nrow = K, ncol = NPER)
  pos <- match(NL$edge_row, rows); keep <- !is.na(pos)
  pos <- pos[keep]; v <- as.numeric(NL[[valcol]][keep])
  ordv <- ave(seq_along(pos), pos, FUN = seq_along)
  sel <- ordv <= NPER
  mm[cbind(pos[sel], ordv[sel])] <- v[sel]
  mm
}
pack <- function(real_vals, null_mat) {
  allm <- cbind(real_vals, null_mat); K <- nrow(allm)
  navail <- rowSums(!is.na(allm))
  P <- matrix(NA_real_, K, ncol(allm))
  for (ii in seq_len(K)) { v <- allm[ii, ][!is.na(allm[ii, ])]; if (length(v)) P[ii, seq_along(v)] <- v }
  list(P = P, navail = navail, K = K, TOT = sum(allm, na.rm = TRUE), CNT = sum(navail))
}
strat_perm <- function(real_vals, null_mat, B = 10000) {
  pk <- pack(real_vals, null_mat); K <- pk$K
  obs <- mean(real_vals, na.rm = TRUE) - mean(null_mat, na.rm = TRUE)
  ge <- 0L; le <- 0L; chunk <- 1000L
  done <- 0L
  while (done < B) {
    bb <- min(chunk, B - done)
    idx <- matrix(ceiling(runif(K * bb) * pk$navail), nrow = K, ncol = bb)
    lin <- (idx - 1L) * K + seq_len(K)
    prs <- colSums(matrix(pk$P[lin], nrow = K, ncol = bb))
    d <- prs / K - (pk$TOT - prs) / (pk$CNT - K)
    ge <- ge + sum(d >= obs); le <- le + sum(d <= obs); done <- done + bb
  }
  list(obs = obs, p_greater = (ge + 1) / (B + 1), p_less = (le + 1) / (B + 1),
       p_two = min(1, 2 * min((ge + 1) / (B + 1), (le + 1) / (B + 1))))
}
auc_fun <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  r <- rank(c(x, y)); (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

## ------------------------------------------------- 6. NULL COMPARISON
say("=== B) NULL COMPARISON ===")
ann <- E$annotation
strata <- list(
  "ALL_sign_annotated"    = which(!is.na(E$rho) & !is.na(E$predicted_sign)),
  "miRNA_target_ALL"      = which(!is.na(E$rho) & E$edge_type == "miRNA_target"),
  "miRNA_target_strong"   = which(!is.na(E$rho) & E$edge_type == "miRNA_target" & E$evidence_tier == "strong"),
  "miRNA_target_weak"     = which(!is.na(E$rho) & E$edge_type == "miRNA_target" & E$evidence_tier == "weak"),
  "miRNA_target_predonly" = which(!is.na(E$rho) & E$edge_type == "miRNA_target" & E$evidence_tier == "predicted_only"),
  "TF_target_Activation"  = which(!is.na(E$rho) & E$edge_type == "TF_target" & ann == "Activation"),
  "TF_target_Repression"  = which(!is.na(E$rho) & E$edge_type == "TF_target" & ann == "Repression"),
  "TF_target_unannotated" = which(!is.na(E$rho) & E$edge_type == "TF_target" & !ann %in% c("Activation", "Repression")),
  "TF_miRNA_Activation"   = which(!is.na(E$rho) & E$edge_type == "TF_miRNA" & ann == "Activation"),
  "TF_miRNA_Repression"   = which(!is.na(E$rho) & E$edge_type == "TF_miRNA" & ann == "Repression"),
  "TF_miRNA_unannotated"  = which(!is.na(E$rho) & E$edge_type == "TF_miRNA" & !ann %in% c("Activation", "Repression")),
  "miRNA_miRNA_cluster"   = which(!is.na(E$rho) & E$edge_type == "miRNA_miRNA"),
  "gene_gene"             = which(!is.na(E$rho) & E$edge_type == "gene_gene")
)
res <- list(); B_PERM <- 10000
for (nm in names(strata)) {
  rows <- strata[[nm]]; if (length(rows) < 2) next
  rr <- E$rho[rows]; nmat <- null_mat_for(rows, "rho")
  nn <- as.vector(nmat); nn <- nn[!is.na(nn)]
  pr <- strat_perm(rr, nmat, B_PERM)
  ps <- E$predicted_sign[rows]
  has_sign <- all(!is.na(ps))
  if (has_sign) {
    prs <- strat_perm(ps * rr, nmat * ps, B_PERM)
    sn <- as.vector(nmat * ps); sn <- sn[!is.na(sn)]
    a_sgn <- auc_fun(ps * rr, sn)
    msr <- mean(ps * rr); msn <- mean(sn)
  } else { prs <- list(obs = NA, p_greater = NA); a_sgn <- NA; msr <- NA; msn <- NA }
  res[[nm]] <- data.frame(stratum = nm, n_edges = length(rows), n_null_pairs = length(nn),
    mean_rho_real = mean(rr), mean_rho_null = mean(nn),
    median_rho_real = median(rr), median_rho_null = median(nn),
    delta_mean_rho = pr$obs, perm_p_two_sided = pr$p_two,
    perm_p_real_lower = pr$p_less, perm_p_real_higher = pr$p_greater,
    mean_signed_rho_real = msr, mean_signed_rho_null = msn,
    delta_signed_rho = prs$obs, perm_p_signed_greater = prs$p_greater,
    AUC_absrho = auc_fun(abs(rr), abs(nn)), AUC_signed_rho = a_sgn, stringsAsFactors = FALSE)
  say(sprintf("%-22s n=%5d mean rho real=%+.4f null=%+.4f delta=%+.4f permP2=%.4g AUC|rho|=%.3f AUCsgn=%s",
      nm, length(rows), mean(rr), mean(nn), pr$obs, pr$p_two, auc_fun(abs(rr), abs(nn)),
      ifelse(is.na(a_sgn), "NA", sprintf("%.3f", a_sgn))))
}
NULLTAB <- do.call(rbind, res)
write.csv(NULLTAB, file.path(REV, "results", "null_comparison_summary.csv"), row.names = FALSE)
say("WROTE results/null_comparison_summary.csv rows=", nrow(NULLTAB))

## ------------------------------------------------- 7. SIGN CONCORDANCE
say("=== C) SIGN-CONCORDANCE TEST (replaces the ANOVA) ===")
cres <- list()
for (nm in names(strata)) {
  rows <- strata[[nm]]
  rows <- rows[!is.na(E$concordant[rows])]
  if (length(rows) < 2) next
  cc <- as.numeric(E$concordant[rows])
  nmat <- null_mat_for(rows, "concordant")
  nn <- as.vector(nmat); nn <- nn[!is.na(nn)]
  p0 <- mean(nn); k <- sum(cc); n <- length(cc)
  bt  <- binom.test(k, n, p = p0,  alternative = "greater")
  bt2 <- binom.test(k, n, p = 0.5, alternative = "two.sided")
  pm <- strat_perm(cc, nmat, 10000)
  sig_conc <- sum(cc == 1 & E$fdr[rows] < 0.05, na.rm = TRUE)
  cres[[nm]] <- data.frame(stratum = nm, n_edges_tested = n, n_concordant = k,
    concordance_rate = k / n, null_concordance_rate = p0, n_null_pairs = length(nn),
    excess_over_null = k / n - p0, binom_p_vs_null = bt$p.value, binom_p_vs_0.5 = bt2$p.value,
    perm_p_vs_null = pm$p_greater,
    n_concordant_and_FDR05 = sig_conc, frac_concordant_and_FDR05 = sig_conc / n,
    stringsAsFactors = FALSE)
  say(sprintf("%-22s n=%5d conc=%.4f null=%.4f excess=%+.4f binomP=%.4g permP=%.4g conc&FDR05=%d (%.4f)",
      nm, n, k / n, p0, k / n - p0, bt$p.value, pm$p_greater, sig_conc, sig_conc / n))
}
CONC <- do.call(rbind, cres)
write.csv(CONC, file.path(REV, "results", "sign_concordance_summary.csv"), row.names = FALSE)
say("WROTE results/sign_concordance_summary.csv rows=", nrow(CONC))

saveRDS(list(E = E, NL = NL, samples = sm, gene_ok = gene_ok, mirna_ok = mirna_ok),
        file.path(REV, "data", "edge_corr_workspace.rds"))
say("=== 08 DONE ===")
close(con)
