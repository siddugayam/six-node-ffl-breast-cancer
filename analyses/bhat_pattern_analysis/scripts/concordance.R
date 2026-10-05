## Section X of SETTINGS.md: sign concordance of the twelve Bhat networks' edges in TCGA-BRCA and CPTAC, against matched
## random-pair nulls, with the methods of scripts/07_expression_validation/edge_correlations/08_expression_validation.R (TCGA) and scripts/07_expression_validation/protein_cptac_and_rppa/21_cptac_protein_validation.R
## (CPTAC).
##   Analysed-network edges use the paper's stored per-edge results and stored null pairs.
##   Census-layer TRRUST arcs have no stored result; they are computed here with the same method (08 / 21 code) and
##   reported separately.
##   Machinery checks: the whole-network strata recomputed from the stored files must equal the paper's summaries.
## Output: concordance/X1_tcga.csv, X1_trrust_layer_edges.csv, X1_trrust_layer_null_pairs.csv, X2_cptac.csv,
##         X2_trrust_layer_edges.csv, X_checks.txt, concordance.log (edge- and cohort-level values only).
## usage: Rscript concordance.R <analysis root> <network deposit folder> <original SIF folder> <output folder> [variant networks]
suppressPackageStartupMessages({ library(data.table) })
a <- commandArgs(TRUE); REV <- a[1]; NET <- a[2]; SIF <- a[3]; OUT <- a[4]; NSRC <- if (length(a) >= 5) a[5] else NET
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(system2("python3", c(file.path(HERE, "bhat_common.py"), "check", REV, NET, SIF)) == 0)
LOG <- file(file.path(OUT, "concordance.log"), "w")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG); flush(LOG) }
set.seed(1234)
NETS <- as.vector(outer(c("miRNA_FFL", "TF_FFL", "composite_FFL"), 3:6, function(k, n) paste0(n, "node_", k)))
key <- function(s, t) paste(s, t, sep = "\r")
ED <- lapply(setNames(NETS, NETS), function(n) fread(file.path(NSRC, paste0(n, "_edges.tsv")), colClasses = list(character = c("sign", "source", "target", "layer", "edge_type"))))
edges <- read.delim(file.path(REV, "data", "canonical_edges.tsv"), stringsAsFactors = FALSE)
nodes <- read.delim(file.path(REV, "data", "canonical_nodes.tsv"), stringsAsFactors = FALSE)
chk <- c()

## ---------------- stratified permutation machinery (verbatim from 08)
NPER <- 10L
null_mat_for_NL <- function(NL, rows, valcol) {
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
  list(obs = obs, p_greater = (ge + 1) / (B + 1), p_less = (le + 1) / (B + 1))
}
## one stratum of 08 section C; E has columns rho, s (predicted sign), conc; NL has edge_row, rho, and conc for that sign
conc_test <- function(E, NL, rows) {
  rows <- rows[!is.na(E$conc[rows])]
  if (length(rows) < 2) return(NULL)
  cc <- as.numeric(E$conc[rows]); nmat <- null_mat_for_NL(NL, rows, "conc")
  nn <- as.vector(nmat); nn <- nn[!is.na(nn)]
  p0 <- mean(nn); k <- sum(cc); n <- length(cc)
  pm <- strat_perm(cc, nmat, 10000)
  data.table(n_edges_tested = n, n_concordant = k, concordance_rate = k / n, null_concordance_rate = p0, n_null_pairs = length(nn),
             excess_over_null = k / n - p0, binom_p_vs_null = binom.test(k, n, p = p0, alternative = "greater")$p.value,
             binom_p_vs_0.5 = binom.test(k, n, p = 0.5)$p.value, perm_p_vs_null = pm$p_greater)
}
strata_of <- function(E, rows) {
  list(ALL_sign_annotated = rows,
       miRNA_target_ALL = rows[E$edge_type[rows] == "miRNA_target"],
       miRNA_target_strong = rows[E$edge_type[rows] == "miRNA_target" & E$tier[rows] %in% "strong"],
       miRNA_target_weak = rows[E$edge_type[rows] == "miRNA_target" & E$tier[rows] %in% "weak"],
       miRNA_target_predonly = rows[E$edge_type[rows] == "miRNA_target" & E$tier[rows] %in% "predicted_only"],
       TF_target_Activation = rows[E$edge_type[rows] == "TF_target" & E$s[rows] == 1],
       TF_target_Repression = rows[E$edge_type[rows] == "TF_target" & E$s[rows] == -1],
       TF_miRNA_Activation = rows[E$edge_type[rows] == "TF_miRNA" & E$s[rows] == 1],
       TF_miRNA_Repression = rows[E$edge_type[rows] == "TF_miRNA" & E$s[rows] == -1])
}

## ================= X1 TCGA =================
EC <- fread(file.path(REV, "results", "edge_correlation.csv"))
stopifnot(nrow(EC) == nrow(edges), identical(EC$source, edges$source), identical(EC$target, edges$target))
NL <- fread(file.path(REV, "results", "edge_correlation_null.csv"))
ecrow <- setNames(seq_len(nrow(EC)), key(EC$source, EC$target))
## the network files' signs against the stored predicted signs
an <- unique(rbindlist(lapply(ED, function(d) d[layer == "analysed network", .(source, target, edge_type, sign, evidence_tier)])))
an[, row := ecrow[key(source, target)]]
stopifnot(!anyNA(an$row))
an[, s_file := suppressWarnings(as.integer(sign))]
an[, pred := EC$predicted_sign[row]]
dis <- an[xor(is.na(s_file), is.na(pred)) | (!is.na(s_file) & !is.na(pred) & s_file != pred)]
chk <- c(chk, sprintf("CHECK X1 signs: %d distinct analysed-network edges in the twelve networks; sign in the network files differs from results/edge_correlation.csv predicted_sign for %d",
                      nrow(an), nrow(dis)))
## the per-edge table used below: all canonical edges, with the network sign where the files give one
E1 <- data.table(source = EC$source, target = EC$target, edge_type = EC$edge_type, rho = EC$rho, tier = EC$evidence_tier,
                 s = EC$predicted_sign)
E1[an$row, s := an$s_file]
E1[, conc := ifelse(is.na(rho) | is.na(s), NA, sign(rho) == s)]
NL1 <- data.table(edge_row = NL$edge_row, rho = NL$rho)
NL1[, conc := ifelse(is.na(E1$s[edge_row]), NA, sign(rho) == E1$s[edge_row])]
## machinery check: the whole network, 08's strata, against results/sign_concordance_summary.csv
PS <- fread(file.path(REV, "results", "sign_concordance_summary.csv"))
allrows <- which(!is.na(E1$rho) & !is.na(E1$s))
W <- rbindlist(lapply(names(strata_of(E1, allrows)), function(nm) {
  r <- conc_test(E1, NL1, strata_of(E1, allrows)[[nm]]); if (is.null(r)) NULL else cbind(stratum = nm, r) }))
mw <- merge(W, PS, by = "stratum", suffixes = c("", ".paper"))
chk <- c(chk, sprintf("CHECK X1 whole network: %d strata matched; n equal %s; concordance equal %s; max |null rate difference| %.3g; max |binomial p difference| %.3g",
                      nrow(mw), all(mw$n_edges_tested == mw$n_edges_tested.paper), all(mw$n_concordant == mw$n_concordant.paper),
                      max(abs(mw$null_concordance_rate - mw$null_concordance_rate.paper)), max(abs(mw$binom_p_vs_null - mw$binom_p_vs_null.paper))))
X1 <- list()
for (net in NETS) {
  rows <- ecrow[key(ED[[net]][layer == "analysed network"]$source, ED[[net]][layer == "analysed network"]$target)]
  rows <- rows[!is.na(E1$rho[rows]) & !is.na(E1$s[rows])]
  st <- strata_of(E1, rows)
  for (nm in names(st)) { r <- conc_test(E1, NL1, st[[nm]]); if (!is.null(r)) X1[[length(X1) + 1]] <- cbind(network = net, layer = "analysed network", stratum = nm, r) }
}
## ---- new TRRUST-layer arcs: 08's correlation and null, computed here
TL <- unique(rbindlist(lapply(ED, function(d) d[layer == "census layer: TRRUST TF-target", .(source, target, edge_type, sign)])))
say("TRRUST-layer arcs in the networks: ", nrow(TL), " (signed ", sum(TL$sign %in% c("1", "-1")), ")")
G0 <- readRDS(file.path(REV, "data", "brca_gene_expr.rds")); M0 <- readRDS(file.path(REV, "data", "brca_mirna_expr_canonical.rds"))
ph <- readRDS(file.path(REV, "data", "brca_pheno.rds"))
tumour <- ph$sample[ph$sample_type == "Primary Tumor"]
sm <- sort(intersect(intersect(colnames(G0), tumour), intersect(colnames(M0), tumour))); N_PAIRED <- length(sm)
G <- G0[, sm, drop = FALSE]; M <- M0[, sm, drop = FALSE]; rm(G0, M0); invisible(gc())
sdG <- apply(G, 1, sd); sdM <- apply(M, 1, sd)
node_type <- setNames(nodes$type, nodes$name)
gene_ok <- intersect(nodes$name[nodes$type %in% c("TF", "Gene")], rownames(G)[sdG > 0])
mirna_ok <- intersect(nodes$name[nodes$type == "miRNA"], rownames(M)[sdM > 0])
zrow <- function(X) { mu <- rowMeans(X); s <- sqrt(rowSums((X - mu)^2)); (X - mu) / s }
Z <- rbind(zrow(t(apply(G[gene_ok, , drop = FALSE], 1, rank))), zrow(t(apply(M[mirna_ok, , drop = FALSE], 1, rank))))
rho_of <- function(a, b) {
  ok <- a %in% rownames(Z) & b %in% rownames(Z); out <- rep(NA_real_, length(a))
  if (any(ok)) out[ok] <- rowSums(Z[a[ok], , drop = FALSE] * Z[b[ok], , drop = FALSE])
  pmax(pmin(out, 1), -1)
}
p_of <- function(rho, n) { tt <- rho * sqrt((n - 2) / pmax(1 - rho^2, .Machine$double.xmin)); 2 * pt(-abs(tt), df = n - 2) }
chk <- c(chk, sprintf("CHECK X1 recomputed rho: largest |difference| from results/edge_correlation.csv over its %d rhos: %.3g",
                      sum(!is.na(EC$rho)), max(abs(rho_of(EC$source, EC$target) - EC$rho), na.rm = TRUE)))
TL[, rho := rho_of(source, target)]; TL[, p := p_of(rho, N_PAIRED)]
TL[, s := suppressWarnings(as.integer(sign))]; TL[, conc := ifelse(is.na(rho) | is.na(s), NA, sign(rho) == s)]
## 08's matched null for the new arcs
deg_out <- table(factor(edges$source, levels = nodes$name)); deg_out <- setNames(as.integer(deg_out), names(deg_out))
deg_in <- table(factor(edges$target, levels = nodes$name)); deg_in <- setNames(as.integer(deg_in), names(deg_in))
pool <- c(gene_ok, mirna_ok)
mean_expr <- c(rowMeans(G[gene_ok, , drop = FALSE]), rowMeans(M[mirna_ok, , drop = FALSE]))
platform <- setNames(c(rep("gene", length(gene_ok)), rep("mirna", length(mirna_ok))), pool)
dec <- setNames(rep(NA_integer_, length(pool)), pool)
for (pl in c("gene", "mirna")) { nn <- pool[platform == pl]; v <- mean_expr[nn]; br <- unique(quantile(v, probs = seq(0, 1, 0.1)))
  dec[nn] <- as.integer(cut(v, breaks = br, include.lowest = TRUE)) }
ntype <- node_type[pool]; names(ntype) <- pool
Eset <- new.env(hash = TRUE, parent = emptyenv())
for (r in seq_len(nrow(edges))) assign(key(edges$source[r], edges$target[r]), TRUE, envir = Eset)
cand_cache <- new.env(hash = TRUE, parent = emptyenv())
get_cand <- function(nd, role) {
  ck <- paste(nd, role, sep = "\r")
  if (exists(ck, envir = cand_cache, inherits = FALSE)) return(get(ck, envir = cand_cache))
  ty <- ntype[[nd]]; d0 <- dec[[nd]]; dvec <- if (role == "out") deg_out else deg_in
  w <- 0L; cands <- character(0)
  repeat { cands <- setdiff(pool[ntype == ty & !is.na(dec) & abs(dec - d0) <= w], nd); if (length(cands) >= 8L || w >= 9L) break; w <- w + 1L }
  if (length(cands) > 8L) { dd <- abs(log2(dvec[cands] + 1) - log2(dvec[[nd]] + 1)); cands <- cands[order(dd)][seq_len(max(8L, ceiling(length(cands) / 2)))] }
  res <- list(cands = cands, w = w); assign(ck, res, envir = cand_cache); res
}
set.seed(1234)
er <- which(!is.na(TL$rho))
capS <- vector("list", length(er)); capT <- vector("list", length(er)); capI <- vector("list", length(er))
for (j in seq_along(er)) {
  r <- er[j]; cs <- get_cand(TL$source[r], "out"); ct2 <- get_cand(TL$target[r], "in")
  got <- 0L; tries <- 0L; seen <- character(0); ss <- character(0); tt2 <- character(0)
  while (got < NPER && tries < 400L) {
    tries <- tries + 1L
    aa <- cs$cands[sample.int(length(cs$cands), 1L)]; bb <- ct2$cands[sample.int(length(ct2$cands), 1L)]
    if (aa == bb) next
    kk <- key(aa, bb); if (kk %in% seen) next
    if (exists(kk, envir = Eset, inherits = FALSE)) next
    seen <- c(seen, kk); ss <- c(ss, aa); tt2 <- c(tt2, bb); got <- got + 1L
  }
  capS[[j]] <- ss; capT[[j]] <- tt2; capI[[j]] <- rep(r, length(ss))
}
NLt <- data.table(edge_row = unlist(capI), rnd_source = unlist(capS), rnd_target = unlist(capT))
NLt[, rho := rho_of(rnd_source, rnd_target)]
NLt[, conc := ifelse(is.na(TL$s[edge_row]), NA, sign(rho) == TL$s[edge_row])]
fwrite(TL[, .(source, target, edge_type, sign, rho, p, n_samples = ifelse(is.na(rho), NA, N_PAIRED), concordant = conc)], file.path(OUT, "X1_trrust_layer_edges.csv"))
fwrite(NLt[, .(edge_row, rnd_source, rnd_target, rho, concordant = conc)], file.path(OUT, "X1_trrust_layer_null_pairs.csv"))
tlrow <- setNames(seq_len(nrow(TL)), key(TL$source, TL$target))
for (net in NETS) {
  d <- ED[[net]][layer == "census layer: TRRUST TF-target"]
  if (!nrow(d)) next
  rows <- tlrow[key(d$source, d$target)]; rows <- rows[!is.na(TL$rho[rows]) & !is.na(TL$s[rows])]
  for (nm in c("TF_target_Activation", "TF_target_Repression")) {
    rr <- rows[TL$s[rows] == if (nm == "TF_target_Activation") 1 else -1]
    r <- conc_test(TL, NLt, rr); if (!is.null(r)) X1[[length(X1) + 1]] <- cbind(network = net, layer = "TRRUST layer (computed here)", stratum = nm, r)
  }
}
X1 <- rbindlist(X1); fwrite(X1, file.path(OUT, "X1_tcga.csv"))
say("X1 rows: ", nrow(X1))

## ================= X2 CPTAC =================
MIN_N <- 30L
EDGE <- fread(file.path(REV, "results", "multiomics", "cptac_miRNA_mRNA_vs_protein.csv"))
A <- EDGE[edge_type == "miRNA_target"]; B <- EDGE[edge_type == "TF_target"]
NC <- fread(file.path(REV, "results", "multiomics", "cptac_null_pairs.csv"))
auc_fun <- function(x, y) { n1 <- length(x); n2 <- length(y); if (n1 == 0 || n2 == 0) return(NA_real_); r <- rank(c(x, y)); (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2) }
summ21 <- function(EDF, rows, assay, NLl) {                 # 21's summarise(): concordance part
  rho_col <- paste0("rho_", assay); conc_col <- paste0("conc_", assay)
  rows <- rows[!is.na(EDF[[rho_col]][rows])]
  if (length(rows) < 3) return(NULL)
  rr <- EDF[[rho_col]][rows]
  nl <- NLl[NLl$edge_row %in% rows & !is.na(NLl$rho), ]
  ps <- EDF$predicted_sign[rows]
  nulconc <- if (all(!is.na(ps))) { psmap <- setNames(EDF$predicted_sign[rows], as.character(rows)); mean(sign(nl$rho) == psmap[as.character(nl$edge_row)], na.rm = TRUE) } else NA_real_
  cc <- EDF[[conc_col]][rows]; cc <- cc[!is.na(cc)]; k <- sum(cc); n <- length(cc)
  bt <- if (n > 0 && !is.na(nulconc) && nulconc > 0 && nulconc < 1) binom.test(k, n, p = nulconc, alternative = "greater")$p.value else NA_real_
  data.table(assay = assay, n_edges = length(rows), n_null_pairs = nrow(nl), mean_rho_real = mean(rr), mean_rho_null = mean(nl$rho),
             n_edges_conc_tested = n, n_concordant = k, concordance_rate = if (n > 0) k / n else NA_real_, null_concordance_rate = nulconc,
             excess_over_null = if (n > 0) k / n - nulconc else NA_real_, binom_p_vs_null = bt,
             binom_p_vs_0.5 = if (n > 0) binom.test(k, n, p = 0.5)$p.value else NA_real_)
}
strat21 <- function(Arows, Brows) list(
  miRNA_target_ALL = list("A", Arows), miRNA_target_strong = list("A", Arows[A$evidence_tier[Arows] %in% "strong"]),
  miRNA_target_weak = list("A", Arows[A$evidence_tier[Arows] %in% "weak"]),
  miRNA_target_predicted_only = list("A", Arows[A$evidence_tier[Arows] %in% "predicted_only"]),
  TF_target_Activation = list("B", Brows[B$trrust_mode[Brows] %in% "Activation"]),
  TF_target_Repression = list("B", Brows[B$trrust_mode[Brows] %in% "Repression"]),
  TF_target_ALL_signed = list("B", Brows[B$trrust_mode[Brows] %in% c("Activation", "Repression")]))
NLof <- list(A_mRNA = NC[layer == "A_miRNA_mRNA"], A_protein = NC[layer == "A_miRNA_protein"], B_mRNA = NC[layer == "B_TF_mRNA"], B_protein = NC[layer == "B_TF_protein"])
run21 <- function(Arows, Brows) {
  st <- strat21(Arows, Brows); out <- list()
  for (nm in names(st)) for (assay in c("mRNA", "protein")) {
    ab <- st[[nm]][[1]]; EDF <- if (ab == "A") A else B
    r <- summ21(EDF, st[[nm]][[2]], assay, NLof[[paste0(ab, "_", assay)]])
    if (!is.null(r)) out[[length(out) + 1]] <- cbind(stratum = nm, r)
  }
  rbindlist(out)
}
## machinery check: the whole network against results/multiomics/cptac_concordance_summary.csv
PC <- fread(file.path(REV, "results", "multiomics", "cptac_concordance_summary.csv"))
Wc <- run21(seq_len(nrow(A)), seq_len(nrow(B)))
mc <- merge(Wc, PC, by = c("stratum", "assay"), suffixes = c("", ".paper"))
chk <- c(chk, sprintf("CHECK X2 whole network: %d strata x assays matched; concordance rate max |difference| %.3g; null rate max |difference| %.3g; binomial p max |difference| %.3g",
                      nrow(mc), max(abs(mc$concordance_rate - mc$concordance_rate.paper)), max(abs(mc$null_concordance_rate - mc$null_concordance_rate.paper)),
                      max(abs(mc$binom_p_vs_null - mc$binom_p_vs_null.paper), na.rm = TRUE)))
## sign agreement of the stored CPTAC edges with the network files
arow <- setNames(seq_len(nrow(A)), key(A$source, A$target)); brow <- setNames(seq_len(nrow(B)), key(B$source, B$target))
X2 <- list()
for (net in NETS) {
  d <- ED[[net]][layer == "analysed network"]
  Ar <- arow[key(d$source, d$target)]; Ar <- sort(unname(Ar[!is.na(Ar)]))
  Br <- brow[key(d$source, d$target)]; Br <- sort(unname(Br[!is.na(Br)]))
  sB <- setNames(suppressWarnings(as.integer(d$sign)), key(d$source, d$target))
  bad <- sum(!is.na(B$predicted_sign[Br]) & B$predicted_sign[Br] != sB[key(B$source[Br], B$target[Br])], na.rm = TRUE) +
         sum(is.na(B$predicted_sign[Br]) != is.na(sB[key(B$source[Br], B$target[Br])]))
  if (bad) chk <- c(chk, sprintf("CHECK X2 signs %s: %d TF -> target edges whose CPTAC predicted sign differs from the network file", net, bad))
  r <- run21(Ar, Br); if (nrow(r)) X2[[length(X2) + 1]] <- cbind(network = net, layer = "analysed network", r)
}
chk <- c(chk, "CHECK X2 signs: every stored CPTAC TF -> target predicted sign equals the network files' sign unless listed above")
## ---- new TRRUST-layer arcs at CPTAC: 21's code
b <- readRDS(file.path(REV, "results", "multiomics", "cptac_bundle.rds"))
S3 <- b$s_all
Pm <- b$Pg[, S3, drop = FALSE]; Rm <- b$Rg[, S3, drop = FALSE]; Mm <- b$MIc[, S3, drop = FALSE]
nz <- function(X) apply(X, 1, function(v) { v <- v[!is.na(v)]; length(v) >= MIN_N && sd(v) > 0 })
prot_ok <- rownames(Pm)[nz(Pm)]; mrna_ok <- rownames(Rm)[nz(Rm)]
mir_all <- intersect(nodes$name[nodes$type == "miRNA"], rownames(Mm)); mir_ok <- mir_all[nz(Mm[mir_all, , drop = FALSE])]
gene_both <- sort(intersect(prot_ok, mrna_ok))
rho_pair <- function(X, a, Y, bnm) {
  n <- length(a); rr <- rep(NA_real_, n); nn <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    if (!(a[i] %in% rownames(X)) || !(bnm[i] %in% rownames(Y))) next
    x <- X[a[i], ]; y <- Y[bnm[i], ]; k <- !is.na(x) & !is.na(y)
    if (sum(k) < MIN_N) next
    if (sd(x[k]) == 0 || sd(y[k]) == 0) next
    rr[i] <- suppressWarnings(cor(x[k], y[k], method = "spearman")); nn[i] <- sum(k)
  }
  list(rho = rr, n = nn)
}
chk_r <- rho_pair(Rm, B$source[1:20], Rm, B$target[1:20])$rho
chk <- c(chk, sprintf("CHECK X2 recomputed CPTAC rho (first 20 TF -> target edges, mRNA): largest |difference| %.3g", max(abs(chk_r - B$rho_mRNA[1:20]), na.rm = TRUE)))
BT <- TL[source %in% gene_both & target %in% gene_both, .(source, target, edge_type, sign)]
BT[, predicted_sign := suppressWarnings(as.integer(sign))]
rm_ <- rho_pair(Rm, BT$source, Rm, BT$target); rp_ <- rho_pair(Pm, BT$source, Pm, BT$target)
BT[, `:=`(rho_mRNA = rm_$rho, n_mRNA = rm_$n, rho_protein = rp_$rho, n_protein = rp_$n)]
BT[, `:=`(conc_mRNA = ifelse(is.na(predicted_sign), NA, sign(rho_mRNA) == predicted_sign), conc_protein = ifelse(is.na(predicted_sign), NA, sign(rho_protein) == predicted_sign))]
say("CPTAC: TRRUST-layer arcs with both partners measured at mRNA and protein: ", nrow(BT))
if (nrow(BT)) {
  deg_out2 <- deg_out; deg_in2 <- deg_in
  pool2 <- c(gene_both, mir_ok); platform2 <- setNames(c(rep("gene", length(gene_both)), rep("mirna", length(mir_ok))), pool2)
  ntype2 <- node_type[pool2]; names(ntype2) <- pool2
  mk_dec <- function(v) { d <- setNames(rep(NA_integer_, length(pool2)), pool2)
    for (pl in c("gene", "mirna")) { nn <- pool2[platform2 == pl]; x <- v[nn]; br <- unique(quantile(x, probs = seq(0, 1, 0.1))); d[nn] <- as.integer(cut(x, breaks = br, include.lowest = TRUE)) }
    d }
  DEC <- list(mRNA = mk_dec(c(rowMeans(Rm[gene_both, , drop = FALSE], na.rm = TRUE), rowMeans(Mm[mir_ok, , drop = FALSE], na.rm = TRUE))),
              protein = mk_dec(c(rowMeans(Pm[gene_both, , drop = FALSE], na.rm = TRUE), rowMeans(Mm[mir_ok, , drop = FALSE], na.rm = TRUE))))
  make_null <- function(EDF, dec) {                     # 21's make_null
    cache <- new.env(hash = TRUE, parent = emptyenv())
    gc2 <- function(nd, role) {
      ck <- paste(nd, role, sep = "\r"); if (exists(ck, envir = cache, inherits = FALSE)) return(get(ck, envir = cache))
      ty <- ntype2[[nd]]; d0 <- dec[[nd]]; dvec <- if (role == "out") deg_out2 else deg_in2
      w <- 0L; cands <- character(0)
      repeat { cands <- setdiff(pool2[ntype2 == ty & !is.na(dec) & abs(dec - d0) <= w], nd); if (length(cands) >= 8L || w >= 9L) break; w <- w + 1L }
      if (length(cands) > 8L) { dd <- abs(log2(dvec[cands] + 1) - log2(dvec[[nd]] + 1)); cands <- cands[order(dd)][seq_len(max(8L, ceiling(length(cands) / 2)))] }
      res <- list(cands = cands, w = w); assign(ck, res, envir = cache); res
    }
    capS <- vector("list", nrow(EDF)); capT <- vector("list", nrow(EDF)); capI <- vector("list", nrow(EDF))
    for (j in seq_len(nrow(EDF))) {
      cs <- gc2(EDF$source[j], "out"); ct <- gc2(EDF$target[j], "in")
      got <- 0L; tries <- 0L; seen <- character(0); ss <- character(0); tt <- character(0)
      while (got < NPER && tries < 400L) { tries <- tries + 1L
        aa <- cs$cands[sample.int(length(cs$cands), 1L)]; bb <- ct$cands[sample.int(length(ct$cands), 1L)]
        if (aa == bb) next; kk <- key(aa, bb); if (kk %in% seen) next
        if (exists(kk, envir = Eset, inherits = FALSE)) next
        seen <- c(seen, kk); ss <- c(ss, aa); tt <- c(tt, bb); got <- got + 1L }
      capS[[j]] <- ss; capT[[j]] <- tt; capI[[j]] <- rep(j, length(ss))
    }
    data.frame(edge_row = unlist(capI), rnd_source = unlist(capS), rnd_target = unlist(capT), stringsAsFactors = FALSE)
  }
  set.seed(1234)
  nm_ <- make_null(BT, DEC$mRNA); np_ <- make_null(BT, DEC$protein)
  nm_$rho <- rho_pair(Rm, nm_$rnd_source, Rm, nm_$rnd_target)$rho; np_$rho <- rho_pair(Pm, np_$rnd_source, Pm, np_$rnd_target)$rho
  fwrite(BT[, .(source, target, edge_type, sign, rho_mRNA, n_mRNA, rho_protein, n_protein, conc_mRNA, conc_protein)], file.path(OUT, "X2_trrust_layer_edges.csv"))
  btrow <- setNames(seq_len(nrow(BT)), key(BT$source, BT$target))
  for (net in NETS) {
    d <- ED[[net]][layer == "census layer: TRRUST TF-target"]; rr <- btrow[key(d$source, d$target)]; rr <- sort(unname(rr[!is.na(rr)]))
    if (!length(rr)) next
    for (nm in c("TF_target_Activation", "TF_target_Repression", "TF_target_ALL_signed")) {
      sel <- rr[BT$predicted_sign[rr] %in% switch(nm, TF_target_Activation = 1L, TF_target_Repression = -1L, TF_target_ALL_signed = c(1L, -1L))]
      for (assay in c("mRNA", "protein")) {
        r <- summ21(BT, sel, assay, if (assay == "mRNA") nm_ else np_)
        if (!is.null(r)) X2[[length(X2) + 1]] <- cbind(network = net, layer = "TRRUST layer (computed here)", stratum = nm, r)
      }
    }
  }
}
X2 <- rbindlist(X2, fill = TRUE); fwrite(X2, file.path(OUT, "X2_cptac.csv"))
say("X2 rows: ", nrow(X2))
writeLines(chk, file.path(OUT, "X_checks.txt")); for (x in chk) say(x)
close(LOG)
