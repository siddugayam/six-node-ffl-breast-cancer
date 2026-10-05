## Section S of SETTINGS.md: each circuit of the four-, five- and six-node Bhat networks against its three-node core.
## A generalised copy of S8/s8_circuit_survival.R (analyses/six_node_pattern). Only its input changes: it reads the observed
## instance listings of section O (written by bhat_null.c, in the S1 program's order). Data, scores, models and
## correction are S8's.
##   circuits  the distinct node sets of each network; for a node set with several role assignments the first one
##             listed is used. Core = (TF1, miR1, G1); extension = the other roles (G2 at four nodes; miR2, G2 at five;
##             TF2, miR2, G2 at six). A one-gene extension (four nodes) is scored by the gene's z-scaled expression,
##             because GSVA cannot score a one-member set (minSize 2).
##   data      TCGA-BRCA primary tumours with gene + miRNA expression and survival; endpoints OS and PFI
##   PRIMARY   GSVA scores (kcdf Gaussian, minSize 2), z-scaled; Cox adjusted for age and AJCC stage (complete cases):
##             core ~ core + age + stage; circuit ~ circuit + age + stage; nested ~ core + extension + age + stage.
##             dC = C(circuit model) - C(core model) (Harrell); LRT core vs nested; BH across the circuits of each
##             network, within method x endpoint.
##   SECONDARY mean-z scores, univariable Cox, the same dC and LRT.
## GSVA scores are computed in chunks of gene sets in parallel (a gene set's GSVA score does not depend on the other
## sets scored with it); the Cox models are fitted in parallel over circuits.
## Output: circuits.csv (one row per circuit x endpoint x method; cohort-level n and events only), summary.csv,
## significant.csv, node_sharing.txt, check_s8.txt, survival.log.
## usage: Rscript survival.R <analysis root> <listing folder> <node names file> <S8 circuits file or "none"> <output folder> <cores>
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA); library(parallel) })
a <- commandArgs(TRUE)
REV <- a[1]; OBS <- a[2]; NAMES <- a[3]; S8 <- a[4]; OUT <- a[5]; NC <- as.integer(a[6])
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(OUT, "survival.log"), "w")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG); flush(LOG) }
set.seed(1234)
names587 <- readLines(NAMES)[1:587]
ROLES <- list(`4` = c("miR1", "TF1", "G1", "G2"), `5` = c("miR1", "miR2", "TF1", "G1", "G2"),
              `6` = c("miR1", "miR2", "TF1", "TF2", "G1", "G2"))
EXT <- list(`4` = "G2", `5` = c("miR2", "G2"), `6` = c("TF2", "miR2", "G2"))
NETS <- as.vector(outer(c("miRNA_FFL", "TF_FFL", "composite_FFL"), 4:6, function(k, n) paste0(n, "node_", k)))

first_per_set <- function(D, fam, roles) {
  D[, key := apply(.SD, 1, function(r) paste(sort(r), collapse = ";")), .SDcols = roles]
  D <- D[!duplicated(key)]; D[, family := fam]; D
}
CIR <- rbindlist(lapply(NETS, function(net) {
  n <- substr(net, 1, 1); k <- sub("^[0-9]node_", "", net)
  fn <- file.path(OBS, sprintf("observed_bhat%s_%s_instances.tsv", n, if (n == "6" && k == "composite_FFL") "composite" else k))
  D <- fread(fn, header = FALSE, col.names = ROLES[[n]]); D <- D[, lapply(.SD, function(i) names587[i + 1])]
  D <- first_per_set(D, net, ROLES[[n]])
  D[, core := paste(TF1, miR1, G1, sep = ";")]
  D[, ext := do.call(paste, c(.SD, sep = ";")), .SDcols = EXT[[n]]]
  D[, size := as.integer(n)]
  D[, members := key]
  D[, c("family", "size", "members", "core", "ext")]
}))
say("node sets: ", paste(CIR[, .N, by = family][, paste(family, N)], collapse = "; "))

gex <- readRDS(file.path(REV, "data/brca_gene_expr.rds")); mex <- readRDS(file.path(REV, "data/brca_mirna_expr_canonical.rds"))
pheno <- readRDS(file.path(REV, "data/brca_pheno.rds")); surv <- as.data.table(readRDS(file.path(REV, "data/brca_survival_clean.rds")))
tumour <- pheno$sample[pheno$sample_type == "Primary Tumor"]
common <- Reduce(intersect, list(colnames(gex), colnames(mex), tumour, surv$sample_id))
X <- rbind(gex[, common, drop = FALSE], mex[, common, drop = FALSE])
v <- apply(X, 1, function(z) stats::sd(z, na.rm = TRUE)); X <- X[is.finite(v) & v > 0, , drop = FALSE]
say("samples: ", length(common), "; features: ", nrow(X))
cl <- surv[sample_id %in% common]; setkey(cl, sample_id); cl <- cl[common]; stopifnot(identical(cl$sample_id, common))

have <- rownames(X)
CIR[, covered := vapply(strsplit(members, ";"), function(s) all(s %in% have), logical(1))]
say("circuits with all members measured: ", paste(CIR[, .(n = sum(covered), of = .N), by = family][, paste0(family, " ", n, " of ", of)], collapse = "; "))
CIR <- CIR[covered == TRUE]
setsl <- unique(c(CIR$members, CIR$core, CIR$ext))
gsl <- setNames(lapply(setsl, function(s) strsplit(s, ";")[[1]]), setsl)
multi <- gsl[lengths(gsl) >= 2]; single <- names(gsl)[lengths(gsl) == 1]
say("distinct gene sets: ", length(gsl), " (", length(multi), " scored by GSVA; ", length(single), " one-gene extensions)")
t0 <- Sys.time()
chunks <- split(names(multi), cut(seq_along(multi), NC, labels = FALSE))
GSl <- mclapply(chunks, function(ch) gsva(gsvaParam(X, multi[ch], kcdf = "Gaussian", minSize = 2, maxSize = Inf), verbose = FALSE),
                mc.cores = NC)
stopifnot(!any(vapply(GSl, inherits, logical(1), "try-error")))
GS <- do.call(rbind, GSl); stopifnot(setequal(rownames(GS), names(multi)))
say("GSVA done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min; ", paste(dim(GS), collapse = " x "))
Z <- t(scale(t(X)))
if (length(single)) GS <- rbind(GS, Z[single, , drop = FALSE])      # one-gene extension: the gene's z-scaled expression
MZ <- t(vapply(gsl, function(s) colMeans(Z[s, , drop = FALSE]), numeric(ncol(Z))))
sc <- function(M, s) as.numeric(scale(M[s, ]))

fit_one <- function(i, M, base, ok, adj) {
  d <- base; d$core <- sc(M, CIR$core[i])[ok]; d$ext <- sc(M, CIR$ext[i])[ok]; d$circ <- sc(M, CIR$members[i])[ok]
  if (adj) {
    f3 <- coxph(Surv(time, ev) ~ core + age + stage, data = d); f6 <- coxph(Surv(time, ev) ~ circ + age + stage, data = d)
    fn <- coxph(Surv(time, ev) ~ core + ext + age + stage, data = d)
  } else {
    f3 <- coxph(Surv(time, ev) ~ core, data = d); f6 <- coxph(Surv(time, ev) ~ circ, data = d)
    fn <- coxph(Surv(time, ev) ~ core + ext, data = d)
  }
  c3 <- summary(f3)$concordance[1]; c6 <- summary(f6)$concordance[1]
  data.table(family = CIR$family[i], size = CIR$size[i], members = CIR$members[i], core = CIR$core[i], extension = CIR$ext[i],
             n = nrow(d), events = sum(d$ev == 1), HR_circuit = unname(exp(coef(f6)["circ"])),
             p_circuit = summary(f6)$coefficients["circ", "Pr(>|z|)"], C_core = c3, C_circuit = c6, deltaC = c6 - c3,
             LRT_p_extension_adds = anova(f3, fn)[2, "Pr(>|Chi|)"])
}
res <- list()
for (ep in c("OS", "PFI")) {
  tt <- cl[[paste0(ep, ".time")]]; ee <- cl[[ep]]
  okb <- is.finite(tt) & tt > 0 & is.finite(ee)
  okc <- okb & !is.na(cl$age) & !is.na(cl$stage_group)
  say(ep, ": univariable n = ", sum(okb), " events ", sum(ee[okb] == 1), "; adjusted n = ", sum(okc), " events ", sum(ee[okc] == 1))
  for (meth in c("GSVA_adjusted_age_stage", "meanZ_univariable")) {
    adj <- meth == "GSVA_adjusted_age_stage"; M <- if (adj) GS else MZ; ok <- if (adj) okc else okb
    base <- data.frame(time = tt[ok], ev = ee[ok], age = cl$age[ok], stage = droplevels(cl$stage_group[ok]))
    t0 <- Sys.time()
    idx <- split(seq_len(nrow(CIR)), cut(seq_len(nrow(CIR)), NC * 4, labels = FALSE))
    R <- rbindlist(mclapply(idx, function(ii) rbindlist(lapply(ii, fit_one, M = M, base = base, ok = ok, adj = adj)), mc.cores = NC))
    R[, `:=`(method = meth, endpoint = ep)]
    res[[length(res) + 1]] <- R
    say(ep, " ", meth, ": ", nrow(R), " circuits in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
}
R <- rbindlist(res); setcolorder(R, c("method", "endpoint"))
R[, q_BH := p.adjust(LRT_p_extension_adds, "BH"), by = .(method, endpoint, family)]
fwrite(R, file.path(OUT, "circuits.csv"))
S <- R[, .(circuits = .N, median_deltaC = median(deltaC), pct_circuit_beats_core = 100 * mean(deltaC > 0),
           n_LRT_nominal_p05 = sum(LRT_p_extension_adds < 0.05), n_q_lt_005 = sum(q_BH < 0.05), min_q = min(q_BH)),
       by = .(method, endpoint, family)]
S[, family := factor(family, levels = NETS)]; setorder(S, method, endpoint, family)
fwrite(S, file.path(OUT, "summary.csv")); print(S)
G <- R[q_BH < 0.05][order(method, endpoint, family, LRT_p_extension_adds)]
fwrite(G, file.path(OUT, "significant.csv"))
ns <- c()
for (m in unique(R$method)) for (ep in c("OS", "PFI")) for (f in NETS) {
  g <- G[method == m & endpoint == ep & family == f]
  if (!nrow(g)) next
  tb <- sort(table(unlist(strsplit(g$members, ";"))), decreasing = TRUE)
  ns <- c(ns, sprintf("%s %s %s: %d circuits at q < 0.05; q %.3g-%.3g; deltaC %+.4f to %+.4f; nodes by the number of circuits containing them: %s",
                      m, ep, f, nrow(g), min(g$q_BH), max(g$q_BH), min(g$deltaC), max(g$deltaC),
                      paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = "; ")))
}
writeLines(if (length(ns)) ns else "no circuit at q < 0.05", file.path(OUT, "node_sharing.txt"))
## check against S8: the six-node composite network (main run; "none" for the variant networks of A2 and A3)
if (S8 == "none") { close(LOG); quit(save = "no") }
s8 <- fread(S8)[family == "BHAT6_composite"]
me <- R[family == "6node_composite_FFL"]
mm <- merge(me, s8, by = c("method", "endpoint", "members"), suffixes = c("", ".s8"))
num <- c(HR_circuit = "HR_six", p_circuit = "p_six", C_core = "C_core", C_circuit = "C_six", deltaC = "deltaC",
         LRT_p_extension_adds = "LRT_p_extension_adds", q_BH = "q_BH")
dif <- vapply(names(num), function(x) max(abs(mm[[x]] - mm[[paste0(num[[x]], if (num[[x]] %in% names(me)) ".s8" else "")]])), numeric(1))
chk <- c(sprintf("CHECK S8: rows here %d, rows in S8 %d, matched on method, endpoint and members %d; core and extension equal: %s; n and events equal: %s",
                 nrow(me), nrow(s8), nrow(mm), all(mm$core == mm$core.s8 & mm$extension == mm$extension.s8),
                 all(mm$n == mm$n.s8 & mm$events == mm$events.s8)),
         sprintf("CHECK S8: largest absolute difference: %s", paste(sprintf("%s %.3g", names(dif), dif), collapse = "; ")),
         sprintf("CHECK S8: GSVA OS circuits at q < 0.05 here %d, in S8 %d", nrow(me[method == "GSVA_adjusted_age_stage" & endpoint == "OS" & q_BH < 0.05]),
                 nrow(s8[method == "GSVA_adjusted_age_stage" & endpoint == "OS" & q_BH < 0.05])))
writeLines(chk, file.path(OUT, "check_s8.txt")); for (x in chk) say(x)
close(LOG)
