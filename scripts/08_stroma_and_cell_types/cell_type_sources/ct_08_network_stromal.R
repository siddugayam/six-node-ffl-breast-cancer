#!/usr/bin/env Rscript
# Network-wide: does the sign concordance of the edges survive adjustment for CAF content?
set.seed(42)
RES  <- "/path/to/revision/results/multiomics/"
DATA <- "/path/to/revision/data/"
CACHE<- "/path/to/revision/cache/celltype/"
expr <- readRDS(paste0(DATA,"brca_gene_expr.rds")); mir <- readRDS(paste0(DATA,"brca_mirna_expr_canonical.rds"))
sc <- readRDS(paste0(CACHE,"tcga_stromal_scores.rds")); rownames(sc) <- sc$sample
smp <- intersect(intersect(colnames(expr), colnames(mir)), sc$sample)
cat("samples:", length(smp), "\n")
G <- expr[, smp]; M <- mir[, smp]; z <- sc[smp, "score_CAF_scRNA"]
G <- G[apply(G,1,sd)>0, ]; M <- M[apply(M,1,sd)>0, ]
rz <- function(X){ R <- t(apply(X,1,rank)); t(scale(t(R))) }
Gz <- rz(G); Mz <- rz(M); zz <- as.numeric(scale(rank(z)))
resid_rows <- function(Xz, zz){ b <- as.numeric(Xz %*% zz)/sum(zz^2); R <- Xz - b %o% zz; t(scale(t(R))) }
Gr <- resid_rows(Gz, zz); Mr <- resid_rows(Mz, zz)
n <- length(smp)
pairrho <- function(A, B, s, t) rowSums(A[s,,drop=FALSE]*B[t,,drop=FALSE])/(n-1)

ec <- read.csv("/path/to/revision/results/edge_correlation.csv", stringsAsFactors=FALSE)
cat("edges in file:", nrow(ec), "\n")
mt <- ec[ec$edge_type=="miRNA_target" & !is.na(ec$predicted_sign) &
         ec$source %in% rownames(M) & ec$target %in% rownames(G), ]
tt <- ec[ec$edge_type=="TF_target" & ec$annotation %in% c("Activation","Repression") &
         ec$source %in% rownames(G) & ec$target %in% rownames(G), ]
cat("usable miRNA_target edges:", nrow(mt), " TF_target signed edges:", nrow(tt), "\n")
mt$rho_raw <- pairrho(Mz, Gz, mt$source, mt$target); mt$rho_adj <- pairrho(Mr, Gr, mt$source, mt$target)
tt$rho_raw <- pairrho(Gz, Gz, tt$source, tt$target); tt$rho_adj <- pairrho(Gr, Gr, tt$source, tt$target)
tt$predicted_sign <- ifelse(tt$annotation=="Activation", 1, -1)

# reuse the EXACT matched-random null pair set built in the transcriptomic phase
# (results/edge_correlation_null.csv, 10 pairs per edge) so the null is identical
nullf <- read.csv("/path/to/revision/results/edge_correlation_null.csv",
                  stringsAsFactors=FALSE)
cat("null pairs available:", nrow(nullf), "\n")
nullf$key <- paste(nullf$real_source, nullf$real_target, sep="|")
mt$key <- paste(mt$source, mt$target, sep="|"); tt$key <- paste(tt$source, tt$target, sep="|")
mtn <- nullf[nullf$key %in% mt$key & nullf$edge_type=="miRNA_target" &
             nullf$rnd_source %in% rownames(M) & nullf$rnd_target %in% rownames(G), ]
ttn <- nullf[nullf$key %in% tt$key & nullf$edge_type=="TF_target" &
             nullf$rnd_source %in% rownames(G) & nullf$rnd_target %in% rownames(G), ]
mtn$predicted_sign <- mt$predicted_sign[match(mtn$key, mt$key)]
ttn$predicted_sign <- tt$predicted_sign[match(ttn$key, tt$key)]
ttn$annotation <- tt$annotation[match(ttn$key, tt$key)]
mtn$evidence_tier <- mt$evidence_tier[match(mtn$key, mt$key)]
cat("null pairs used: miRNA_target", nrow(mtn), " TF_target", nrow(ttn), "\n")
mtn$rho_raw <- pairrho(Mz, Gz, mtn$rnd_source, mtn$rnd_target)
mtn$rho_adj <- pairrho(Mr, Gr, mtn$rnd_source, mtn$rnd_target)
ttn$rho_raw <- pairrho(Gz, Gz, ttn$rnd_source, ttn$rnd_target)
ttn$rho_adj <- pairrho(Gr, Gr, ttn$rnd_source, ttn$rnd_target)

conc <- function(d, col) mean(sign(d[[col]]) == d$predicted_sign, na.rm=TRUE)
mkrow <- function(lab, d, dn) {
  k_raw <- sum(sign(d$rho_raw)==d$predicted_sign); k_adj <- sum(sign(d$rho_adj)==d$predicted_sign)
  data.frame(stratum=lab, n_edges=nrow(d),
    concordance_unadjusted=100*conc(d,"rho_raw"), null_unadjusted=100*conc(dn,"rho_raw"),
    concordance_CAFadjusted=100*conc(d,"rho_adj"), null_CAFadjusted=100*conc(dn,"rho_adj"),
    binom_p_unadjusted=binom.test(k_raw, nrow(d), conc(dn,"rho_raw"), alternative="greater")$p.value,
    binom_p_CAFadjusted=binom.test(k_adj, nrow(d), conc(dn,"rho_adj"), alternative="greater")$p.value,
    mean_rho_unadjusted=mean(d$rho_raw), mean_rho_CAFadjusted=mean(d$rho_adj),
    stringsAsFactors=FALSE)
}
res <- rbind(
  mkrow("miRNA_target_ALL", mt, mtn),
  mkrow("miRNA_target_strong", mt[mt$evidence_tier %in% "strong",], mtn[mtn$evidence_tier %in% "strong",]),
  mkrow("TF_target_signed_ALL", tt, ttn),
  mkrow("TF_target_Activation", tt[tt$annotation=="Activation",], ttn[ttn$annotation=="Activation",]),
  mkrow("TF_target_Repression", tt[tt$annotation=="Repression",], ttn[ttn$annotation=="Repression",]))
res[,3:10] <- lapply(res[,3:10], function(x) signif(x,4))
write.csv(res, paste0(RES,"stromal_adjusted_network_concordance.csv"), row.names=FALSE)
print(res, row.names=FALSE)
cl <- c("source","target","edge_type","predicted_sign","annotation","evidence_tier","rho_raw","rho_adj")
alledges <- rbind(mt[,cl], tt[,cl])
write.csv(alledges, paste0(RES,"stromal_adjusted_edge_rho.csv"), row.names=FALSE)
cat("edge-level rows written:", nrow(alledges), "\n")
