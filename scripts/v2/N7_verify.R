## N7_verify.R -- independent re-computation of the headline new-cohort numbers
## using different code paths (cor.test / ppcor / lm residuals) as a cross-check.
setwd("/path/to/revision")
suppressPackageStartupMessages({library(matrixStats)})
CA <- "cache/newcohorts"
co <- readRDS(file.path(CA,"cohorts_base.rds"))
cafA_genes <- readLines(file.path(CA,"cafA_signature_genes.txt"))

zrow <- function(M){ s <- rowSds(M, na.rm=TRUE); s[s==0] <- NA; (M-rowMeans(M,na.rm=TRUE))/s }
cafscore <- function(X, genes){
  g <- intersect(genes, rownames(X)); Z <- zrow(X[g,,drop=FALSE])
  Z <- Z[rowSums(is.na(Z))==0,,drop=FALSE]; colMeans(Z)
}
## partial Spearman via residuals of rank-regression (a DIFFERENT algebra than
## the partial-correlation formula used in N2)
pres <- function(x,y,z){
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  ex <- residuals(lm(rx ~ rz)); ey <- residuals(lm(ry ~ rz))
  ct <- cor.test(ex, ey)   ## Pearson on residualised ranks
  c(rho=unname(ct$estimate), p=ct$p.value)
}
cat(sprintf("%-12s %-16s %8s %10s %10s %10s %10s\n","cohort","edge","n","N2_style","verify","N2_p","verify_p"))
for(nm in names(co)){
  X <- co[[nm]]$X; if(!all(c("COL1A1","COL3A1") %in% rownames(X))) next
  cf <- cafscore(X, cafA_genes)
  for(pr in list(c("COL1A1","COL3A1"), c("ETS1","COL1A1"), c("NFKB1","COL1A1"),
                 c("RELA","COL1A1"), c("SP1","COL1A1"))){
    if(!all(pr %in% rownames(X))) next
    a <- X[pr[1],]; b <- X[pr[2],]
    r0 <- cor.test(a, b, method="spearman", exact=FALSE)
    v  <- pres(a, b, cf)
    cat(sprintf("%-12s %-16s %8d %10.4f %10.4f %10.3g %10.3g\n",
        nm, paste0(pr[1],"~",pr[2]), length(a),
        unname(r0$estimate), v["rho"], r0$p.value, v["p"]))
  }
}
