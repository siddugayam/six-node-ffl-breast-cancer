## ==========================================================================
## N9_gtex_sex.R -- GTEx "breast - mammary tissue" contains MALE donors.
## Re-run the (b)/(c) correlations in females only, males only and all, and
## append the rows to newcohorts_summary.csv. All values computed in this run.
## ==========================================================================
setwd("/path/to/revision")
suppressPackageStartupMessages(library(matrixStats))
OUT <- "results/v2"; CA <- "cache/newcohorts"
co <- readRDS(file.path(CA,"cohorts_base.rds"))
cafA_genes <- readLines(file.path(CA,"cafA_signature_genes.txt"))
X <- co$GTEx_breast$X; ph <- co$GTEx_breast$pheno
cat("GTEx donors by SEX (1=male, 2=female):\n"); print(table(ph$sex, useNA="ifany"))

zrow <- function(M){ s <- rowSds(M, na.rm=TRUE); s[s==0] <- NA; (M-rowMeans(M,na.rm=TRUE))/s }
rk <- function(v){ r <- rank(v); (r-mean(r))/sd(r) }
sp <- function(x,y,z=NULL){
  n <- length(x); a <- rk(x); b <- rk(y)
  if(is.null(z)){ r <- cor(a,b); k <- 0 } else {
    c0 <- rk(z); rxy<-cor(a,b); rxz<-cor(a,c0); ryz<-cor(b,c0)
    r <- (rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2)); k <- 1 }
  df <- n-2-k; tt <- r*sqrt(df/(1-r^2))
  c(rho=r, n=n, p=2*pt(-abs(tt), df), se=1/sqrt(n-3-k))
}
strata <- list(all=rep(TRUE,ncol(X)), female=ph$sex==2, male=ph$sex==1)
rows <- list()
for(st in names(strata)){
  keep <- strata[[st]]
  Xs <- X[, keep, drop=FALSE]
  g <- intersect(cafA_genes, rownames(Xs)); Z <- zrow(Xs[g,,drop=FALSE])
  Z <- Z[rowSums(is.na(Z))==0,,drop=FALSE]; caf <- colMeans(Z)
  cat("\n--", st, " n =", ncol(Xs), " CAF-A genes =", nrow(Z), "--\n")
  for(pr in list(c("COL1A1","COL3A1"),c("ETS1","COL1A1"),c("NFKB1","COL1A1"),
                 c("RELA","COL1A1"),c("SP1","COL1A1"))){
    r0 <- sp(Xs[pr[1],], Xs[pr[2],]); rA <- sp(Xs[pr[1],], Xs[pr[2],], caf)
    cat(sprintf("  %-6s ~ %-6s  raw=%+.3f (p=%.3g) | CAF-A=%+.3f (p=%.3g)\n",
                pr[1], pr[2], r0["rho"], r0["p"], rA["rho"], rA["p"]))
    rows[[length(rows)+1]] <- data.frame(
      cohort=paste0("GTEx_breast_",st), accession="GTEx v8 breast mammary tissue",
      platform="Illumina TrueSeq RNA-seq (log2 TPM+1)", n=ncol(Xs),
      analysis=paste0("f_gtex_sexstratified_", gsub("[^A-Za-z0-9]","_",paste0(pr[1],"_",pr[2]))),
      n_tested=nrow(Z), n_concordant=NA_integer_,
      pct=unname(r0["rho"]), binom_p=unname(r0["p"]),
      extra=sprintf("rho_raw=%+.4f p_raw=%.3g; rho_CAFA=%+.4f p_CAFA=%.3g; se_raw=%.4f se_cafA=%.4f",
                    r0["rho"], r0["p"], rA["rho"], rA["p"], r0["se"], rA["se"]),
      stringsAsFactors=FALSE)
  }
}
S <- read.csv(file.path(OUT,"newcohorts_summary.csv"), stringsAsFactors=FALSE)
S <- S[!grepl("^f_gtex_sexstratified", S$analysis), ]
S2 <- rbind(S, do.call(rbind, rows))
write.csv(S2, file.path(OUT,"newcohorts_summary.csv"), row.names=FALSE)
cat("\nappended", length(rows), "sex-stratified GTEx rows ->", nrow(S2), "summary rows\n")
