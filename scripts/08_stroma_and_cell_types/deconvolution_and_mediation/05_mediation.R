## Task D(ii): formal causal-mediation analysis -- how much of TF -> COL1A1 is
## carried by CAF (stromal) content?  Non-parametric bootstrap CIs.
suppressPackageStartupMessages({library(mediation); library(matrixStats)})
setwd("/path/to/revision")
SIMS <- 1000L
cat("### mediation bootstrap sims =", SIMS, "(boot=TRUE, perc CI)\n")
gexp <- readRDS("data/brca_gene_expr.rds"); mexp <- readRDS("data/brca_mirna_expr_canonical.rds")
ph <- readRDS("data/brca_pheno.rds"); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
sG <- sort(intersect(colnames(gexp), tum)); cat("gene-assay primary tumours n =", length(sG), "\n")
sB <- sort(intersect(sG, colnames(mexp)));  cat("both-assay primary tumours n =", length(sB), "\n")
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors=FALSE)$name
gmt <- readLines(system.file("extdata","SI_geneset.gmt", package="estimate"))
strom <- strsplit(gmt[grep("StromalSignature",gmt)],"\t")[[1]][-c(1,2)]
sigA <- intersect(setdiff(setdiff(strom, grep("^COL",strom,value=TRUE)), nodes), rownames(gexp))
sigB <- intersect(setdiff(c("DCN","LUM","FAP","PDGFRB","THY1","POSTN"), nodes), rownames(gexp))
cat("CAF-A genes:", length(sigA), "  CAF-B genes:", paste(sigB, collapse=","), "\n")
stopifnot(!any(grepl("^COL", c(sigA,sigB))), !any(c(sigA,sigB) %in% nodes))
caf <- function(S, ss){ X <- gexp[S, ss, drop=FALSE]; colMeans((X-rowMeans(X))/rowSds(X)) }
ns <- function(v){ n <- length(v); qnorm((rank(v)-0.5)/n) }   # rank -> normal scores

run_med <- function(xname, yname, cafname, ss, xassay="gene"){
  X <- if (xassay=="miRNA") mexp[xname, ss] else gexp[xname, ss]
  Y <- gexp[yname, ss]
  Cc <- if (cafname=="A") caf(sigA, ss) else caf(sigB, ss)
  d <- data.frame(x=ns(X), y=ns(Y), m=ns(Cc))
  mm <- lm(m ~ x, data=d); om <- lm(y ~ x + m, data=d)
  set.seed(7)
  fit <- mediate(mm, om, treat="x", mediator="m", boot=TRUE, sims=SIMS, boot.ci.type="perc")
  s <- summary(fit)
  data.frame(x=xname, y=yname, caf=cafname, n=length(ss),
    rho_raw = cor(X, Y, method="spearman"),
    total   = fit$tau.coef,  total_lo=fit$tau.ci[1],  total_hi=fit$tau.ci[2],  total_p=fit$tau.p,
    ACME    = fit$d0,        acme_lo=fit$d0.ci[1],    acme_hi=fit$d0.ci[2],    acme_p=fit$d0.p,
    ADE     = fit$z0,        ade_lo=fit$z0.ci[1],     ade_hi=fit$z0.ci[2],     ade_p=fit$z0.p,
    prop_med= fit$n0,        prop_lo=fit$n0.ci[1],    prop_hi=fit$n0.ci[2],    prop_p=fit$n0.p,
    a_path  = coef(mm)["x"], b_path=coef(om)["m"], cprime=coef(om)["x"],
    stringsAsFactors=FALSE)
}
## bootstrap Sobel / percentile CI on a*b, independent implementation
boot_ab <- function(xname, yname, cafname, ss, xassay="gene", B=SIMS){
  X <- if (xassay=="miRNA") mexp[xname, ss] else gexp[xname, ss]
  Y <- gexp[yname, ss]; Cc <- if(cafname=="A") caf(sigA, ss) else caf(sigB, ss)
  d <- data.frame(x=ns(X), y=ns(Y), m=ns(Cc)); n <- nrow(d)
  set.seed(13); ab <- pm <- numeric(B)
  for(b in seq_len(B)){ i <- sample.int(n, n, replace=TRUE); dd <- d[i,]
    a <- coef(lm(m~x,dd))["x"]; o <- coef(lm(y~x+m,dd)); ab[b] <- a*o["m"]
    pm[b] <- (a*o["m"])/(a*o["m"]+o["x"]) }
  mm <- lm(m~x,d); om <- lm(y~x+m,d)
  a <- coef(mm)["x"]; b <- coef(om)["m"]; sa <- summary(mm)$coef["x",2]; sb <- summary(om)$coef["m",2]
  z <- a*b/sqrt(b^2*sa^2 + a^2*sb^2)
  data.frame(x=xname,y=yname,caf=cafname,n=n,B=B, ab=a*b,
    ab_lo=quantile(ab,.025), ab_hi=quantile(ab,.975),
    propmed_boot=median(pm), pm_lo=quantile(pm,.025), pm_hi=quantile(pm,.975),
    sobel_z=z, sobel_p=2*pnorm(-abs(z)), row.names=NULL)
}

jobs <- list(c("ETS1","COL1A1","A"), c("ETS1","COL1A1","B"),
             c("NFKB1","COL1A1","A"), c("NFKB1","COL1A1","B"),
             c("SP1","COL1A1","A"), c("RELA","COL1A1","A"))
res <- do.call(rbind, lapply(jobs, function(j) run_med(j[1],j[2],j[3], sG)))
## miRNA axes (need both assays)
res <- rbind(res,
  do.call(rbind, lapply(list(c("hsa-miR-29a","COL3A1","A"), c("hsa-miR-29a","COL1A1","A"),
                             c("hsa-miR-29b","COL3A1","A")),
                        function(j) run_med(j[1],j[2],j[3], sB, xassay="miRNA"))))
print(res[,c("x","y","caf","n","rho_raw","total","ACME","acme_lo","acme_hi","ADE","ade_lo","ade_hi",
             "prop_med","prop_lo","prop_hi","prop_p")], digits=3, row.names=FALSE)
sob <- do.call(rbind, lapply(list(c("ETS1","COL1A1","A"),c("ETS1","COL1A1","B"),
                                  c("NFKB1","COL1A1","A"),c("NFKB1","COL1A1","B")),
                             function(j) boot_ab(j[1],j[2],j[3], sG)))
cat("\n--- independent bootstrap Baron-Kenny / Sobel cross-check (B =",SIMS,") ---\n")
print(sob, digits=3, row.names=FALSE)
write.csv(res, "results/v2/verify_mediation.csv", row.names=FALSE)
write.csv(sob, "results/v2/v2_mediation_sobel.csv", row.names=FALSE)
cat("DONE 05\n")
