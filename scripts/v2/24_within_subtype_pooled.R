## ============================================================================
## Addendum to TASK A: decompose each axis into its WITHIN-subtype and
## BETWEEN-subtype components. Adds a WITHIN_SUBTYPE_POOLED row per axis to
## subtype_stratified.csv. Motivated by miR-29c -> COL1A1, where the pooled
## correlation (-0.201) is WEAKER than every large subtype's (LumA -0.411).
## ============================================================================
setwd("/path/to/revision")
OUT <- "results/v2"

g  <- readRDS("data/brca_gene_expr.rds"); m <- readRDS("data/brca_mirna_expr_canonical.rds")
ph <- readRDS("data/brca_pheno.rds"); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
s  <- sort(intersect(intersect(colnames(g), colnames(m)), tum))
cl <- read.delim("/path/to/home/Desktop/DD/R_GPR/ESIA_BRCA/brca_eisa_pilot/clinical/BRCA_clinicalMatrix",
                 stringsAsFactors=FALSE)
p  <- setNames(cl$PAM50Call_RNAseq, cl$sampleID)[s]
k  <- !is.na(p) & p != ""; s <- s[k]
p  <- factor(p[k], levels=c("LumA","LumB","Her2","Basal","Normal"))
levels(p)[levels(p)=="Normal"] <- "Normal-like"
cat("n =", length(s), "\n"); print(table(p))

cafdef <- readRDS(file.path(OUT,"v2_caf_scores.rds")); sigA <- cafdef$sigA; sigB <- cafdef$sigB
zsc <- function(v)(v-mean(v))/sd(v)
mk  <- function(sig,X){gg<-intersect(sig,rownames(X)); colMeans(t(apply(X[gg,s,drop=FALSE],1,zsc)))}
cafA <- mk(sigA,g); cafB <- mk(sigB,g)

getv <- function(nm) if (nm %in% rownames(m)) m[nm,s] else g[nm,s]
pcor <- function(x,y,z){rx<-rank(x);ry<-rank(y);rz<-rank(z)
  a<-cor(rx,ry);b<-cor(rx,rz);c<-cor(ry,rz);(a-b*c)/sqrt((1-b^2)*(1-c^2))}

AX <- c(paste0(rep(c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"),each=2)," -> ",c("COL1A1","COL3A1")),
        paste0(rep(c("ETS1","NFKB1","SP1","RELA"),each=2)," -> ",c("COL1A1","COL3A1")),
        "COL1A1 ~ COL3A1")

out <- list()
cat("\n axis                    pooled   within-subtype   between-subtype(k=5 means)   ratio\n")
for (a in AX){
  pr <- strsplit(a, " ~ | -> ")[[1]]
  x <- getv(trimws(pr[1])); y <- getv(trimws(pr[2]))
  rx <- rank(x); ry <- rank(y)
  pooled  <- cor(rx, ry)
  cx <- rx - ave(rx, p); cy <- ry - ave(ry, p)     # remove subtype means from the ranks
  within  <- cor(cx, cy)
  betw    <- suppressWarnings(cor(tapply(x, p, mean), tapply(y, p, mean)))
  ## partial correlations of the within-subtype-centred variables vs CAF
  cz <- rank(cafA) - ave(rank(cafA), p)
  wpa <- { A<-cor(cx,cy); B<-cor(cx,cz); C<-cor(cy,cz); (A-B*C)/sqrt((1-B^2)*(1-C^2)) }
  czb <- rank(cafB) - ave(rank(cafB), p)
  wpb <- { A<-cor(cx,cy); B<-cor(cx,czb); C<-cor(cy,czb); (A-B*C)/sqrt((1-B^2)*(1-C^2)) }
  ## p-value for the within-subtype correlation: t test with df = n - k - 1
  df <- length(s) - nlevels(p) - 1
  tt <- within*sqrt(df/(1-within^2)); pv <- 2*pt(-abs(tt), df)
  cat(sprintf(" %-22s %+7.3f   %+7.3f          %+7.3f              %6.2f\n",
              a, pooled, within, betw, within/pooled))
  out[[length(out)+1]] <- data.frame(cohort="TCGA-BRCA", stratum="WITHIN_SUBTYPE_POOLED",
    n=length(s), axis=a, source=trimws(pr[1]), target=trimws(pr[2]),
    rho=within, p=pv, rho_partial_cafA=wpa, rho_partial_cafB=wpb,
    analysis="within-subtype-centred", stringsAsFactors=FALSE)
  out[[length(out)+1]] <- data.frame(cohort="TCGA-BRCA", stratum="BETWEEN_SUBTYPE_MEANS",
    n=nlevels(p), axis=a, source=trimws(pr[1]), target=trimws(pr[2]),
    rho=betw, p=NA, rho_partial_cafA=NA, rho_partial_cafB=NA,
    analysis="between-subtype (correlation of the 5 subtype means)", stringsAsFactors=FALSE)
}
ADD <- do.call(rbind, out)

ST <- read.csv(file.path(OUT,"subtype_stratified.csv"), stringsAsFactors=FALSE)
ST <- ST[ST$stratum != "WITHIN_SUBTYPE_POOLED" & ST$stratum != "BETWEEN_SUBTYPE_MEANS",]
for (cn in setdiff(colnames(ST), colnames(ADD))) ADD[[cn]] <- NA
ADD <- ADD[, colnames(ST)]
ST2 <- rbind(ST, ADD)
write.csv(ST2, file.path(OUT,"subtype_stratified.csv"), row.names=FALSE)
cat("\nsubtype_stratified.csv now has", nrow(ST2), "rows (added", nrow(ADD), ")\n")
cat("DONE ADDENDUM\n")
