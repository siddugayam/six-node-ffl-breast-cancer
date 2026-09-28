#!/usr/bin/env Rscript
## Quantify the circularity in PART 2H: do miR-29 and miR-130a-3p have OPPOSITE bulk
## associations with stromal/CAF content in TCGA-BRCA?  If so, the apparent compartment
## segregation of their bulk-anti-correlated target sets is inherited from that, not from
## compartment-specific target repression.
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
mir <- readRDS(file.path(REV,"data/brca_mirna_expr_canonical.rds"))
ge  <- readRDS(file.path(REV,"data/brca_gene_expr.rds"))
msg("miRNA matrix:", nrow(mir), "x", ncol(mir), "| gene matrix:", nrow(ge), "x", ncol(ge))
stopifnot(sum(duplicated(rownames(ge)))==0)
tum <- function(x) x[substr(x,14,15)=="01"]
cs <- intersect(tum(colnames(mir)), tum(colnames(ge)))
msg("paired primary tumours:", length(cs))
D <- fread(file.path(OUT,"deconv_all_celltype_estimates.csv"))
msg("deconvolution methods:", paste(sort(unique(D$method)), collapse=", "))
msg("features available:", paste(sort(unique(D$feature))[1:20], collapse=", "))
## use every method that reports a fibroblast/CAF feature
fib <- D[grepl("CAF|ibroblast|stromal|Stromal", feature)]
msg("fibroblast-like features:", paste(unique(paste0(fib$method,":",fib$feature)), collapse=" | "))
rows <- list()
armnames <- grep("hsa-miR-130a|hsa-miR-29", rownames(mir), value=TRUE)
msg("miRNA rows matched:", paste(armnames, collapse=", "))
mir29 <- colMeans(mir[grep("hsa-miR-29[abc]$|hsa-miR-29[abc]-3p$", rownames(mir)), , drop=FALSE], na.rm=TRUE)
m130  <- if ("hsa-miR-130a-3p" %in% rownames(mir)) mir["hsa-miR-130a-3p",] else colMeans(mir[grep("hsa-miR-130a",rownames(mir)),,drop=FALSE],na.rm=TRUE)
msg("miR-29 feature built from:", paste(grep("hsa-miR-29[abc]$|hsa-miR-29[abc]-3p$", rownames(mir), value=TRUE), collapse=","))
for (k in unique(paste0(fib$method,"||",fib$feature))) {
  mm <- strsplit(k,"\\|\\|")[[1]]
  sub <- fib[method==mm[1] & feature==mm[2]]
  v <- setNames(sub$value, sub$sample)
  ov <- intersect(names(v), cs); if (length(ov) < 100) next
  for (mn in c("miR29","miR130a_3p")) {
    x <- if (mn=="miR29") mir29[ov] else m130[ov]
    ct <- suppressWarnings(cor.test(x, v[ov], method="spearman", exact=FALSE))
    rows[[length(rows)+1]] <- data.table(method=mm[1], feature=mm[2], mirna=mn, n=length(ov),
                                         rho=unname(ct$estimate), p=ct$p.value)
  }
}
R <- rbindlist(rows); setorder(R, mirna, -rho)
print(R, nrows=200)
fwrite(R, file.path(OUT,"sc_verify_mirna_vs_stroma_circularity.csv"))
msg("SUMMARY: miR-29 median rho vs fibroblast content =", round(median(R[mirna=="miR29", rho]),4),
    "| miR-130a-3p median rho =", round(median(R[mirna=="miR130a_3p", rho]),4))
msg("sign agreement: miR-29 negative in", sum(R[mirna=="miR29", rho]<0), "of", nrow(R[mirna=="miR29"]),
    "features; miR-130a-3p positive in", sum(R[mirna=="miR130a_3p", rho]>0), "of", nrow(R[mirna=="miR130a_3p"]))
