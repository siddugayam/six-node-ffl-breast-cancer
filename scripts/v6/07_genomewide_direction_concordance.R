# Is the direction disagreement between HPA's breast prognostic calls and our own
# TCGA/METABRIC Cox specific to our 26 genes, or genome-wide? Test concordance
# separately among genes HPA itself calls prognostic (p<0.001) and unprognostic.
suppressPackageStartupMessages({library(survival)})
setwd("/path/to/revision")

pa <- read.delim("data/hpa/pathology.tsv", sep="\t", quote="", check.names=FALSE,
                 stringsAsFactors=FALSE, colClasses="character")
pa <- pa[pa$Cancer == "breast cancer", ]
dirn <- ifelse(pa[["prognostic - favorable"]]   != "", "favourable",
        ifelse(pa[["unprognostic - favorable"]] != "", "favourable",
        ifelse(pa[["prognostic - unfavorable"]] != "", "unfavourable",
        ifelse(pa[["unprognostic - unfavorable"]] != "", "unfavourable", NA))))
call <- ifelse(pa[["prognostic - favorable"]] != "" | pa[["prognostic - unfavorable"]] != "",
               "prognostic", ifelse(!is.na(dirn), "unprognostic", NA))
pval <- suppressWarnings(as.numeric(ifelse(pa[["prognostic - favorable"]]   != "", pa[["prognostic - favorable"]],
                 ifelse(pa[["unprognostic - favorable"]] != "", pa[["unprognostic - favorable"]],
                 ifelse(pa[["prognostic - unfavorable"]] != "", pa[["prognostic - unfavorable"]],
                        pa[["unprognostic - unfavorable"]])))))
hpa <- data.frame(gene=pa$`Gene name`, hpa_dir=dirn, hpa_call=call, hpa_p=pval,
                  stringsAsFactors=FALSE)
hpa <- hpa[!is.na(hpa$hpa_dir) & !duplicated(hpa$gene), ]
cat("HPA breast genes with a direction:", nrow(hpa),
    " prognostic:", sum(hpa$hpa_call=="prognostic"), "\n")

fastcox <- function(X, time, event){
  ok <- is.finite(time) & is.finite(event) & time > 0
  X <- X[, ok, drop=FALSE]; time <- time[ok]; event <- event[ok]
  z <- numeric(nrow(X)); p <- numeric(nrow(X))
  S <- Surv(time, event)
  for (i in seq_len(nrow(X))) {
    x <- as.numeric(X[i, ])
    if (length(unique(x)) < 5) { z[i] <- NA; p[i] <- NA; next }
    fit <- try(coxph(S ~ scale(x)), silent=TRUE)
    if (inherits(fit,"try-error")) { z[i] <- NA; p[i] <- NA; next }
    s <- summary(fit)$coefficients
    z[i] <- s[1,"z"]; p[i] <- s[1,"Pr(>|z|)"]
  }
  data.frame(gene=rownames(X), z=z, p=p, stringsAsFactors=FALSE)
}

out <- list()

ge <- readRDS("data/brca_gene_expr_symbol.rds")
sv <- readRDS("data/brca_survival.rds"); rownames(sv) <- sv$sample
cm <- intersect(colnames(ge), rownames(sv)); sv <- sv[cm, ]
ge <- ge[rownames(ge) %in% hpa$gene, cm, drop=FALSE]
ge <- ge[apply(ge, 1, function(r) sum(is.finite(r)) > 100), , drop=FALSE]
cat("TCGA genes tested:", nrow(ge), "samples:", ncol(ge), "\n")
r1 <- fastcox(ge, as.numeric(sv$OS.time), as.numeric(sv$OS)); r1$cohort <- "TCGA-BRCA"; r1$endpoint <- "OS"
out[[1]] <- r1

mb <- readRDS("data/metabric.rds"); M <- mb$M; cl <- as.data.frame(mb$cl)
rownames(cl) <- cl$SAMPLE_ID
cm2 <- intersect(colnames(M), rownames(cl)); cl <- cl[cm2, ]
M <- M[rownames(M) %in% hpa$gene, cm2, drop=FALSE]
cat("METABRIC genes tested:", nrow(M), "samples:", ncol(M), "\n")
r2 <- fastcox(M, as.numeric(cl$OS_time), as.numeric(cl$OS_event)); r2$cohort <- "METABRIC"; r2$endpoint <- "OS"
out[[2]] <- r2

res <- do.call(rbind, out)
res <- merge(res, hpa, by="gene")
res$our_dir <- ifelse(res$z > 0, "unfavourable", "favourable")
res$concordant <- res$our_dir == res$hpa_dir
write.csv(res, "results/v6/hpa_direction_concordance_genomewide.csv", row.names=FALSE)

tab <- do.call(rbind, lapply(split(res, list(res$cohort, res$hpa_call), drop=TRUE), function(d){
  d <- d[is.finite(d$z), ]
  bt <- binom.test(sum(d$concordant), nrow(d), 0.5)
  data.frame(cohort=d$cohort[1], hpa_call=d$hpa_call[1], n=nrow(d),
             n_concordant=sum(d$concordant), pct=round(100*mean(d$concordant),1),
             binom_p=bt$p.value)
}))
# additionally restrict to genes where OUR test is significant
tab2 <- do.call(rbind, lapply(split(res, list(res$cohort, res$hpa_call), drop=TRUE), function(d){
  d <- d[is.finite(d$z) & d$p < 0.01, ]
  if (nrow(d) < 10) return(NULL)
  bt <- binom.test(sum(d$concordant), nrow(d), 0.5)
  data.frame(cohort=d$cohort[1], hpa_call=paste0(d$hpa_call[1], " & our p<0.01"), n=nrow(d),
             n_concordant=sum(d$concordant), pct=round(100*mean(d$concordant),1),
             binom_p=bt$p.value)
}))
tab <- rbind(tab, tab2)
write.csv(tab, "results/v6/hpa_direction_concordance_summary.csv", row.names=FALSE)
print(tab)
