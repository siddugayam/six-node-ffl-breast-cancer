#!/usr/bin/env Rscript
# C) Survival replication in METABRIC: hub set + collagen/TF module scores,
#    with and without adjustment for CAF (stromal) content.
suppressMessages({library(survival)})
set.seed(42)
RES  <- "/path/to/revision/results/multiomics/"
DATA <- "/path/to/revision/data/"
CACHE<- "/path/to/revision/cache/celltype/"
mb <- readRDS(paste0(DATA,"metabric.rds")); M <- mb$M; cl <- as.data.frame(mb$cl)
rownames(cl) <- cl$PATIENT_ID
cl <- cl[colnames(M), ]
cat("METABRIC expression:", dim(M), " clinical rows matched:", sum(!is.na(cl$PATIENT_ID)), "\n")
cl$OS_ev <- as.numeric(cl$OS_STATUS == "1:DECEASED"); cl$OS_t <- as.numeric(cl$OS_MONTHS)
cl$RFS_ev <- as.numeric(cl$RFS_STATUS == "1:Recurred"); cl$RFS_t <- as.numeric(cl$RFS_MONTHS)
cat("OS events:", sum(cl$OS_ev, na.rm=TRUE), " RFS events:", sum(cl$RFS_ev, na.rm=TRUE), "\n")

CAFSIG <- readLines(paste0(CACHE,"caf_signature_scrna.txt"))
gcaf <- intersect(CAFSIG, rownames(M))
caf <- colMeans(t(scale(t(M[gcaf,]))), na.rm=TRUE)
cat("CAF score built from", length(gcaf), "genes\n")

hubs <- read.csv("/path/to/revision/results/network_topology_hubs.csv", stringsAsFactors=FALSE)
hub_genes <- hubs$name[hubs$is_hub & hubs$type != "miRNA"]
named_ms  <- hubs$name[hubs$hub_named_in_MS & hubs$type != "miRNA"]
hub_genes <- intersect(hub_genes, rownames(M))
cat("hub genes tested in METABRIC:", length(hub_genes), " (of", sum(hubs$is_hub & hubs$type!="miRNA"), " protein-coding hubs)\n")

modules <- list(
  collagen_module = c("COL1A1","COL3A1"),
  TF_module       = c("NFKB1","RELA","SP1","ETS1"),
  collagen_TF_6gene_module = c("COL1A1","COL3A1","NFKB1","RELA","SP1","ETS1"),
  miR29_target_module = c("COL1A1","COL3A1","COL1A2","COL4A1","COL5A1","ELN","FBN1"),
  CAF_score = gcaf)
modscore <- sapply(modules, function(g){ g <- intersect(g, rownames(M)); colMeans(t(scale(t(M[g,,drop=FALSE]))), na.rm=TRUE) })
cat("module gene counts:", sapply(modules, function(g) length(intersect(g, rownames(M)))), "\n")

fitcox <- function(x, ev, tt, extra=NULL, label, endpoint, model) {
  d <- data.frame(x=as.numeric(scale(x)), ev=ev, tt=tt)
  if (!is.null(extra)) d <- cbind(d, extra)
  d <- d[complete.cases(d) & d$tt > 0, ]
  f <- if (is.null(extra)) Surv(tt, ev) ~ x else as.formula(paste("Surv(tt, ev) ~ x +", paste(setdiff(names(d), c("x","ev","tt")), collapse=" + ")))
  m <- try(coxph(f, data=d), silent=TRUE)
  if (inherits(m,"try-error")) return(NULL)
  s <- summary(m)
  data.frame(feature=label, endpoint=endpoint, model=model, n=nrow(d), n_event=sum(d$ev),
             HR_per_SD=s$coefficients["x","exp(coef)"],
             CI_low=s$conf.int["x","lower .95"], CI_high=s$conf.int["x","upper .95"],
             z=s$coefficients["x","z"], p=s$coefficients["x","Pr(>|z|)"],
             C_index=unname(s$concordance[1]), stringsAsFactors=FALSE)
}
out <- list()
for (endp in c("OS","RFS")) {
  ev <- cl[[paste0(endp,"_ev")]]; tt <- cl[[paste0(endp,"_t")]]
  for (g in hub_genes) out[[length(out)+1]] <- fitcox(M[g,], ev, tt, NULL, g, endp, "univariate")
  for (g in hub_genes) out[[length(out)+1]] <- fitcox(M[g,], ev, tt, data.frame(caf=as.numeric(scale(caf))), g, endp, "adjusted_for_CAF")
  for (mm in colnames(modscore)) {
    out[[length(out)+1]] <- fitcox(modscore[,mm], ev, tt, NULL, mm, endp, "univariate")
    if (mm != "CAF_score")
      out[[length(out)+1]] <- fitcox(modscore[,mm], ev, tt, data.frame(caf=as.numeric(scale(caf))), mm, endp, "adjusted_for_CAF")
    out[[length(out)+1]] <- fitcox(modscore[,mm], ev, tt,
      data.frame(age=as.numeric(cl$AGE_AT_DIAGNOSIS), grade=as.numeric(cl$GRADE),
                 size=as.numeric(cl$TUMOR_SIZE), nodes=as.numeric(cl$LYMPH_NODES_EXAMINED_POSITIVE)),
      mm, endp, "adjusted_for_clinical")
  }
}
res <- do.call(rbind, out)
res$q_BH <- NA
for (endp in unique(res$endpoint)) for (mo in unique(res$model)) {
  i <- res$endpoint==endp & res$model==mo & res$feature %in% hub_genes
  if (any(i)) res$q_BH[i] <- p.adjust(res$p[i], "BH")
}
write.csv(res, paste0(RES,"metabric_survival_replication.csv"), row.names=FALSE)
cat("\nmetabric_survival_replication.csv rows:", nrow(res), "\n")
for (endp in c("OS","RFS")) for (mo in c("univariate","adjusted_for_CAF")) {
  i <- res$endpoint==endp & res$model==mo & res$feature %in% hub_genes
  cat(sprintf("%s %-18s hubs tested %3d | nominal p<0.05: %3d | BH q<0.05: %3d\n",
      endp, mo, sum(i), sum(res$p[i]<0.05), sum(res$q_BH[i]<0.05, na.rm=TRUE)))
}
cat("\nHubs surviving BH q<0.05 (univariate):\n")
i <- res$model=="univariate" & res$feature %in% hub_genes & !is.na(res$q_BH) & res$q_BH<0.05
print(res[i, c("feature","endpoint","n","n_event","HR_per_SD","p","q_BH","C_index")], row.names=FALSE, digits=3)
cat("\nMODULE SCORES:\n")
print(res[res$feature %in% colnames(modscore), c("feature","endpoint","model","n","n_event","HR_per_SD","CI_low","CI_high","p","C_index")], row.names=FALSE, digits=3)
cat("\nNamed-in-manuscript hubs present in METABRIC:", paste(intersect(named_ms, rownames(M)), collapse=", "), "\n")
print(res[res$feature %in% named_ms & res$model=="univariate", c("feature","endpoint","HR_per_SD","p","q_BH")], row.names=FALSE, digits=3)
