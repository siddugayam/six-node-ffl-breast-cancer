## ==========================================================================
## N8_inventory.R -- append an explicit per-cohort inventory to
## results/v2/newcohorts_summary.csv: accession, n, platform, genes measured,
## survival endpoints and event counts, covariates usable, and an explicit
## statement of which of analyses (a)-(e) could and could not be run.
## All counts computed in this run.
## ==========================================================================
setwd("/path/to/revision")
OUT <- "results/v2"; CA <- "cache/newcohorts"
co <- readRDS(file.path(CA,"cohorts_base.rds"))
S  <- read.csv(file.path(OUT,"newcohorts_summary.csv"), stringsAsFactors=FALSE)
S  <- S[!grepl("^z_cohort_inventory", S$analysis), ]   ## idempotent
E  <- read.csv(file.path(OUT,"newcohorts_per_edge.csv"), stringsAsFactors=FALSE)
hubs <- read.csv("results/network_topology_hubs.csv", stringsAsFactors=FALSE)
HUBS <- hubs$name[hubs$is_hub & hubs$type!="miRNA"]

rows <- list()
for(nm in names(co)){
  o <- co[[nm]]; ph <- o$pheno
  eps <- c()
  for(ep in c("OS","RFS","DFS","DMFS","DRFS")){
    tc <- paste0(ep,"_time"); ec <- paste0(ep,"_event")
    if(all(c(tc,ec) %in% names(ph))){
      ok <- is.finite(ph[[tc]]) & ph[[tc]]>0 & !is.na(ph[[ec]])
      if(sum(ok)>30 && sum(ph[[ec]][ok])>=10)
        eps <- c(eps, sprintf("%s(n=%d,events=%d)", ep, sum(ok), sum(ph[[ec]][ok])))
    }
  }
  cvs <- c()
  for(cv in c("age","grade","stage","ER","node","size"))
    if(cv %in% names(ph)){
      v <- ph[[cv]]
      if(mean(!is.na(v))>0.7 && length(unique(na.omit(v)))>1) cvs <- c(cvs, cv)
    }
  done <- c()
  if(any(E$cohort==nm & E$class=="TF_target")) done <- c(done,"a,b,c")
  if(any(E$cohort==nm & E$class=="cox_hub"))   done <- c(done,"d")
  if(nm != "GTEx_breast") done <- c(done,"e(meta)")
  cant <- c("miRNA layer (no miRNA assay in this cohort)")
  if(!length(eps)) cant <- c(cant, "(d) survival: no usable endpoint")
  if(nm=="GTEx_breast") cant <- c(cant, "non-tumour: excluded from the tumour meta-analysis pool")
  if(nm=="GSE96058") cant <- c(cant,
      "expression restricted to the 173 pre-selected symbols streamed from the 30,865-gene matrix")
  rows[[length(rows)+1]] <- data.frame(
    cohort=nm, accession=o$accession, platform=o$platform, n=ncol(o$X),
    analysis="z_cohort_inventory",
    n_tested=nrow(o$X), n_concordant=length(intersect(HUBS, rownames(o$X))),
    pct=NA_real_, binom_p=NA_real_,
    extra=sprintf("genes=%d; hubs_measurable=%d/50; endpoints=%s; covariates=%s; analyses_run=%s; NOT_possible=%s; %s",
      nrow(o$X), length(intersect(HUBS, rownames(o$X))),
      if(length(eps)) paste(eps, collapse=" ") else "NONE",
      if(length(cvs)) paste(cvs, collapse="+") else "NONE",
      paste(done, collapse=" "), paste(cant, collapse=" | "), o$notes),
    stringsAsFactors=FALSE)
}
INV <- do.call(rbind, rows)
print(INV[, c("cohort","accession","platform","n","n_concordant")], row.names=FALSE)
cat("\n")
for(i in seq_len(nrow(INV))) cat(INV$cohort[i], ": ", INV$extra[i], "\n\n", sep="")
write.csv(rbind(S, INV), file.path(OUT,"newcohorts_summary.csv"), row.names=FALSE)
cat("WROTE", file.path(OUT,"newcohorts_summary.csv"), nrow(S)+nrow(INV), "rows\n")
