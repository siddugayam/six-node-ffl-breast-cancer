#!/usr/bin/env Rscript
# 08a: build the node pools for the DISEASE-UNRELATED matched control networks and
#      retrieve their multiMiR *validated* miRNA-target edges (genome-wide per gene).
#
# The control design needs, for genes/TFs that are NOT in the breast-cancer network:
#   - TRRUST TF->target edges          (local file, no fetch needed)
#   - TransmiR TF->miRNA edges         (local file, no fetch needed)
#   - multiMiR validated miRNA->target (MUST be fetched; querying by TARGET returns the
#                                       complete genome-wide miRNA list for that gene, so
#                                       gene in-degrees are directly comparable to BRCA's)
suppressPackageStartupMessages({library(multiMiR); library(data.table)})
BASE  <- "/path/to/revision"
CACHE <- file.path(BASE,"cache","control_multimir"); dir.create(CACHE, showWarnings=FALSE, recursive=TRUE)
LOG   <- file.path(BASE,"logs","motif_control_fetch.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse=""))
                       cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE); flush.console() }
cat("", file=LOG)
set.seed(20260908)

nodes <- fread(file.path(BASE,"data","canonical_nodes.tsv"), sep="\t")
brca_all <- unique(nodes$name)
n_tf   <- sum(nodes$type=="TF"); n_gene <- sum(nodes$type=="Gene"); n_mir <- sum(nodes$type=="miRNA")
say("BRCA network node sizes: TF=",n_tf," Gene=",n_gene," miRNA=",n_mir)

## ---------- universes ----------
sym <- readLines(file.path(BASE,"cache","tcga_gene_symbols.txt"))
sym <- sym[!grepl("^[0-9]+$", sym)]                 # drop unmapped entrez-only ids
prot_universe <- setdiff(sym, brca_all)
say("protein-coding universe (TCGA-BRCA measured symbols, minus BRCA network nodes) = ", length(prot_universe))

tr <- fread(file.path(BASE,"data","db","trrust_human.tsv"), header=FALSE, sep="\t",
            col.names=c("TF","target","mode","pmid"))
tf_universe <- setdiff(unique(tr$TF), brca_all)
say("TRRUST TF universe minus BRCA TFs = ", length(tf_universe), " (TRRUST total ", uniqueN(tr$TF), ")")

## ---------- gene / TF pool to query ----------
POOL_N <- 2000
pool_genes <- sample(prot_universe, POOL_N)
query_targets <- sort(unique(c(pool_genes, tf_universe)))
say("pool_genes=", length(pool_genes), "  tf_universe=", length(tf_universe),
    "  distinct targets to query = ", length(query_targets))
writeLines(pool_genes, file.path(CACHE,"pool_genes.txt"))
writeLines(tf_universe, file.path(CACHE,"pool_tfs.txt"))

## ---------- multiMiR validated retrieval, batched + cached ----------
B <- 25
batches <- split(query_targets, ceiling(seq_along(query_targets)/B))
say("multiMiR validated query: ", length(batches), " batches of <= ", B)
t0 <- Sys.time()
for (i in seq_along(batches)) {
  f <- file.path(CACHE, sprintf("ctrl_batch_%04d.rds", i))
  if (file.exists(f)) next
  ok <- FALSE; att <- 0
  while (!ok && att < 4) { att <- att + 1
    d <- tryCatch({ r <- get_multimir(org="hsa", target=batches[[i]], table="validated",
                                      summary=FALSE, legacy.out=FALSE); as.data.table(r@data) },
                  error=function(e){ say("  batch ",i," attempt ",att," ERROR: ",conditionMessage(e)); NULL })
    if (!is.null(d)) ok <- TRUE else Sys.sleep(15) }
  if (!ok) { say("FAILED batch ", i, " after 4 attempts -- writing empty placeholder"); d <- data.table() }
  saveRDS(d, f)
  if (i %% 10 == 0 || i == length(batches)) {
    el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
    say("  batch ", i, "/", length(batches), "  rows=", nrow(d),
        "  elapsed=", round(el,1), " min  eta=", round(el/i*(length(batches)-i),1), " min")
  }
}
files <- list.files(CACHE, pattern="^ctrl_batch_.*rds$", full.names=TRUE)
V <- rbindlist(lapply(files, readRDS), fill=TRUE)
say("total control validated rows retrieved = ", nrow(V), " from ", length(files), " batch files")
saveRDS(V, file.path(CACHE,"control_validated_all.rds"))
say("DONE 08a")
