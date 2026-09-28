#!/usr/bin/env Rscript
# 10a: export, as flat TSVs, the three edge sources used to rebuild BOTH the breast-cancer
#      network and the disease-unrelated matched control networks under an IDENTICAL pipeline:
#         TRRUST   TF -> target        (local)
#         TransmiR TF -> miRNA         (local)
#         multiMiR validated miRNA -> target  (BRCA cache + control cache)
#      Querying multiMiR BY TARGET returns the complete genome-wide miRNA list for that gene,
#      so gene in-degrees are directly comparable between the BRCA and the control gene sets.
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
source(file.path(BASE,"scripts","00_mirna_canon.R"))
LOG <- file.path(BASE,"logs","motif_specificity.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse=""))
                       cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE); flush.console() }
cat("", file=LOG)
D <- file.path(BASE,"data","control"); dir.create(D, showWarnings=FALSE, recursive=TRUE)

## ---- TRRUST ----
tr <- fread(file.path(BASE,"data","db","trrust_human.tsv"), header=FALSE, sep="\t",
            col.names=c("TF","target","mode","pmid"))
tr <- unique(tr[TF != target, .(TF, target)])
fwrite(tr, file.path(D,"trrust_pairs.tsv"), sep="\t")
say("TRRUST unique directed TF->target pairs = ", nrow(tr),
    "  (TFs ", uniqueN(tr$TF), ", targets ", uniqueN(tr$target), ")")

## ---- TransmiR ----
tm <- fread(file.path(BASE,"data","db","transmir_hsa.tsv"), header=FALSE, sep="\t")
setnames(tm, 1:2, c("TF","mirna_raw"))
tm[, mirna := canon_mirna(mirna_raw)]
tm <- unique(tm[, .(TF, mirna)])
fwrite(tm, file.path(D,"transmir_pairs.tsv"), sep="\t")
say("TransmiR unique TF->miRNA pairs (canonical miRNA names) = ", nrow(tm),
    "  (TFs ", uniqueN(tm$TF), ", miRNAs ", uniqueN(tm$mirna), ")")

## ---- multiMiR validated : BRCA cache + control cache ----
grab <- function(files, tag) {
  if (!length(files)) return(data.table())
  V <- rbindlist(lapply(files, readRDS), fill=TRUE)
  if (!nrow(V)) return(data.table())
  V[, st_l := tolower(trimws(support_type))]
  V <- V[!grepl("non-functional", st_l)]          # only miRTarBase Non-Functional MTI is refuted
  V[, mirna := canon_mirna(mature_mirna_id)]
  V[, gene  := toupper(target_symbol)]
  out <- unique(V[nzchar(mirna) & nzchar(gene), .(mirna, gene)])
  out[, pool := tag]
  say(tag, ": ", length(files), " cache files -> ", nrow(V), " supporting records -> ",
      nrow(out), " unique miRNA-gene pairs (", uniqueN(out$mirna), " miRNAs, ",
      uniqueN(out$gene), " genes)")
  out
}
brca_f <- list.files(file.path(BASE,"cache","multimir"), pattern="^validated_batch_.*rds$", full.names=TRUE)
ctrl_f <- list.files(file.path(BASE,"cache","control_multimir"), pattern="^ctrl_batch_.*rds$", full.names=TRUE)
V <- unique(rbindlist(list(grab(brca_f,"BRCA"), grab(ctrl_f,"CONTROL")), fill=TRUE)[, .(mirna, gene)])
fwrite(V, file.path(D,"validated_mirna_target_pairs.tsv"), sep="\t")
say("UNION validated miRNA->target pairs = ", nrow(V), " over ", uniqueN(V$gene),
    " queried genes and ", uniqueN(V$mirna), " miRNAs")

## ---- the exact gene universe that was queried (needed so that a gene with zero validated
##      miRNAs is distinguished from a gene that was never queried) ----
qb <- unique(fread(file.path(BASE,"data","canonical_edges.tsv"))[edge_type=="miRNA_target", target])
qc <- c(readLines(file.path(BASE,"cache","control_multimir","pool_genes.txt")),
        readLines(file.path(BASE,"cache","control_multimir","pool_tfs.txt")))
fwrite(data.table(gene=sort(unique(c(qb,qc)))), file.path(D,"queried_gene_universe.tsv"), sep="\t")
say("queried gene universe (BRCA targets ", length(qb), " + control pool ", length(unique(qc)),
    ") = ", length(unique(c(qb,qc))))

## ---- miRBase mature-stem universe (for uniform random miRNA sampling) ----
g <- readLines(file.path(BASE,"data","db","hsa.gff3"))
g <- g[!startsWith(g,"#")]
f <- tstrsplit(g, "\t", fixed=TRUE)
nm <- sub(".*Name=([^;]+).*", "\\1", f[[9]])
mir_all <- sort(unique(canon_mirna(nm[f[[3]]=="miRNA"])))
writeLines(mir_all, file.path(D,"mirbase_mature_stems.txt"))
say("miRBase mature stems (canonicalised, unique) = ", length(mir_all))
say("DONE 10a")
