#!/usr/bin/env Rscript
# A) TF-TF and TF-gene layer with literature-derived signs, from TRRUST v2 (human).
# Restricted to the 587-node canonical universe.
suppressPackageStartupMessages(library(data.table))
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs", "02_layer_TF_target.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S"), " | ", paste0(..., collapse=""))
                       cat(m, "\n"); cat(m, "\n", file=LOG, append=TRUE) }
cat("", file=LOG)

nodes <- fread(file.path(BASE,"data","canonical_nodes.tsv"), sep="\t", header=TRUE)
edges <- fread(file.path(BASE,"data","canonical_edges.tsv"), sep="\t", header=TRUE)
say("nodes=", nrow(nodes), " edges=", nrow(edges))

universe   <- nodes$name
prot_nodes <- nodes[type %in% c("TF","Gene"), name]     # protein-coding universe
tf_nodes   <- nodes[type == "TF", name]
say("protein-coding nodes=", length(prot_nodes), " TF nodes=", length(tf_nodes))

tr <- fread(file.path(BASE,"data","db","trrust_human.tsv"), sep="\t", header=FALSE,
            col.names=c("TF","target","mode","pmids"), quote="")
say("TRRUST rows read=", nrow(tr))

# --- keep rows where BOTH TF and target are in the node universe ---
sub <- tr[TF %in% prot_nodes & target %in% prot_nodes]
say("TRRUST rows with both partners in node universe=", nrow(sub))

# --- aggregate duplicate (TF,target) rows ---
sub[, pmids := gsub(" ", "", pmids)]
agg <- sub[, {
  md <- unique(mode)
  inf <- setdiff(md, "Unknown")
  m <- if (length(inf) == 0) "Unknown"
       else if (length(inf) == 1) inf
       else "Ambiguous"
  p <- unique(unlist(strsplit(paste(pmids, collapse=";"), "[;,]")))
  p <- p[nzchar(p)]
  .(mode = m, pmids = paste(sort(p), collapse=";"), n_pmid = length(p))
}, by=.(source=TF, target)]
agg[, sign := fifelse(mode=="Activation", 1L, fifelse(mode=="Repression", -1L, 0L))]
setcolorder(agg, c("source","target","mode","sign","pmids","n_pmid"))
setorder(agg, source, target)

out <- file.path(BASE,"data","layer_TF_target.tsv")
fwrite(agg, out, sep="\t", quote=FALSE)
say("WROTE ", out, " rows=", nrow(agg))

say("--- mode breakdown (all TRRUST edges inside universe) ---")
print(agg[, .N, by=mode][order(-N)])
sink(LOG, append=TRUE); print(agg[, .N, by=mode][order(-N)]); sink()

# how many involve a TF->TF edge vs TF->Gene
agg2 <- merge(agg, nodes[, .(source=name, src_type=type)], by="source", all.x=TRUE)
agg2 <- merge(agg2, nodes[, .(target=name, tgt_type=type)], by="target", all.x=TRUE)
say("--- layer split by node types ---")
print(agg2[, .N, by=.(src_type, tgt_type)][order(-N)])
sink(LOG, append=TRUE); print(agg2[, .N, by=.(src_type,tgt_type)][order(-N)]); sink()
# TF-TF sub-layer specifically
say("TF->TF edges in TRRUST layer = ", nrow(agg2[src_type=="TF" & tgt_type=="TF"]))
say("TF->Gene edges in TRRUST layer = ", nrow(agg2[src_type=="TF" & tgt_type=="Gene"]))
say("non-TF source edges (TRRUST source not labelled TF by authors) = ",
    nrow(agg2[src_type!="TF"]))

# --- CRUCIAL: coverage of the 770 published TF_target edges ---
pub <- edges[edge_type=="TF_target", .(source, target)]
say("published TF_target edges = ", nrow(pub), " (unique pairs=", nrow(unique(pub)), ")")
pub <- unique(pub)
cov <- merge(pub, agg[, .(source, target, mode, sign)], by=c("source","target"), all.x=TRUE)
cov[, covered := !is.na(mode)]
say("published TF_target edges COVERED by TRRUST = ", sum(cov$covered),
    " (", round(100*mean(cov$covered),1), "%)")
say("--- of the covered ones, mode breakdown ---")
tab <- cov[covered==TRUE, .N, by=mode][order(-N)]
print(tab); sink(LOG, append=TRUE); print(tab); sink()
say("published TF_target edges annotated REPRESSION by TRRUST = ",
    nrow(cov[covered==TRUE & mode=="Repression"]))
say("published TF_target edges annotated ACTIVATION by TRRUST = ",
    nrow(cov[covered==TRUE & mode=="Activation"]))
say("published TF_target edges annotated UNKNOWN by TRRUST = ",
    nrow(cov[covered==TRUE & mode=="Unknown"]))
say("published TF_target edges annotated AMBIGUOUS by TRRUST = ",
    nrow(cov[covered==TRUE & mode=="Ambiguous"]))
say("=> manuscript assigned sign=+1 to ALL ", nrow(pub),
    "; sign is demonstrably WRONG for ", nrow(cov[covered==TRUE & mode=="Repression"]),
    " and unsupported for ", nrow(cov[covered==TRUE & mode %in% c("Unknown","Ambiguous")]),
    " covered edges; ", sum(!cov$covered), " have no TRRUST record at all.")

fwrite(cov, file.path(BASE,"results","TF_target_TRRUST_coverage.tsv"), sep="\t", quote=FALSE)
say("WROTE results/TF_target_TRRUST_coverage.tsv rows=", nrow(cov))
say("DONE")
