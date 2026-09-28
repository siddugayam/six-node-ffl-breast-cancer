#!/usr/bin/env Rscript
# B) TF-miRNA layer with literature-derived signs, from TransmiR v2 (human).
suppressPackageStartupMessages(library(data.table))
BASE <- "/path/to/revision"
source(file.path(BASE,"scripts","01_network_assembly","01_mirna_canon.R"))
LOG <- file.path(BASE,"logs","04_layer_TF_miRNA.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse=""))
                       cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE) }
cat("", file=LOG)

nodes <- fread(file.path(BASE,"data","canonical_nodes.tsv"), sep="\t", header=TRUE)
edges <- fread(file.path(BASE,"data","canonical_edges.tsv"), sep="\t", header=TRUE)
prot  <- nodes[type %in% c("TF","Gene"), name]
mirs  <- nodes[type=="miRNA", name]
say("protein nodes=",length(prot)," miRNA nodes=",length(mirs))

tm <- fread(file.path(BASE,"data","db","transmir_hsa.tsv"), sep="\t", header=FALSE, quote="")
setnames(tm, paste0("V",1:ncol(tm)))
say("TransmiR rows read=", nrow(tm))
tm[, TF := trimws(V1)]
tm[, mirna_raw := trimws(V2)]
tm[, mirna := canon_mirna(mirna_raw)]
tm[, mode_raw := trimws(V5)]
tm[, pmid := trimws(V6)]
tm[, evidence := trimws(V7)]

# normalise regulation type -> mode / sign
norm_mode <- function(m) {
  ml <- tolower(m)
  fifelse(grepl("^activation", ml) | grepl("^activation\\(", ml), "Activation",
   fifelse(grepl("^repression", ml), "Repression", "Regulation"))
}
tm[, mode := norm_mode(mode_raw)]
# feedback-flagged entries keep sign but are marked
tm[, feedback := grepl("feedback|loop|circuit|network", tolower(mode_raw))]
# a few odd strings: "Activation(a negative feedback loop)" -> Activation ; treat
# "Regulation(Double-Negative Feedback Loop)" / "auto-regulatory..." as Regulation (sign 0)
tm[, sign := fifelse(mode=="Activation", 1L, fifelse(mode=="Repression", -1L, 0L))]

say("--- TransmiR mode normalisation over ALL rows ---")
print(tm[, .N, by=.(mode_raw, mode)][order(-N)])
sink(LOG,append=TRUE); print(tm[, .N, by=.(mode_raw,mode)][order(-N)]); sink()

sub <- tm[TF %in% prot & mirna %in% mirs]
say("TransmiR rows with TF and miRNA both in node universe=", nrow(sub))
say("  distinct TF-miRNA pairs=", nrow(unique(sub[, .(TF, mirna)])))

agg <- sub[, {
  md <- unique(mode); inf <- setdiff(md, "Regulation")
  m <- if (length(inf)==0) "Regulation" else if (length(inf)==1) inf else "Ambiguous"
  p <- unique(unlist(strsplit(paste(pmid, collapse=";"), "[;,]"))); p <- p[nzchar(p) & p!="n/a"]
  e <- paste(sort(unique(evidence)), collapse=";")
  .(mode=m, pmid=paste(sort(p), collapse=";"), evidence=e)
}, by=.(source=TF, target=mirna)]
agg[, sign := fifelse(mode=="Activation",1L, fifelse(mode=="Repression",-1L,0L))]
setcolorder(agg, c("source","target","mode","sign","pmid","evidence"))
setorder(agg, source, target)
out <- file.path(BASE,"data","layer_TF_miRNA.tsv")
fwrite(agg, out, sep="\t", quote=FALSE)
say("WROTE ",out," rows=", nrow(agg))

say("--- mode breakdown of the TF_miRNA layer (node universe) ---")
print(agg[, .N, by=mode][order(-N)]); sink(LOG,append=TRUE); print(agg[,.N,by=mode][order(-N)]); sink()

# coverage of published TF_miRNA edges
pub <- unique(edges[edge_type=="TF_miRNA", .(source, target)])
say("published TF_miRNA edges=", nrow(edges[edge_type=="TF_miRNA"]), " unique pairs=", nrow(pub))
cov <- merge(pub, agg[, .(source,target,mode,sign)], by=c("source","target"), all.x=TRUE)
cov[, covered := !is.na(mode)]
say("published TF_miRNA edges COVERED by TransmiR=", sum(cov$covered),
    " (", round(100*mean(cov$covered),1), "%)")
tab <- cov[covered==TRUE, .N, by=mode][order(-N)]
print(tab); sink(LOG,append=TRUE); print(tab); sink()
say("published TF_miRNA edges annotated REPRESSION=", nrow(cov[covered==TRUE & mode=="Repression"]))
say("published TF_miRNA edges annotated ACTIVATION=", nrow(cov[covered==TRUE & mode=="Activation"]))
say("published TF_miRNA edges annotated REGULATION(unsigned)=", nrow(cov[covered==TRUE & mode=="Regulation"]))
say("published TF_miRNA edges annotated AMBIGUOUS=", nrow(cov[covered==TRUE & mode=="Ambiguous"]))
say("published TF_miRNA edges NOT in TransmiR=", sum(!cov$covered))
fwrite(cov, file.path(BASE,"results","TF_miRNA_TransmiR_coverage.tsv"), sep="\t", quote=FALSE)
say("WROTE results/TF_miRNA_TransmiR_coverage.tsv rows=", nrow(cov))
say("DONE")
