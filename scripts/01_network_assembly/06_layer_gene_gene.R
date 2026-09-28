#!/usr/bin/env Rscript
# D) Gene-gene layer: two clearly separated evidence tiers.
#    Tier 1 (TRRUST)  : directed, literature-curated transcriptional regulation onto a network
#                       Gene node, including sources the manuscript never classified as TFs.
#    Tier 2 (STRING)  : UNDIRECTED STRING v12 high-confidence (combined_score >= 900)
#                       interactions; physical-subnetwork and full functional network both run.
#    Identifiers are resolved EXACTLY via the STRING v12 flat files (see 05a) because the
#    /get_string_ids API fuzzy-matches and mis-resolves VEGFA -> COL18A1, VDR -> CYP27B1.
suppressPackageStartupMessages(library(data.table))
BASE <- "/path/to/revision"
LOG  <- file.path(BASE,"logs","06_layer_gene_gene.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse=""))
                       cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE) }
pr <- function(x){ print(x); sink(LOG,append=TRUE); print(x); sink() }
cat("", file=LOG)

nodes <- fread(file.path(BASE,"data","canonical_nodes.tsv"), sep="\t", header=TRUE)
edges <- fread(file.path(BASE,"data","canonical_edges.tsv"), sep="\t", header=TRUE)
prot  <- nodes[type %in% c("TF","Gene"), name]
genes <- nodes[type=="Gene", name]
tfs   <- nodes[type=="TF", name]
say("protein-coding nodes=", length(prot), " Gene=", length(genes), " TF=", length(tfs))

## ---------- Tier 1 : TRRUST ----------
tr <- fread(file.path(BASE,"data","db","trrust_human.tsv"), sep="\t", header=FALSE,
            col.names=c("TF","target","mode","pmids"), quote="")
t1 <- tr[TF %in% prot & target %in% genes]
say("TRRUST rows with source in universe and target a network Gene = ", nrow(t1))
t1a <- t1[, { md <- unique(mode); inf <- setdiff(md,"Unknown")
              .(mode = if(length(inf)==0) "Unknown" else if(length(inf)==1) inf else "Ambiguous")
            }, by=.(source=TF, target)]
say("  -> unique directed TRRUST gene-regulatory pairs = ", nrow(t1a))
t1a[, source_is_manuscript_TF := source %in% tfs]
say("  source is a manuscript-declared TF : ", sum(t1a$source_is_manuscript_TF))
say("  source is a manuscript-declared GENE (strict Gene->Gene) : ", sum(!t1a$source_is_manuscript_TF))
say("  strict Gene->Gene TRRUST edges:"); pr(t1a[source_is_manuscript_TF==FALSE][order(source,target)])
say("  TRRUST tier mode breakdown:"); pr(t1a[, .N, by=mode][order(-N)])
tier1 <- t1a[, .(source, target, evidence="TRRUST", directed=TRUE, score_or_mode=mode,
                 string_network_type=NA_character_, source_is_manuscript_TF)]

## ---------- Tier 2 : STRING v12 (exact ID mapping) ----------
smap <- fread(file.path(BASE,"cache","string_symbol_map_exact.tsv"), sep="\t", header=TRUE)
smap <- smap[!is.na(ensp) & ensp != ""]
say("STRING-mapped node symbols = ", nrow(smap), " / ", length(prot),
    " (absent from STRING v12 human proteome: ",
    paste(setdiff(prot, smap$symbol), collapse=","), ")")
say("  mapping route: "); pr(smap[, .N, by=how])
ensp2sym <- setNames(smap$symbol, smap$ensp)

read_string <- function(f, nettype) {
  d <- fread(f, sep="\t", header=TRUE, quote="")
  d[, a := ensp2sym[stringId_A]]
  d[, b := ensp2sym[stringId_B]]
  bad <- sum(is.na(d$a) | is.na(d$b))
  say("  ", nettype, ": ", nrow(d), " returned rows, ", bad, " dropped (endpoint outside node universe)")
  d <- d[!is.na(a) & !is.na(b) & a != b]
  d[, `:=`(s = pmin(a,b), t = pmax(a,b))]
  d <- d[order(-score)][, .SD[1], by=.(s,t)]
  d[, .(source=s, target=t, score=score, net=nettype)]
}
sp <- read_string(file.path(BASE,"cache","string_physical_900_ensp.tsv"), "physical")
sf <- read_string(file.path(BASE,"cache","string_functional_900_ensp.tsv"), "functional")
say("STRING physical  network, score>=900 : ", nrow(sp), " undirected pairs in node universe")
say("STRING functional network, score>=900 : ", nrow(sf), " undirected pairs in node universe")
say("  physical pairs also in functional set : ",
    sum(paste(sp$source,sp$target) %in% paste(sf$source,sf$target)))
sall <- rbind(sp, sf)[, .(score=max(score), net=paste(sort(unique(net)),collapse="+")), by=.(source,target)]
say("STRING union (physical OR functional, score>=900) : ", nrow(sall), " undirected pairs")
nt <- nodes[, .(name, type)]
sall <- merge(sall, nt[, .(source=name, st=type)], by="source", all.x=TRUE)
sall <- merge(sall, nt[, .(target=name, tt=type)], by="target", all.x=TRUE)
sall[, pairtype := ifelse(st=="Gene" & tt=="Gene", "Gene-Gene",
                   ifelse(st=="TF" & tt=="TF", "TF-TF", "TF-Gene"))]
say("  STRING pairs by node-type combination:"); pr(sall[, .N, by=pairtype][order(-N)])
say("  STRING physical-only subset by node-type:"); 
spx <- merge(sp, nt[,.(source=name, st=type)], by="source"); spx <- merge(spx, nt[,.(target=name, tt=type)], by="target")
spx[, pairtype := ifelse(st=="Gene" & tt=="Gene","Gene-Gene", ifelse(st=="TF" & tt=="TF","TF-TF","TF-Gene"))]
pr(spx[, .N, by=pairtype][order(-N)])
tier2 <- sall[, .(source, target, evidence="STRING", directed=FALSE,
                  score_or_mode=sprintf("%.3f", score), string_network_type=net,
                  source_is_manuscript_TF = source %in% tfs)]

## ---------- emit ----------
out <- rbind(tier1, tier2); setorder(out, evidence, source, target)
fp <- file.path(BASE,"data","layer_gene_gene.tsv")
fwrite(out, fp, sep="\t", quote=FALSE)
say("WROTE ", fp, " rows=", nrow(out), "  (TRRUST tier=", nrow(tier1), ", STRING tier=", nrow(tier2), ")")

## ---------- COL1A1 -> COL3A1 ----------
say("=== the authors' single deposited gene-gene edge: COL1A1 -> COL3A1 ===")
say("TRRUST rows COL1A1->COL3A1 : ", nrow(tr[TF=="COL1A1" & target=="COL3A1"]))
say("TRRUST rows COL3A1->COL1A1 : ", nrow(tr[TF=="COL3A1" & target=="COL1A1"]))
say("TRRUST rows with COL1A1 or COL3A1 as the regulator : ",
    nrow(tr[TF %in% c("COL1A1","COL3A1")]))
pr(out[(source=="COL1A1" & target=="COL3A1") | (source=="COL3A1" & target=="COL1A1")])
say("VERDICT: no directed transcriptional evidence in TRRUST; STRING supports only an",
    " UNDIRECTED high-confidence association.")
say("DONE")
