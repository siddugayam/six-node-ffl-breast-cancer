#!/usr/bin/env Rscript
# C) miRNA-miRNA layer defined as POLYCISTRONIC CO-TRANSCRIPTION (genomic clustering),
#    NOT as invented direct regulatory interactions.
#    Definition: two miRNA precursors on the same chromosome and strand whose primary
#    transcripts are separated by <= T bp belong to one polycistronic cluster
#    (single-linkage chaining). Altuvia 2005; Baskerville & Bartel 2005; Saini 2007.
suppressPackageStartupMessages(library(data.table))
BASE <- "/path/to/revision"
source(file.path(BASE,"scripts","01_network_assembly","01_mirna_canon.R"))
LOG <- file.path(BASE,"logs","05_layer_miRNA_miRNA.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse=""))
                       cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE) }
pr <- function(x){ print(x); sink(LOG,append=TRUE); print(x); sink() }
cat("", file=LOG)

nodes <- fread(file.path(BASE,"data","canonical_nodes.tsv"), sep="\t", header=TRUE)
edges <- fread(file.path(BASE,"data","canonical_edges.tsv"), sep="\t", header=TRUE)
mirs  <- nodes[type=="miRNA", name]
say("network miRNA nodes=", length(mirs))

## ---- parse miRBase GFF3 primary transcripts ----
raw <- fread(file.path(BASE,"data","db","hsa.gff3"), sep="\t", header=FALSE,
             skip="chr", quote="", col.names=c("chrom","src","feature","start","end","score","strand","frame","attr"))
say("gff3 feature rows=", nrow(raw))
loci <- raw[feature=="miRNA_primary_transcript"]
say("miRNA_primary_transcript loci=", nrow(loci))
loci[, precursor := sub(".*Name=([^;]+).*", "\\1", attr)]
loci[, mirbase_id := sub(".*ID=([^;]+).*", "\\1", attr)]
loci[, canon := canon_mirna(precursor)]
loci <- loci[, .(chrom, start, end, strand, precursor, mirbase_id, canon)]
setorder(loci, chrom, strand, start)
say("distinct canonical names from loci=", uniqueN(loci$canon))

mapped <- intersect(mirs, unique(loci$canon))
say("network miRNAs mapped to >=1 miRBase locus = ", length(mapped),
    " / ", length(mirs), " (unmapped: ", paste(setdiff(mirs, mapped), collapse=","), ")")

## ---- single-linkage chaining clustering over ALL miRBase loci ----
cluster_at <- function(dt, thr) {
  d <- copy(dt); setorder(d, chrom, strand, start)
  d[, gap := start - shift(end), by=.(chrom, strand)]
  d[, newcl := is.na(gap) | gap > thr]
  d[, cid := cumsum(newcl)]
  d[, .(chrom, strand, start, end, precursor, mirbase_id, canon,
        cluster_id = paste0("cl", thr/1000, "kb_", cid))]
}
THR <- c("3kb"=3000, "10kb"=10000, "50kb"=50000)
cl <- lapply(THR, function(t) cluster_at(loci, t))
for (nm in names(cl)) {
  d <- cl[[nm]]
  sz <- d[, .N, by=cluster_id]
  say("threshold ", nm, ": ", nrow(sz), " clusters over ", nrow(d), " loci; ",
      nrow(sz[N>=2]), " multi-locus clusters; largest=", max(sz$N))
}

## ---- restrict to node universe, build undirected pairs ----
pairs_at <- function(d) {
  dn <- d[canon %in% mirs]
  out <- dn[, {
    if (.N < 2) NULL else {
      cmb <- as.data.table(t(combn(.N, 2)))
      setnames(cmb, c("i","j"))
      cmb[, .(a = canon[i], b = canon[j],
              chrom = chrom[i], strand = strand[i],
              dist = pmax(0L, as.integer(pmax(start[i], start[j]) - pmin(end[i], end[j]))),
              cluster_id = cluster_id[1])]
    }
  }, by=cluster_id]
  out <- out[a != b]
  if (nrow(out)==0) return(out)
  out[, `:=`(m1 = pmin(a,b), m2 = pmax(a,b))]
  out <- out[order(dist)][, .SD[1], by=.(m1,m2)]
  out[, .(miRNA_1=m1, miRNA_2=m2, chrom, strand, distance_bp=dist, cluster_id)]
}
p3  <- pairs_at(cl[["3kb"]]);  p10 <- pairs_at(cl[["10kb"]]); p50 <- pairs_at(cl[["50kb"]])
say("network-miRNA pairs co-clustered:  3kb=", nrow(p3), "  10kb=", nrow(p10), "  50kb=", nrow(p50))

## union of pairs (50kb is the loosest and is a superset)
u <- unique(rbind(p3[,.(miRNA_1,miRNA_2)], p10[,.(miRNA_1,miRNA_2)], p50[,.(miRNA_1,miRNA_2)]))
key3  <- paste(p3$miRNA_1,  p3$miRNA_2)
key10 <- paste(p10$miRNA_1, p10$miRNA_2)
key50 <- paste(p50$miRNA_1, p50$miRNA_2)
u[, k := paste(miRNA_1, miRNA_2)]
u <- merge(u, p10[, .(k=paste(miRNA_1,miRNA_2), chrom, distance_bp, cluster_id)], by="k", all.x=TRUE)
u <- merge(u, p50[, .(k=paste(miRNA_1,miRNA_2), chrom50=chrom, dist50=distance_bp, cid50=cluster_id)],
           by="k", all.x=TRUE)
u[is.na(chrom), chrom := chrom50]
u[is.na(distance_bp), distance_bp := dist50]
u[, `:=`(threshold_10kb = k %in% key10, threshold_3kb = k %in% key3, threshold_50kb = k %in% key50)]
u[threshold_10kb==FALSE, cluster_id := NA_character_]   # cluster_id refers to the 10kb definition
out <- u[, .(miRNA_1, miRNA_2, chrom, distance_bp, cluster_id, threshold_10kb, threshold_3kb, threshold_50kb)]
setorder(out, -threshold_3kb, -threshold_10kb, chrom, distance_bp)
fp <- file.path(BASE,"data","layer_miRNA_miRNA.tsv")
fwrite(out, fp, sep="\t", quote=FALSE)
say("WROTE ", fp, " rows=", nrow(out))
say("  rows with threshold_3kb TRUE=", sum(out$threshold_3kb),
    "  10kb TRUE=", sum(out$threshold_10kb), "  50kb TRUE=", sum(out$threshold_50kb))

## ---- network miRNAs in multi-miRNA clusters ----
for (nm in names(cl)) {
  d <- cl[[nm]][canon %in% mirs]
  # clusters containing >=2 DISTINCT network miRNAs
  sz <- d[, .(n_distinct = uniqueN(canon)), by=cluster_id]
  keep <- sz[n_distinct>=2, cluster_id]
  say("threshold ", nm, ": network miRNAs in a cluster with >=2 distinct network miRNAs = ",
      uniqueN(d[cluster_id %in% keep, canon]), " / ", length(mirs))
}

## largest 10kb clusters (network miRNAs only)
d10 <- cl[["10kb"]][canon %in% mirs]
top <- d10[, .(n=uniqueN(canon), members=paste(sort(unique(canon)), collapse=",")), by=.(cluster_id)][n>=2][order(-n)]
say("--- largest 10kb polycistronic clusters (network miRNAs) : ", nrow(top), " multi-miRNA clusters ---")
pr(head(top, 15))
fwrite(top, file.path(BASE,"results","miRNA_clusters_10kb.tsv"), sep="\t", quote=FALSE)

## also full cluster tables for all 3 thresholds
for (nm in names(cl)) {
  dd <- cl[[nm]][canon %in% mirs]
  tt <- dd[, .(n=uniqueN(canon), members=paste(sort(unique(canon)),collapse=",")), by=cluster_id][n>=2][order(-n)]
  fwrite(tt, file.path(BASE,"results",paste0("miRNA_clusters_",nm,".tsv")), sep="\t", quote=FALSE)
}

## ---- audit the 30 author exemplar miRNA_miRNA edges ----
au <- edges[edge_type=="miRNA_miRNA", .(source, target)]
au[, `:=`(m1=pmin(source,target), m2=pmax(source,target))]
au <- unique(au[, .(m1,m2)])
say("author exemplar miRNA_miRNA edges=", nrow(edges[edge_type=="miRNA_miRNA"]), " unique undirected pairs=", nrow(au))
au[, k := paste(m1,m2)]
au[, `:=`(coclust_3kb = k %in% key3, coclust_10kb = k %in% key10, coclust_50kb = k %in% key50)]
say("author pairs genomically co-clustered:  3kb=", sum(au$coclust_3kb),
    "  10kb=", sum(au$coclust_10kb), "  50kb=", sum(au$coclust_50kb))
say("--- author pairs, co-clustering status ---")
pr(au[, .(m1,m2,coclust_3kb,coclust_10kb,coclust_50kb)])
fwrite(au[, .(miRNA_1=m1, miRNA_2=m2, coclust_3kb, coclust_10kb, coclust_50kb)],
       file.path(BASE,"results","author_miRNA_miRNA_audit.tsv"), sep="\t", quote=FALSE)
say("DONE")
