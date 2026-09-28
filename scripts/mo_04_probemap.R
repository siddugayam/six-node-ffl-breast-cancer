#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
RV<-"/path/to/revision"; CA<-file.path(RV,"cache/multiomics")
ann <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
cat("Illumina 450k annotation:", nrow(ann), "probes\n")
a <- as.data.table(ann[,c("chr","pos","strand","Name","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island")])
fwrite(a, file.path(CA,"anno450k.tsv.gz"), sep="\t")
cat("groups present (first probe of each):\n")
print(head(sort(table(unlist(strsplit(a$UCSC_RefGene_Group,";"))), decreasing=TRUE)))

nodes <- fread(file.path(RV,"data/canonical_nodes.tsv"))
prot  <- nodes[type %in% c("Gene","TF"), name]
bg    <- readLines(file.path(CA,"background_genes.txt"))
PROM  <- c("TSS1500","TSS200","5'UTR","1stExon")

## explode gene/group pairs
ex <- a[UCSC_RefGene_Name!="", .(Name, chr, pos,
        g=strsplit(UCSC_RefGene_Name,";"), gr=strsplit(UCSC_RefGene_Group,";"))]
ex <- ex[, .(gene=unlist(g), grp=unlist(gr)), by=.(Name,chr,pos)]
cat("gene-probe annotations:", nrow(ex), "\n")
pm <- unique(ex[grp %in% PROM, .(Name, chr, pos, gene, grp)])
cat("promoter-proximal gene-probe pairs (TSS1500/TSS200/5'UTR/1stExon):", nrow(pm),
    " unique probes:", uniqueN(pm$Name), "\n")

g_net <- pm[gene %in% prot]; g_bg <- pm[gene %in% bg]
cat("network protein-coding genes with >=1 promoter probe:", uniqueN(g_net$gene), "/", length(prot), "\n")
cat("background genes with >=1 promoter probe:", uniqueN(g_bg$gene), "/", length(bg), "\n")
fwrite(g_net, file.path(CA,"promoter_probes_network_genes.tsv"), sep="\t")
fwrite(g_bg,  file.path(CA,"promoter_probes_background_genes.tsv"), sep="\t")

## ---- miRNA promoter probes: window upstream of precursor, strand aware ----
ml <- fread(file.path(CA,"mirna_loci_hg19.tsv"))
cat("miRNA loci:", nrow(ml), "for", uniqueN(ml$mature), "mature miRNAs\n")
UP <- 2000; DN <- 500
ml[, wstart := ifelse(strand=="+", start-UP, end-DN)]
ml[, wend   := ifelse(strand=="+", start+DN, end+UP)]
setkey(a, chr, pos)
mp <- rbindlist(lapply(seq_len(nrow(ml)), function(i){
  r <- ml[i]; s <- a[chr==r$chrom & pos>=r$wstart & pos<=r$wend]
  if(!nrow(s)) return(NULL)
  data.table(mature=r$mature, precursor=r$precursor, sym=r$sym, Name=s$Name, chr=s$chr, pos=s$pos,
             src="window")
}))
cat("probes in +/-", UP, "/", DN, "bp precursor window:", nrow(mp), "\n")
## plus probes annotated by Illumina to the MIR gene symbol (any group, since MIR genes are short)
ms <- ex[gene %in% ml$sym]
cat("probes annotated to a MIR* symbol:", nrow(ms), "\n")
ms2 <- merge(ms[,.(sym=gene,Name,chr,pos)], unique(ml[,.(mature,precursor,sym)]), by="sym", allow.cartesian=TRUE)
ms2[, src := "mir_symbol"]
mall <- unique(rbind(mp[,.(mature,precursor,sym,Name,chr,pos,src)],
                     ms2[,.(mature,precursor,sym,Name,chr,pos,src)]), by=c("mature","Name"))
cat("union miRNA promoter probes:", nrow(mall), " unique probes:", uniqueN(mall$Name),
    " miRNAs covered:", uniqueN(mall$mature), "/", uniqueN(ml$mature), "\n")
fwrite(mall, file.path(CA,"promoter_probes_mirnas.tsv"), sep="\t")

feat <- c("hsa-miR-130a","hsa-miR-124","hsa-miR-101","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c",
          "hsa-let-7b","hsa-let-7e","hsa-miR-34a","hsa-miR-200b","hsa-miR-200c","hsa-miR-145")
cat("\nprobe counts for featured miRNAs:\n")
print(mall[mature %in% feat, .N, by=.(mature)][order(mature)])
cat("featured miRNAs with zero probes:", paste(setdiff(feat, mall$mature), collapse=", "), "\n")

allprobes <- unique(c(g_net$Name, g_bg$Name, mall$Name))
writeLines(allprobes, file.path(CA,"probes_final.txt"))
cat("\nTOTAL probes to extract:", length(allprobes), "\n")
