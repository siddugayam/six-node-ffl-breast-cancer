#!/usr/bin/env Rscript
# Exact symbol -> STRING v12 protein-ID mapping using the official STRING
# info/aliases flat files (the /get_string_ids API does FUZZY text matching and
# mis-resolves e.g. VEGFA -> COL18A1, VDR -> CYP27B1; we must not use it).
suppressPackageStartupMessages({library(data.table)})
BASE <- "/path/to/revision"
LOG <- file.path(BASE,"logs","05a_string_map.log")
say <- function(...) { m <- paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse=""))
                       cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE) }
cat("", file=LOG)
nodes <- fread(file.path(BASE,"data","canonical_nodes.tsv"), sep="\t", header=TRUE)
prot  <- nodes[type %in% c("TF","Gene"), name]
say("protein-coding node symbols to map: ", length(prot))

info <- fread(cmd=paste0("zcat ", file.path(BASE,"cache","9606.protein.info.v12.0.txt.gz")),
              sep="\t", header=TRUE, quote="")
setnames(info, 1, "string_protein_id")
ali  <- fread(cmd=paste0("zcat ", file.path(BASE,"cache","9606.protein.aliases.v12.0.txt.gz")),
              sep="\t", header=TRUE, quote="")
setnames(ali, c("string_protein_id","alias","source"))
say("STRING v12 human proteins=", nrow(info), " alias rows=", nrow(ali))

map <- data.table(symbol=prot, ensp=NA_character_, how=NA_character_)
# tier a: exact preferred_name
m <- match(map$symbol, info$preferred_name)
map[!is.na(m), `:=`(ensp=info$string_protein_id[m[!is.na(m)]], how="preferred_name")]
say("mapped by exact preferred_name: ", sum(!is.na(map$ensp)))
# tier b: primary official symbol alias sources
prim <- c("Ensembl_HGNC_symbol","Ensembl_HGNC","BioMart_HUGO","UniProt_GN_Name",
          "Ensembl_EntrezGene","Ensembl_WikiGene")
ab <- unique(ali[source %in% prim, .(alias, string_protein_id)])
ab <- ab[, .N, by=alias][N==1][, alias]
abm <- unique(ali[source %in% prim & alias %in% ab, .(alias, string_protein_id)])
need <- is.na(map$ensp); m2 <- match(map$symbol[need], abm$alias)
map[need][!is.na(m2)]  # no-op, use explicit indices
idx <- which(need)[!is.na(m2)]
map[idx, `:=`(ensp=abm$string_protein_id[m2[!is.na(m2)]], how="official_alias")]
say("mapped after official-symbol aliases: ", sum(!is.na(map$ensp)))
# tier c: HGNC alias symbols / synonyms, only if they resolve uniquely
sec <- c("Ensembl_HGNC_alias_symbol","Ensembl_external_synonym_HGNC","UniProt_GN_Synonyms")
ac <- unique(ali[source %in% sec, .(alias, string_protein_id)])
ac <- ac[, .N, by=alias][N==1][, alias]
acm <- unique(ali[source %in% sec & alias %in% ac, .(alias, string_protein_id)])
need <- is.na(map$ensp); m3 <- match(map$symbol[need], acm$alias)
idx <- which(need)[!is.na(m3)]
map[idx, `:=`(ensp=acm$string_protein_id[m3[!is.na(m3)]], how="hgnc_alias_symbol")]
say("mapped after HGNC alias symbols: ", sum(!is.na(map$ensp)))
say("UNMAPPED symbols (absent from STRING v12 human proteome): ",
    paste(map[is.na(ensp), symbol], collapse=","))
# collisions: two node symbols mapping to the same STRING protein
dup <- map[!is.na(ensp)][, .N, by=ensp][N>1]
if (nrow(dup)) { say("WARNING - symbol collisions on same STRING protein:")
  print(map[ensp %in% dup$ensp][order(ensp)]) }
say("distinct STRING proteins to query: ", uniqueN(map[!is.na(ensp), ensp]))
fwrite(map, file.path(BASE,"cache","string_symbol_map_exact.tsv"), sep="\t", quote=FALSE)
writeLines(unique(map[!is.na(ensp), ensp]), file.path(BASE,"cache","string_ensp_ids.txt"))
say("DONE")
