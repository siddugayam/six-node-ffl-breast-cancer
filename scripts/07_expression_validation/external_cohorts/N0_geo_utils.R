## ==========================================================================
## N0_geo_utils.R -- self-contained GEO series-matrix parser + probe collapse
## Written for the new-cohort replication task (v2). No GEOquery dependency
## for parsing; annotation read directly from GEO .annot.gz / GPL soft files.
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})

## ---- parse a series_matrix.txt.gz into list(X = probe x sample, meta = df) --
parse_series_matrix <- function(path){
  ## fast path: pull the metadata lines and the numeric table with two
  ## separate streaming passes instead of readLines() over the whole file.
  meta_raw <- system(paste0("zcat ", shQuote(path), " | grep '^!Sample_'"),
                     intern=TRUE)
  stopifnot(length(meta_raw) > 0)
  tab <- data.table::fread(
    cmd=paste0("zcat ", shQuote(path),
               " | sed -n '/^!series_matrix_table_begin/,/^!series_matrix_table_end/p'",
               " | grep -v '^!series_matrix_table'"),
    sep="\t", header=TRUE, quote='"', na.strings=c("NA","null",""))
  ids <- as.character(tab[[1]])
  vals <- as.matrix(tab[, -1, with=FALSE])
  storage.mode(vals) <- "double"
  rownames(vals) <- ids
  nc <- ncol(vals) + 1L
  mkey <- sub("^!(Sample_[^\t]*)\t.*$", "\\1", meta_raw)
  mrow <- lapply(meta_raw, function(s){
    p <- strsplit(s, "\t", fixed=TRUE)[[1]][-1]; gsub('"','',p)
  })
  ok <- vapply(mrow, length, 1L) == (nc-1L)
  mkey <- mkey[ok]; mrow <- mrow[ok]
  M <- as.data.frame(do.call(rbind, mrow), stringsAsFactors=FALSE)
  colnames(M) <- colnames(vals)
  key2 <- make.unique(mkey, sep="__")
  rownames(M) <- key2
  meta <- as.data.frame(t(M), stringsAsFactors=FALSE)
  meta$geo_accession <- rownames(meta)
  list(X = vals, meta = meta,
       platform = unique(unlist(M[grep("^Sample_platform_id", key2), , drop=TRUE])))
}

## ---- read an NCBI GEO .annot.gz platform annotation -> data.frame(ID, Symbol)
read_gpl_annot <- function(path){
  con <- gzfile(path,"rt"); L <- readLines(con); close(con)
  b <- grep("^!platform_table_begin", L); e <- grep("^!platform_table_end", L)
  tab <- L[(b+1):(e-1)]
  dt <- data.table::fread(text=paste(tab, collapse="\n"), sep="\t",
                          quote="", header=TRUE, colClasses="character")
  nm <- colnames(dt)
  sym_col <- nm[grepl("^Gene symbol$|^GENE_SYMBOL$|^Gene Symbol$|^Symbol$", nm, ignore.case=TRUE)][1]
  stopifnot(!is.na(sym_col))
  data.frame(ID=as.character(dt[[1]]), Symbol=as.character(dt[[sym_col]]),
             stringsAsFactors=FALSE)
}

## ---- read a GPL soft family / annot file downloaded as plain SOFT ----------
read_gpl_soft <- function(path, sym_regex="symbol"){
  con <- if(grepl("\\.gz$",path)) gzfile(path,"rt") else file(path,"rt")
  L <- readLines(con); close(con)
  b <- grep("^!platform_table_begin", L); e <- grep("^!platform_table_end", L)
  tab <- L[(b+1):(e-1)]
  dt <- data.table::fread(text=paste(tab, collapse="\n"), sep="\t",
                          quote="", header=TRUE, colClasses="character")
  nm <- colnames(dt)
  cand <- nm[grepl(sym_regex, nm, ignore.case=TRUE)]
  stopifnot(length(cand)>0)
  data.frame(ID=as.character(dt[[1]]), Symbol=as.character(dt[[cand[1]]]),
             stringsAsFactors=FALSE)
}

## ---- legacy -> current HGNC aliases seen on older GEO platforms -----------
GEO_ALIAS <- c(VEGF="VEGFA", MADH3="SMAD3", MADH4="SMAD4", MADH7="SMAD7",
               MADH2="SMAD2", MADH1="SMAD1", "NFKB3"="RELA", "PPARG1"="PPARG",
               "AML1"="RUNX1", "CBFA2"="RUNX1", "TP53BP2"="TP53BP2",
               "HIF1"="HIF1A", "ARNT2"="ARNT2", "CTNNB"="CTNNB1")

## ---- collapse probes to gene symbols by maximum mean expression ------------
## `priority`: symbols we care about. For a multi-mapping probe ("A///B///C")
## the first token is used unless one of the tokens is in `priority`, in which
## case that token wins -- otherwise genes such as RUNX1 (annotated only as
## "LOC101928269///LOC100506403///RUNX1" on GPL570) are silently lost.
collapse_to_symbol <- function(X, map, priority=character(0)){
  raw <- map$Symbol[match(rownames(X), map$ID)]
  raw <- trimws(ifelse(is.na(raw), "", raw))
  pick <- function(s){
    if(s=="" || s=="---") return(NA_character_)
    tk <- trimws(strsplit(s, "///", fixed=TRUE)[[1]])
    tk <- ifelse(tk %in% names(GEO_ALIAS), GEO_ALIAS[tk], tk)
    hit <- tk[tk %in% priority]
    if(length(hit)) hit[1] else tk[1]
  }
  usym <- unique(raw)
  lut <- vapply(usym, pick, "")
  sym <- unname(lut[match(raw, usym)])
  ok <- !is.na(sym) & sym != "" & sym != "---"
  X <- X[ok, , drop=FALSE]; sym <- sym[ok]
  mu <- rowMeans(X, na.rm=TRUE)
  ord <- order(sym, -mu)
  X <- X[ord, , drop=FALSE]; sym <- sym[ord]
  keep <- !duplicated(sym)
  Y <- X[keep, , drop=FALSE]
  rownames(Y) <- sym[keep]
  stopifnot(!any(duplicated(rownames(Y))))
  Y
}

## ---- helper: pull a "key: value" field out of characteristics columns -----
pull_char <- function(meta, key){
  cols <- grep("^Sample_characteristics_ch1", colnames(meta), value=TRUE)
  out <- rep(NA_character_, nrow(meta))
  pat <- paste0("^\\s*", key, "\\s*:\\s*")
  for(cc in cols){
    v <- meta[[cc]]
    hit <- grepl(pat, v, ignore.case=TRUE)
    out[hit & is.na(out)] <- sub(pat, "", v[hit & is.na(out)], ignore.case=TRUE)
  }
  trimws(out)
}
