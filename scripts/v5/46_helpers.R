## shared GEO series-matrix parsing helpers (extracted from 46_neoadjuvant_response.R)
GEO <- file.path("/path/to/revision", "cache/v5/geo")
read_series <- function(gse) {
  f <- file.path(GEO, paste0(gse, ".txt.gz"))
  L <- readLines(f, warn = FALSE)
  acc <- strsplit(grep("^!Sample_geo_accession", L, value = TRUE)[1], "\t")[[1]][-1]
  acc <- gsub('"', "", acc)
  ch <- grep("^!Sample_characteristics_ch1", L, value = TRUE)
  meta <- vector("list", length(acc)); names(meta) <- acc
  for (i in seq_along(acc)) meta[[i]] <- list()
  for (line in ch) {
    v <- gsub('"', "", strsplit(line, "\t")[[1]][-1])
    for (i in seq_along(v)) {
      if (!grepl(":", v[i], fixed = TRUE)) next
      kk <- trimws(sub(":.*$", "", v[i])); vv <- trimws(sub("^[^:]*:", "", v[i]))
      meta[[i]][[kk]] <- vv
    }
  }
  st <- grep("^!series_matrix_table_begin", L); en <- grep("^!series_matrix_table_end", L)
  tab <- fread(text = L[(st + 1):(en - 1)], sep = "\t", header = TRUE, data.table = FALSE)
  rn <- gsub('"', "", tab[[1]]); tab <- tab[, -1, drop = FALSE]
  colnames(tab) <- gsub('"', "", colnames(tab))
  X <- as.matrix(tab); rownames(X) <- rn
  X <- X[, acc[acc %in% colnames(X)], drop = FALSE]
  list(X = X, meta = meta[colnames(X)])
}
getf <- function(meta, key) sapply(meta, function(m) if (is.null(m[[key]])) NA_character_ else m[[key]])

## ------------------------------------------------------------------ probe -> symbol ----
annot_map <- function(platform) {
  f <- file.path(GEO, paste0(platform, ".annot.gz"))
  L <- readLines(f, warn = FALSE)
  st <- grep("^!platform_table_begin", L); en <- grep("^!platform_table_end", L)
  a <- fread(text = L[(st + 1):(en - 1)], sep = "\t", header = TRUE, data.table = TRUE,
             quote = "", fill = TRUE)
  a <- a[, .(ID, sym = `Gene symbol`)]
  a <- a[sym != "" & !is.na(sym)]
  a[, sym := sub("///.*$", "", sym)]        # first symbol of multi-mapping probes
  a
}
collapse_to_symbol <- function(X, amap) {
  m <- amap[match(rownames(X), ID)]
  keep <- !is.na(m$sym)
  X <- X[keep, , drop = FALSE]; sym <- m$sym[keep]
  mu <- rowMeans(X, na.rm = TRUE)
  ord <- order(sym, -mu)
  X <- X[ord, , drop = FALSE]; sym <- sym[ord]
  first <- !duplicated(sym)
  Y <- X[first, , drop = FALSE]; rownames(Y) <- sym[first]
  Y
}

