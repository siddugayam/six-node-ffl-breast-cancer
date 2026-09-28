# Shared helper: canonicalise a miRBase name (precursor or mature, any case)
# to the canonical mature stem form used by canonical_nodes.tsv:
#   hsa-mir-155        -> hsa-miR-155
#   hsa-mir-29b-1      -> hsa-miR-29b     (strip precursor copy-number suffix)
#   hsa-let-7a-1       -> hsa-let-7a
#   hsa-miR-17-5p      -> hsa-miR-17      (strip arm suffix)
canon_mirna <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("^hsa-mir-", "hsa-miR-", y, ignore.case=TRUE)
  y <- sub("^hsa-let-", "hsa-let-", y, ignore.case=TRUE)
  y <- sub("^HSA-", "hsa-", y)
  # strip arm suffix (-5p / -3p, optionally with .1)
  y <- sub("-[35]p(\\.[0-9]+)?$", "", y)
  # strip precursor copy-number suffix: base = -<digits><optional letters>, then -<digits>
  y <- sub("^(hsa-(miR|let)-[0-9]+[a-zA-Z]*)-[0-9]+$", "\\1", y)
  y
}
