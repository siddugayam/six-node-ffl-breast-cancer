#!/usr/bin/env Rscript
## ============================================================================
## 08_exir_prep_run.R
## ExIR (Experimental-data-based Integrative Ranking) on TCGA-BRCA
## Package: influential 2.3.2 (Salavaty, Ramialison & Currie)
## classification.
##
## NOTE ON THE NETWORK: exir() does NOT accept an external graph. It
## reconstructs its own feature-feature association network internally from
## Exptl_data with fcor() (Spearman + mutual-rank). The canonical FFL network
## is therefore NOT an input here -- which is exactly what makes the
## driver-vs-FFL-hub overlap test in step 09 an INDEPENDENT validation.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(influential)
  library(igraph)
})

ROOT <- "/path/to/revision"
LOG  <- file.path(ROOT, "logs", "exir_agent.log")
dir.create(dirname(LOG), showWarnings = FALSE, recursive = TRUE)

say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = LOG, append = TRUE)
  cat(msg, "\n")
  flush.console()
}

say("=== 08_exir_prep_run.R START ===")
say("R: ", R.version.string, " | influential ", as.character(packageVersion("influential")))

## ---------------------------------------------------------------- 1. INPUTS
gexp <- readRDS(file.path(ROOT, "data", "brca_gene_expr.rds"))
mexp <- readRDS(file.path(ROOT, "data", "brca_mirna_expr_canonical.rds"))
phen <- readRDS(file.path(ROOT, "data", "brca_pheno.rds"))
de   <- fread(file.path(ROOT, "results", "BRCA_DEX_ALL_nodes.csv"))

say("gene matrix  : ", nrow(gexp), " x ", ncol(gexp))
say("miRNA matrix : ", nrow(mexp), " x ", ncol(mexp))
say("DE table rows: ", nrow(de), " (cols: ", paste(names(de), collapse = ","), ")")

## sanity: DE table must be the full 20745-feature version, not a truncated one
stopifnot(nrow(de) == 20745)

## ------------------------------------------------- 2. COMMON SAMPLE UNIVERSE
shared <- intersect(colnames(gexp), colnames(mexp))
say("samples shared between gene and miRNA assays: ", length(shared))

phen <- as.data.table(phen)
setkey(phen, sample)
cond_raw <- phen[shared, sample_type]
stopifnot(!any(is.na(cond_raw)))
condition <- ifelse(cond_raw == "Solid Tissue Normal", "Normal",
             ifelse(cond_raw == "Primary Tumor", "Tumor", NA))
keep <- !is.na(condition)
shared    <- shared[keep]
condition <- condition[keep]
say("after dropping non Tumor/Normal samples: n = ", length(shared),
    "  Tumor = ", sum(condition == "Tumor"), "  Normal = ", sum(condition == "Normal"))

## ------------------------------------------------- 3. COMBINED EXPRESSION
g <- gexp[, shared, drop = FALSE]
m <- mexp[, shared, drop = FALSE]
rm(gexp, mexp); gc()

ov <- intersect(rownames(g), rownames(m))
say("feature-name collisions gene vs miRNA: ", length(ov),
    if (length(ov)) paste0(" [", paste(ov, collapse = ","), "]") else "")
stopifnot(length(ov) == 0)

expr <- rbind(g, m)
rm(g, m); gc()
say("combined Exptl_data: ", nrow(expr), " features x ", ncol(expr), " samples")
say("value range: ", paste(round(range(expr), 3), collapse = " .. "),
    "  NAs: ", sum(is.na(expr)))

## ------------------------------------------------- 4. Diff_data
## significance rule identical to script 01: |logFC| > 1 AND adj.P.Val < 0.05
sigdt <- de[abs(logFC) > 1 & adj.P.Val < 0.05]
say("significant DE features (|logFC|>1 & adj.P<0.05): ", nrow(sigdt),
    "  (genes ", sum(sigdt$class == "gene"), ", miRNAs ", sum(sigdt$class == "miRNA"), ")")

sigdt <- sigdt[Gene %in% rownames(expr)]
say("of which present in the combined expression matrix: ", nrow(sigdt))
stopifnot(!any(duplicated(sigdt$Gene)))

Diff_data <- data.frame(logFC     = sigdt$logFC,
                        adj.P.Val = sigdt$adj.P.Val,
                        row.names = sigdt$Gene)
## exir aborts on any NA and on p-values outside [0,1]
stopifnot(!any(is.na(Diff_data)))
stopifnot(min(Diff_data$adj.P.Val) >= 0, max(Diff_data$adj.P.Val) <= 1)
say("Diff_data: ", nrow(Diff_data), " rows x ", ncol(Diff_data), " cols  ",
    "Diff_value = 1 (logFC), Sig_value = 2 (adj.P.Val)")

saveRDS(list(Diff_data = Diff_data, condition = condition, samples = shared),
        file.path(ROOT, "data", "exir_inputs_meta.rds"))

## ------------------------------------------------- 5. RUN exir()
say("--- calling exir() ---")
t0 <- Sys.time()
set.seed(1234)
res <- exir(
  Desired_list           = NULL,          # unrestricted: drivers may be ANY feature,
                                          # so that 'drivers the FFL analysis missed'
                                          # is answerable
  Diff_data              = Diff_data,
  Diff_value             = 1,
  Regr_value             = NULL,
  Sig_value              = 2,
  Exptl_data             = expr,
  Exptl_data_type        = "bulk",
  condition              = condition,
  Exptl_data_orientation = "features_rows",
  normalize              = FALSE,         # Xena values are already log2(x+1)
  pseudo_sample          = FALSE,         # use all 1156 real samples
  Exptl_data_size_check  = TRUE,
  feature_filter         = TRUE,
  min_feature_variance   = 1e-12,
  always_keep_diff_features = TRUE,
  cor_thresh_method      = "mr",
  r                      = 0.5,
  mr                     = 20,
  max.connections        = 50000,
  alpha                  = 0.05,
  num_trees              = 500,
  mtry                   = NULL,
  num_permutations       = 50,
  inf_const              = 1e10,
  ncores                 = 15,
  seed                   = 1234,
  verbose                = TRUE
)
t1 <- Sys.time()
say("exir() finished in ", round(as.numeric(difftime(t1, t0, units = "mins")), 2), " min")

saveRDS(res, file.path(ROOT, "data", "exir_result.rds"))
say("saved ", file.path(ROOT, "data", "exir_result.rds"))

for (nm in names(res)) {
  if (inherits(res[[nm]], "igraph")) {
    say("  ", nm, ": igraph  V=", vcount(res[[nm]]), " E=", ecount(res[[nm]]))
  } else {
    say("  ", nm, ": ", nrow(res[[nm]]), " rows  cols=[",
        paste(colnames(res[[nm]]), collapse = ","), "]")
  }
}
say("=== 08 DONE ===")
