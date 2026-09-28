#!/usr/bin/env Rscript
## 20_deconv_install.R -- install every deconvolution package we intend to use.
## Each install is wrapped so one failure does not abort the rest; the final block
## prints an honest availability table that downstream scripts key off.
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8, timeout = 3600)
BASE <- "/path/to/revision"
try_inst <- function(label, expr) {
  cat("\n#### INSTALL", label, "####\n")
  ok <- tryCatch({ eval(expr); TRUE }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  cat("#### ", label, " -> ", ok, "\n", sep = "")
}
try_inst("CRAN deps", quote(install.packages(c("nnls","quadprog","snowfall","corrplot","pheatmap","Rfast","gtools"))))
try_inst("Bioc deps", quote(BiocManager::install(c("GSVA","preprocessCore","Biobase","SummarizedExperiment"), ask = FALSE, update = FALSE)))
try_inst("xCell",      quote(remotes::install_github("dviraran/xCell", upgrade = "never")))
try_inst("MCPcounter", quote(remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source", upgrade = "never")))
try_inst("EPIC",       quote(remotes::install_github("GfellerLab/EPIC", upgrade = "never")))
try_inst("ConsensusTME", quote(remotes::install_github("cansysbio/ConsensusTME", upgrade = "never")))
try_inst("quantiseqr", quote(BiocManager::install("quantiseqr", ask = FALSE, update = FALSE)))
try_inst("InstaPrism", quote(remotes::install_github("humengying0907/InstaPrism", upgrade = "never")))
try_inst("immunedeconv", quote(remotes::install_github("omnideconv/immunedeconv", upgrade = "never")))
pk <- c("xCell","MCPcounter","EPIC","ConsensusTME","quantiseqr","InstaPrism","immunedeconv",
        "estimate","GSVA","e1071","nnls","mediation","preprocessCore")
ip <- rownames(installed.packages())
cat("\n===== AVAILABILITY =====\n")
for (p in pk) cat(sprintf("%-14s %s\n", p, p %in% ip))
