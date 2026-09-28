#!/usr/bin/env Rscript
# Installs the R packages of r_packages.tsv at their recorded versions (R 4.5.1, Bioconductor 3.21).
#   Rscript install_r_packages.R                   install
#   Rscript install_r_packages.R --with-original   also the two packages of the original submission code whose
#                                                   versions were not recorded (current CRAN versions)
#   Rscript install_r_packages.R --check           install nothing; compare the installed versions with the list
args <- commandArgs(TRUE)
options(timeout = 600)
if (as.character(getRversion()) != "4.5.1") message("note: the analysis used R 4.5.1; this is R ", getRversion())

# Bioconductor 3.21 (software, annotation and experiment-data packages; the release fixes the versions)
bioc <- c("AnnotationDbi", "Biobase", "BiocParallel", "BioNet", "clariomshumantranscriptcluster.db", "clusterProfiler", "DESeq2", "DOSE", "edgeR", "enrichplot", "EpiDISH", "fgsea", "GO.db", "graph", "GSEABase", "GSVA", "hgu133plus2.db", "hugene20sttranscriptcluster.db", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "limma", "microRNAome", "minfi", "multiMiR", "oligo", "org.Hs.eg.db", "pd.clariom.s.human", "pd.hg.u133.plus.2", "pd.hg.u133a.2", "pd.hta.2.0", "pd.huex.1.0.st.v2", "pd.hugene.1.0.st.v1", "pd.hugene.2.0.st", "pd.hugene.2.1.st", "pd.mogene.1.0.st.v1", "pd.mogene.1.1.st.v1", "pd.mouse430.2", "pd.mouse430a.2", "preprocessCore", "quantiseqr", "reactome.db", "ReactomePA", "SummarizedExperiment")
# CRAN at the recorded versions; remotes::install_version takes versions that are no longer current from the CRAN archive
cran <- c("cowplot" = "1.2.0", "data.table" = "1.18.2.1", "dplyr" = "1.2.0", "dtangle" = "2.0.10", "e1071" = "1.7-17", "forcats" = "1.0.1", "future" = "1.70.0", "ggplot2" = "4.0.2", "ggpubr" = "0.6.3", "ggraph" = "2.2.2", "ggrepel" = "0.9.7", "ggridges" = "0.5.7", "ggtext" = "0.1.2", "igraph" = "2.2.2", "influential" = "2.3.2", "jsonlite" = "2.0.0", "limSolve" = "2.0.3", "matrixStats" = "1.5.0", "mclust" = "6.1.3", "mediation" = "4.5.1", "metafor" = "5.0-1", "msigdbr" = "26.1.0", "nnls" = "1.6", "patchwork" = "1.3.2", "pheatmap" = "1.0.13", "poweRlaw" = "1.0.0", "ppcor" = "1.1", "ProNet" = "1.0.0", "quadprog" = "1.5-8", "ragg" = "1.5.1", "Rcpp" = "1.1.1", "readxl" = "1.4.5", "scales" = "1.4.0", "stringr" = "1.6.0", "survival" = "3.8-6", "survminer" = "0.5.2", "tibble" = "3.3.1", "tidygraph" = "1.3.1", "tidyr" = "1.3.2", "tidyverse" = "2.0.0")
# GitHub at the recorded commits: package = c(user/repo, commit, subdirectory)
github <- list("BayesPrism" = c("Danko-Lab/BayesPrism", "19052052a6f30833b27e2148294459b0ba2f923e", "BayesPrism"), "CellChat" = c("jinworks/CellChat", "75253cd0c9e68410e6e721a6d3a0419a1d7e358f", ""), "ConsensusTME" = c("cansysbio/ConsensusTME", "6e14ba3d09e39d48b5b9a6b5fa9cd10296cce06a", ""), "EPIC" = c("GfellerLab/EPIC", "50a4f404f96c2842b2891b517b4e3bfaa6c64b8f", ""), "immunedeconv" = c("omnideconv/immunedeconv", "e625e6c28ed14a30f9f40f159925cc9f0df4fa49", ""), "InstaPrism" = c("humengying0907/InstaPrism", "7d3b57cde9342c5eac7282f4c2d8bff09f21e693", ""), "MCPcounter" = c("ebecht/MCPcounter", "b6eac73e91c246fcff0bb1a5c68a816cd588fc48", "Source"), "xCell" = c("dviraran/xCell", "20e2919eefd37e15af35f29f4944e30697098a28", ""))
# R-Forge
rforge <- c("estimate" = "1.0.13")
# installed with R 4.5.1 (recommended packages): MASS 7.3-65, Matrix 1.7-4
original_not_recorded <- c("factoextra", "ggupset")
expected <- c("AnnotationDbi" = "1.70.0", "BayesPrism" = "2.2.3", "Biobase" = "2.68.0", "BiocParallel" = "1.42.2", "BioNet" = "1.68.0", "CellChat" = "2.2.0.9001", "clariomshumantranscriptcluster.db" = "8.8.0", "clusterProfiler" = "4.16.0", "ConsensusTME" = "0.0.1.9000", "cowplot" = "1.2.0", "data.table" = "1.18.2.1", "DESeq2" = "1.48.2", "DOSE" = "4.2.0", "dplyr" = "1.2.0", "dtangle" = "2.0.10", "e1071" = "1.7-17", "edgeR" = "4.6.3", "enrichplot" = "1.28.4", "EPIC" = "1.1.7", "EpiDISH" = "2.24.0", "estimate" = "1.0.13", "fgsea" = "1.34.2", "forcats" = "1.0.1", "future" = "1.70.0", "ggplot2" = "4.0.2", "ggpubr" = "0.6.3", "ggraph" = "2.2.2", "ggrepel" = "0.9.7", "ggridges" = "0.5.7", "ggtext" = "0.1.2", "GO.db" = "3.21.0", "graph" = "1.86.0", "GSEABase" = "1.70.1", "GSVA" = "2.2.1", "hgu133plus2.db" = "3.13.0", "hugene20sttranscriptcluster.db" = "8.8.0", "igraph" = "2.2.2", "IlluminaHumanMethylation450kanno.ilmn12.hg19" = "0.6.1", "immunedeconv" = "2.1.4", "influential" = "2.3.2", "InstaPrism" = "0.1.6", "jsonlite" = "2.0.0", "limma" = "3.64.3", "limSolve" = "2.0.3", "MASS" = "7.3-65", "Matrix" = "1.7-4", "matrixStats" = "1.5.0", "mclust" = "6.1.3", "MCPcounter" = "1.2.0", "mediation" = "4.5.1", "metafor" = "5.0-1", "microRNAome" = "1.30.0", "minfi" = "1.54.1", "msigdbr" = "26.1.0", "multiMiR" = "1.30.0", "nnls" = "1.6", "oligo" = "1.72.0", "org.Hs.eg.db" = "3.21.0", "patchwork" = "1.3.2", "pd.clariom.s.human" = "3.14.1", "pd.hg.u133.plus.2" = "3.12.0", "pd.hg.u133a.2" = "3.12.0", "pd.hta.2.0" = "3.12.2", "pd.huex.1.0.st.v2" = "3.14.1", "pd.hugene.1.0.st.v1" = "3.14.1", "pd.hugene.2.0.st" = "3.14.1", "pd.hugene.2.1.st" = "3.14.1", "pd.mogene.1.0.st.v1" = "3.14.1", "pd.mogene.1.1.st.v1" = "3.14.1", "pd.mouse430.2" = "3.12.0", "pd.mouse430a.2" = "3.12.0", "pheatmap" = "1.0.13", "poweRlaw" = "1.0.0", "ppcor" = "1.1", "preprocessCore" = "1.70.0", "ProNet" = "1.0.0", "quadprog" = "1.5-8", "quantiseqr" = "1.16.0", "ragg" = "1.5.1", "Rcpp" = "1.1.1", "reactome.db" = "1.92.0", "ReactomePA" = "1.52.0", "readxl" = "1.4.5", "scales" = "1.4.0", "stringr" = "1.6.0", "SummarizedExperiment" = "1.38.1", "survival" = "3.8-6", "survminer" = "0.5.2", "tibble" = "3.3.1", "tidygraph" = "1.3.1", "tidyr" = "1.3.2", "tidyverse" = "2.0.0", "xCell" = "1.1.0")

check <- function() {
  have <- vapply(names(expected), function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not installed"), "")
  bad <- have != vapply(expected, function(v) as.character(package_version(v)), "")
  cat(sum(!bad), "of", length(expected), "packages at the recorded version\n")
  if (any(bad)) print(data.frame(package = names(expected)[bad], recorded = expected[bad], installed = have[bad]), row.names = FALSE)
  invisible(!any(bad))
}
if ("--check" %in% args) { check(); quit(save = "no") }

for (p in c("BiocManager", "remotes")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
options(repos = BiocManager::repositories(version = "3.21"))   # CRAN plus Bioconductor 3.21, for dependencies
for (p in names(cran)) {
  have <- tryCatch(as.character(packageVersion(p)), error = function(e) "")
  if (have != as.character(package_version(cran[[p]]))) remotes::install_version(p, version = cran[[p]], upgrade = "never")
}
BiocManager::install(bioc, version = "3.21", ask = FALSE, update = FALSE)
for (p in names(github)) {
  g <- github[[p]]
  remotes::install_github(g[1], ref = g[2], subdir = if (nzchar(g[3])) g[3] else NULL, upgrade = "never")
}
for (p in names(rforge)) install.packages(p, repos = c("https://R-Forge.R-project.org", getOption("repos")))
if ("--with-original" %in% args) install.packages(original_not_recorded)
check()
