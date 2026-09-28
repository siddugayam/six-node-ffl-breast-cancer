#!/bin/bash
# Build sandbox_<mode>/ mirroring the paths scripts/14_module_detection/03_igraph_battery.R expects, so that the
# script runs unchanged except for its setwd() line. The project itself is only read.
set -e
REV=/path/to/revision
HERE=$(cd "$(dirname "$0")" && pwd)
for m in full nolegacy; do
  S=$HERE/sandbox_$m; mkdir -p $S/scripts/14_module_detection $S/results/v7 $S/data
  cp $REV/scripts/14_module_detection/01_cm2_mcl.R $S/scripts/14_module_detection/
  sed "s#^setwd(\"/path/to/revision\")#setwd(\"$S\")#" $REV/scripts/14_module_detection/03_igraph_battery.R > $S/scripts/14_module_detection/03_igraph_battery.R
  grep -q "setwd(\"$S\")" $S/scripts/14_module_detection/03_igraph_battery.R
  cp $REV/data/canonical_nodes.tsv $S/data/
  Rscript -e "
    suppressMessages({library(igraph); library(data.table)})
    e <- fread('$REV/data/canonical_edges.tsv'); n <- fread('$REV/data/canonical_nodes.tsv')
    if ('$m' == 'nolegacy') { stopifnot(sum(e\$edge_type == 'miRNA_miRNA') == 30); e <- e[edge_type != 'miRNA_miRNA'] }
    gd <- graph_from_data_frame(e[, .(source, target)], directed = TRUE, vertices = n[, .(name, type)])
    gu <- simplify(as_undirected(gd, mode = 'collapse'), remove.multiple = TRUE, remove.loops = TRUE)
    E(gu)\$weight <- 1; stopifnot(components(gu)\$no == 1)
    if ('$m' == 'full') { g0 <- readRDS('$REV/results/v7/cm2_graph_undirected.rds')
      stopifnot(identical(as_edgelist(g0), as_edgelist(gu)), identical(V(g0)\$name, V(gu)\$name)) }
    saveRDS(gu, '$S/results/v7/cm2_graph_undirected.rds')
    cat('$m:', vcount(gu), 'nodes', ecount(gu), 'undirected edges\n')"
done
