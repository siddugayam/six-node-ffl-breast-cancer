## Exploration: node universe, DE coverage, network connectivity
suppressMessages({library(igraph); library(data.table)})
ROOT <- "/path/to/revision"

E <- fread(file.path(ROOT,"data/canonical_edges.tsv"))
cat("edges:", nrow(E), "\n")
nodes <- unique(rbind(E[,.(name=source,type=source_type)], E[,.(name=target,type=target_type)]))
nodes <- unique(nodes, by="name")
cat("nodes:", nrow(nodes), "\n")
print(table(nodes$type))

## self loops / duplicated undirected pairs
E[, `:=`(a=pmin(source,target), b=pmax(source,target))]
cat("self loops:", sum(E$source==E$target), "\n")
cat("unique undirected pairs:", uniqueN(E[,paste(a,b)]), "\n")

g <- graph_from_data_frame(E[source!=target, .(source,target)], directed=FALSE, vertices=nodes)
g <- simplify(g)
cat("simplified undirected: V=", vcount(g), " E=", ecount(g), "\n")
comp <- components(g)
cat("components:", comp$no, " sizes:", paste(sort(comp$csize,decreasing=TRUE)[1:min(10,comp$no)],collapse=","), "\n")

DEg <- fread(file.path(ROOT,"results/v2/v2_DE_genes.csv"))
DEm <- fread(file.path(ROOT,"results/v2/v2_DE_mirnas.csv"))
cat("DE genes:", nrow(DEg), " DE mirnas:", nrow(DEm), "\n")

ng <- nodes[type!="miRNA"]; nm <- nodes[type=="miRNA"]
cat("gene/TF nodes:", nrow(ng), " with gene p:", sum(ng$name %in% DEg$feature), "\n")
cat("miRNA nodes:", nrow(nm), " with mir p:", sum(nm$name %in% DEm$feature), "\n")
cat("miRNA nodes ALSO in gene DE table:", sum(nm$name %in% DEg$feature), "\n")
miss_g <- setdiff(ng$name, DEg$feature); cat("gene/TF missing:", length(miss_g), paste(head(miss_g,30),collapse=","), "\n")
miss_m <- setdiff(nm$name, DEm$feature); cat("miRNA missing:", length(miss_m), paste(head(miss_m,30),collapse=","), "\n")

## p-value distribution sanity
cat("\ngenome-wide gene p: min", min(DEg$P.Value), " frac<0.05:", mean(DEg$P.Value<0.05), "\n")
cat("network gene/TF p: frac<0.05:", mean(DEg[feature %in% ng$name]$P.Value<0.05), " n=", sum(DEg$feature %in% ng$name), "\n")
cat("miRNA p: frac<0.05:", mean(DEm$P.Value<0.05), "\n")
cat("network miRNA p: frac<0.05:", mean(DEm[feature %in% nm$name]$P.Value<0.05), " n=", sum(DEm$feature %in% nm$name), "\n")
cat("n p-values exactly 0 (genes):", sum(DEg$P.Value==0), " (mirna):", sum(DEm$P.Value==0), "\n")

## the paper's key nodes
key <- c("COL1A1","COL3A1","COL1A2","COL5A2","FN1","ETS1","NFKB1","RELA","SP1","EZH2","hsa-miR-29a","hsa-miR-101")
print(DEg[feature %in% key, .(feature,logFC,P.Value,adj.P.Val)])
print(DEm[feature %in% key, .(feature,logFC,P.Value,adj.P.Val)])
cat("key in network:", paste(key[key %in% nodes$name], collapse=","), "\n")
