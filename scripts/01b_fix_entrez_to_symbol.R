# Repair: the pan-cancer gene matrix is keyed by Entrez Gene ID, not HGNC symbol.
# Re-key brca_gene_expr.rds and the gene DE table to symbols so they join to the network.
suppressMessages({library(data.table); library(org.Hs.eg.db); library(AnnotationDbi); library(limma)})
REV <- "/path/to/revision"
lg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

expr  <- readRDS(file.path(REV,"data/brca_gene_expr.rds"))
pheno <- readRDS(file.path(REV,"data/brca_pheno.rds"))
lg("expr dims", paste(dim(expr), collapse=" x "))

ent <- rownames(expr)
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ent, keytype="ENTREZID", columns="SYMBOL"))
map <- map[!is.na(map$SYMBOL), ]
map <- map[!duplicated(map$ENTREZID), ]
lg("entrez ids:", length(ent), " mapped to symbol:", nrow(map))

expr <- expr[map$ENTREZID, , drop=FALSE]
sym  <- map$SYMBOL
# collapse duplicate symbols by highest mean expression
o    <- order(sym, -rowMeans(expr, na.rm=TRUE))
keep <- !duplicated(sym[o]); idx <- o[keep]
expr <- expr[idx, , drop=FALSE]; rownames(expr) <- sym[idx]
expr <- expr[order(rownames(expr)), , drop=FALSE]
lg("symbol-keyed expr dims", paste(dim(expr), collapse=" x "))
saveRDS(expr, file.path(REV,"data/brca_gene_expr_symbol.rds"))

# --- re-run limma DE on symbol-keyed matrix
st <- pheno$sample_type[match(colnames(expr), pheno$sample)]
grp <- factor(ifelse(st=="Primary Tumor","Tumor", ifelse(st=="Solid Tissue Normal","Normal", NA)))
ok  <- !is.na(grp); expr2 <- expr[, ok, drop=FALSE]; grp <- droplevels(grp[ok])
lg("DE groups:", paste(names(table(grp)), table(grp), collapse="  "))
keepr <- apply(expr2, 1, function(v) sum(!is.na(v))>=3 && sd(v, na.rm=TRUE)>0)
expr2 <- expr2[keepr, , drop=FALSE]
des <- model.matrix(~0+grp); colnames(des) <- levels(grp)
fit <- eBayes(contrasts.fit(lmFit(expr2, des), makeContrasts(Tumor-Normal, levels=des)))
tt  <- topTable(fit, number=Inf, sort.by="P")
tt$feature <- rownames(tt)
tt <- as.data.table(tt)[, .(feature, logFC, AveExpr, t, P.Value, adj.P.Val)]
fwrite(tt, file.path(REV,"results/BRCA_DEX_genes.csv"))
lg("genes tested:", nrow(tt), " sig |logFC|>1 & adj.P<0.05:", sum(abs(tt$logFC)>1 & tt$adj.P.Val<0.05))

mir <- fread(file.path(REV,"results/BRCA_DEX_mirnas.csv"))
all <- rbind(tt[, .(Gene=feature, logFC, adj.P.Val, class="gene")],
             mir[, .(Gene=feature, logFC, adj.P.Val, class="miRNA")])
fwrite(all, file.path(REV,"results/BRCA_DEX_ALL_nodes.csv"))
lg("BRCA_DEX_ALL_nodes.csv rows:", nrow(all))

nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
nodes[, matched := name %in% all$Gene]
print(nodes[, .(n=.N, matched=sum(matched), rate=round(100*mean(matched),1)), by=type])
lg("unmatched genes/TFs:", paste(head(nodes[type!="miRNA" & !matched, name], 20), collapse=", "))
fwrite(nodes, file.path(REV,"results/node_DE_match.tsv"), sep="\t")
lg("DONE")
