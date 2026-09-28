# The pan-cancer gene matrix is keyed by HGNC SYMBOL (a minority of rows are bare Entrez
# IDs for loci with no symbol). Recompute gene-level DE on the native identifiers and
# verify the join to the network node set.
suppressMessages({library(data.table); library(limma)})
REV <- "/path/to/revision"
lg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

expr  <- readRDS(file.path(REV,"data/brca_gene_expr.rds"))
pheno <- readRDS(file.path(REV,"data/brca_pheno.rds"))
lg("expr dims", paste(dim(expr), collapse=" x "))
lg("rownames that are bare numeric IDs:", sum(grepl("^[0-9]+$", rownames(expr))),
   " symbol-like:", sum(!grepl("^[0-9]+$", rownames(expr))))
saveRDS(expr, file.path(REV,"data/brca_gene_expr_symbol.rds"))   # canonical name for downstream

st  <- pheno$sample_type[match(colnames(expr), pheno$sample)]
grp <- factor(ifelse(st=="Primary Tumor","Tumor", ifelse(st=="Solid Tissue Normal","Normal", NA)))
ok  <- !is.na(grp); e2 <- expr[, ok, drop=FALSE]; grp <- droplevels(grp[ok])
lg("DE groups:", paste(paste0(names(table(grp)),"=",table(grp)), collapse="  "))

keepr <- apply(e2, 1, function(v) sum(!is.na(v))>=3 && sd(v, na.rm=TRUE)>0)
e2 <- e2[keepr, , drop=FALSE]
des <- model.matrix(~0+grp); colnames(des) <- levels(grp)
fit <- eBayes(contrasts.fit(lmFit(e2, des), makeContrasts(Tumor-Normal, levels=des)))
tt  <- as.data.table(topTable(fit, number=Inf, sort.by="P"), keep.rownames="feature")
tt  <- tt[, .(feature, logFC, AveExpr, t, P.Value, adj.P.Val)]
fwrite(tt, file.path(REV,"results/BRCA_DEX_genes.csv"))
lg("genes tested:", nrow(tt), " sig |logFC|>1 & adj.P<0.05:",
   sum(abs(tt$logFC)>1 & tt$adj.P.Val<0.05))
lg("top 8 by p:", paste(head(tt$feature,8), collapse=", "))

mir <- fread(file.path(REV,"results/BRCA_DEX_mirnas.csv"))
all <- rbind(tt[, .(Gene=feature, logFC, adj.P.Val, class="gene")],
             mir[, .(Gene=feature, logFC, adj.P.Val, class="miRNA")])
fwrite(all, file.path(REV,"results/BRCA_DEX_ALL_nodes.csv"))
lg("BRCA_DEX_ALL_nodes.csv rows:", nrow(all))

nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
nodes[, matched := name %in% all$Gene]
print(nodes[, .(n=.N, matched=sum(matched), rate=round(100*mean(matched),1)), by=type])
lg("unmatched genes/TFs:", paste(nodes[type!="miRNA" & !matched, name], collapse=", "))
lg("unmatched miRNAs:",   paste(nodes[type=="miRNA"  & !matched, name], collapse=", "))
fwrite(nodes, file.path(REV,"results/node_DE_match.tsv"), sep="\t")
lg("DONE")
