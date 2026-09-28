suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db);library(data.table);library(AnnotationDbi)})
REV<-"/path/to/revision"
de<-fread(file.path(REV,"results/BRCA_DEX_genes.csv")); nodes<-fread(file.path(REV,"data/canonical_nodes.tsv"))
s2e<-function(x){s<-suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=unique(x),keytype="SYMBOL",columns="ENTREZID"));s<-s[!is.na(s$ENTREZID),];unique(s$ENTREZID)}
U<-s2e(unique(c(de$feature[grepl("^[A-Za-z]",de$feature)], nodes[type%in%c("TF","Gene"),name])))
g6<-s2e(c("COL1A1","COL3A1","NFKB1","RELA","SP1"))
cat("universe",length(U),"input",length(g6),"\n")
for(o in c("CC","BP")){
  e<-enrichGO(g6,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont=o,universe=U,pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.2,minGSSize=10,maxGSSize=500,readable=TRUE)
  d<-as.data.frame(e); cat("\n== GO",o,"significant terms:",nrow(d),"; with Count<3:",sum(d$Count<3),"; with Count==1:",sum(d$Count==1),"\n")
  print(head(d[,c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count","geneID")],10))
}
k<-enrichKEGG(g6,organism="hsa",universe=U,pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.2,minGSSize=10,maxGSSize=500)
k<-setReadable(k,org.Hs.eg.db,"ENTREZID"); d<-as.data.frame(k)
cat("\n== KEGG significant terms:",nrow(d),"; with Count<3:",sum(d$Count<3),"\n")
print(head(d[,c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count","geneID")],12))
