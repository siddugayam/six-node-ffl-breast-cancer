suppressMessages(library(data.table))
B <- "/path/to/revision"
s <- data.table(
 layer=c("TF_target(TRRUST)","TF_target(TRRUST)","TF_target(TRRUST)","TF_target(TRRUST)",
         "TF_miRNA(TransmiR)","TF_miRNA(TransmiR)","TF_miRNA(TransmiR)","TF_miRNA(TransmiR)",
         "miRNA_miRNA(polycistron)","miRNA_miRNA(polycistron)","miRNA_miRNA(polycistron)",
         "gene_gene(TRRUST)","gene_gene(TRRUST_strict_GeneGene)","gene_gene(STRING_physical>=900)","gene_gene(STRING_functional>=900)",
         "miRNA_target(evidence)","miRNA_target(evidence)","miRNA_target(evidence)"),
 metric=c("Activation","Repression","Unknown","Ambiguous",
          "Activation","Repression","Regulation_unsigned","Ambiguous",
          "pairs_3kb","pairs_10kb","pairs_50kb",
          "directed_edges_onto_network_genes","directed_Gene_to_Gene","undirected_pairs","undirected_pairs",
          "strong","weak","predicted_only"),
 n=c(474,292,568,57, 804,424,102,43, 212,480,495, 697,16,446,941, 707,2309,1803))
fwrite(s, file.path(B,"results","layer_summary_counts.tsv"), sep="\t", quote=FALSE)
cat("rows:", nrow(s), "\n")
