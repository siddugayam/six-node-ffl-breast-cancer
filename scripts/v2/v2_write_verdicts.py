#!/usr/bin/env python3
import csv,os
OUT="/path/to/revision/results/v2/verify_network.csv"
rows=[
# claim, published_value, my_value, verdict, note
("N1 raw node labels","803","803","CONFIRMED","union of SIF endpoints (803) == union of node_attributes name columns (803); identical sets, 0 asymmetric labels"),
("N1 canonical nodes after harmonisation","587","587","CONFIRMED","215 labels changed by hsa-mir->hsa-miR; 2 labels (hsa-mir-129-2, hsa-miR-129-2) collapsed onto hsa-miR-129; 803-215-1=587"),
("N1 node type breakdown","223 miRNA / 207 Gene / 157 TF","223 miRNA / 207 Gene / 157 TF","CONFIRMED","TF = labelled TF in >=1 of the six attribute CSVs (157 labels); 10 labels carry conflicting Gene/TF Type across files, TF wins"),
("N1 unique directed edges","6859","6860","DIFFERS (6860)","difference is exactly the FOS->FOS self-loop, present in 3-miR/3-TF/3-Comp.sif. 6859 if self-loops are dropped; the deposited canonical_edges.tsv drops it"),
("N1 miRNA_target edges","4819","4819","CONFIRMED",""),
("N1 TF_miRNA edges","1239","1239","CONFIRMED",""),
("N1 TF_target edges","770","771","DIFFERS (771)","771 includes FOS->FOS; 770 excluding the self-loop"),
("N1 miRNA_miRNA edges","30","30","CONFIRMED",""),
("N1 gene_gene edges","1","1","CONFIRMED","COL1A1->COL3A1"),
("N1 reciprocal TF<->miRNA pairs","1223","1223","CONFIRMED","pairs with both TF->miRNA and miRNA->TF present in the deduplicated union"),
("N2 miRNA labels present in BOTH miRBase cases","214 of 439","214 of 439 (canonical nodes) / 215 of 439 (raw case-groups)","CONFIRMED (definition-dependent)","439 raw miRNA labels; 224 distinct case-insensitive labels, 215 of which occur in both cases; those 215 map onto only 214 canonical nodes because hsa-mir-129 and hsa-mir-129-2 both collapse to hsa-miR-129"),
("N2 Jaccard of two copies' target sets, miR-130a","0.00","0.0000","CONFIRMED","hsa-mir-130a 8 targets vs hsa-miR-130a 50 targets, 0 shared"),
("N2 Jaccard, miR-21","0.00","0.0000","CONFIRMED","27 vs 13 targets, 0 shared"),
("N2 Jaccard, miR-124","0.00","0.0000","CONFIRMED","11 vs 46 targets, 0 shared"),
("N2 Jaccard, miR-29a","0.00","0.0217","DIFFERS (0.0217)","hsa-mir-29a 15 targets vs hsa-miR-29a 32 targets; exactly ONE shared target, NFKB1; union 46 -> J=1/46=0.0217. Not zero, but still effectively disjoint"),
("N2 Jaccard, miR-34a","0.00","0.0000","CONFIRMED","23 vs 9 targets, 0 shared"),
("N3 exhaustive 3-node FFL cores (union network)","6037","not reproduced; 1649 (mutually exclusive classes) / 3083 (composites double-counted)","COULD NOT REPRODUCE","exhaustive enumeration over the 587-node/6860-edge union under 144 definitional variants (regulator-pair types x target-node types x mutual-pair handling x raw/harmonised labels x per-network/union) produced no value within +/-1 of 6037; max union value observed 4637"),
("N3 3-node FFL cores, 3-miR.sif","1593","1593","CONFIRMED","target node restricted to Gene|TF; TF-FFL 0 + miRNA-FFL 1593; composites counted in both TF- and miRNA-FFL classes"),
("N3 3-node FFL cores, 3-TF.sif","1370","1370","CONFIRMED","TF-FFL 1303 + miRNA-FFL 67 (of which 4 composite, counted in both)"),
("N3 3-node FFL cores, 3-Comp.sif","2823","2823","CONFIRMED","TF-FFL 1380 + miRNA-FFL 1443; 1380 composite (all TF->miRNA arcs with a shared target are reciprocated)"),
("N3 4-TF.sif contains zero 3-node cores","0","0","CONFIRMED",""),
("N3 5-TF.sif contains zero 3-node cores","0","0","CONFIRMED",""),
("N3 6-TF.sif contains zero 3-node cores","0","0 (Gene|TF target) / 8 (if a miRNA may be the FFL target)","CONFIRMED (with caveat)","8 TF->miRNA->miRNA triangles exist in 6-TF.sif (e.g. NFKB1->let-7b, NFKB1->miR-29a, let-7b->miR-29a); they are excluded only because the spec confines FFL targets to Gene|TF"),
("N7 TRRUST Activation among 770 TF_target edges","301","301","CONFIRMED","resolution rule: Unknown absorbed by a co-occurring signed mode; Activation+Repression together = Ambiguous"),
("N7 TRRUST Repression among 770 TF_target edges","136 (17.7%)","136 (17.7%)","CONFIRMED",""),
("N7 TRRUST Unknown among 770 TF_target edges","300","300","CONFIRMED",""),
("N7 TRRUST Ambiguous among 770 TF_target edges","33","33","CONFIRMED","under a stricter rule (any pair with >1 distinct TRRUST mode = Ambiguous) the split is 239/126/300/105"),
("N7 TRRUST coverage of TF_target edges","not stated","770/770 = 100%","NEW FINDING","every single published TF_target edge is a TRRUST pair; the TF->gene layer is therefore a TRRUST subset, so TRRUST is not an independent validation source for it"),
("N13 miR-130a distinct target count","58 (manuscript says 109)","58","CONFIRMED","hsa-miR-130a 50 targets + hsa-mir-130a 8 targets, disjoint, 58 after harmonisation; all 58 are Gene/TF"),
("N13 VEGFA a miR-130a target","no","no","CONFIRMED","VEGFA is in the network (degree 75 in 3-Comp) but has no incoming edge from either miR-130a copy"),
("N13 MMP2 a miR-130a target","no","no","CONFIRMED","MMP2 is in the network but not a miR-130a target"),
("N13 COL1A1 a miR-130a target","no","no","CONFIRMED","COL1A1 is in the network but not a miR-130a target"),
("N13 HK2 present in the network","absent","absent","CONFIRMED","HK2 is not among the 803 raw labels nor any endpoint of any SIF edge; the manuscript's let-7b-HK2 axis has no support in the deposited files"),
("N13 miRNAs of manuscript section 3.5 absent from network","8 of 13","7 of 13 absent; 6 present; only 4 of 13 actually have an edge to COL1A1","DIFFERS (7 absent / 9 unsupported)","absent: miR-3619, miR-412, miR-4776, miR-4800, miR-516a, miR-516b, miR-7162. Present but with NO COL1A1 edge: miR-124, miR-92a-1 (as hsa-miR-92a). Present with a COL1A1 edge: miR-133a, miR-133b, miR-143, miR-218"),
("N13 VEGFA degree in 3-Comp","75 (manuscript says 111)","75","CONFIRMED","VEGFA degree is 75 in 3-Comp (raw and harmonised and the deposited attribute CSV); 111 is its degree in 3-TF"),
("EXTRA manuscript degrees attributed to 3-Comp","VEGFA 111, CCND2 88, ADAMTS5 70, TGFBR2 70, TP53 111, MYC 104, RELA 98, SP1 98, miR-130a 50, miR-21 48, miR-124 46, miR-34a 41","all twelve match the 3-TF.csv Degree column exactly; 3-Comp values are VEGFA 75, CCND2 88, ADAMTS5 71, TGFBR2 69, TP53 135, MYC 121, RELA 124, SP1 114","DIFFERS","the manuscript's 'composite framework' degree table is the 3-TF network's attribute table, and it silently mixes miRBase cases: 130a/124 are the hsa-miR- copies, 21/34a are the hsa-mir- copies"),
("EXTRA layer counts, 3-TF network","1554 TF-miRNA / 808 TF-Gene / 3379 miRNA-Gene","1554 / 808 (=722 TF->Gene + 86 TF->TF) / 3379 (=3160 miRNA->Gene + 219 miRNA->TF)","CONFIRMED","these are NON-deduplicated SIF line counts (5741 lines); after dedup they are 1223 / 758 / 3379"),
("G  diff vs revision/data/canonical_edges.tsv","6859 edges","6860 edges; 1 edge in mine only, 0 in theirs only, 0 edge_type disagreements","DIFFERS BY 1","the only difference is the FOS->FOS self-loop, which the deposited pipeline drops"),
("G  diff vs revision/data/canonical_nodes.tsv","587 nodes","587 nodes; 0 nodes unique to either side, 0 type disagreements","CONFIRMED","node set and node typing reproduce exactly from the raw files with an independent implementation"),
]
os.makedirs(os.path.dirname(OUT),exist_ok=True)
with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["claim","published_value","my_value","verdict","note"]); w.writerows(rows)
print("wrote",OUT,len(rows),"rows")
