# ChIP-seq, ReMap, ENCODE and TargetScan evidence for network edges

Two passes: `01`–`06` (ChIP-Atlas, ReMap promoters, ENCODE counts and the TargetScan download) and `08`–`23`, which
wrote the deposited tables (`results/v2/chip_evidence.csv`, `targetscan_sites.csv` and the summaries). The branches
for ChIP-Atlas, ReMap, ENCODE and TargetScan are independent of one another.

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_chipatlas_download.sh` | Download ChIP-Atlas hg38 "Target Genes" tables for a list of TF symbols at a given TSS-distance threshold (1, 5 or 10 kb). |
| `02_chipatlas_edge_support.py` | ChIP-seq evidence for the TF->target edges of the breast-cancer miRNA-TF FFL network. |
| `03_chipatlas_edge_null.py` | Is the ChIP-seq support of the 770 published TF->target edges better than chance? |
| `04_remap_collagen_promoters.py` | Reverse query: which TFs are actually bound at the COL1A1 and COL3A1 promoters? |
| `05_encode_experiment_counts.py` | Third resource: ENCODE portal REST API - how many human TF ChIP-seq experiments exist for NFKB1 / RELA / SP1 / ETS1, and in which biosamples (breast highlighted). |
| `06_targetscan_download.sh` | N-way parallel HTTP range download (targetscan.org throttles single connections hard). |
| `08_node_tss_table.py` | Build hg38 TSS table for every gene/TF node of the canonical network from UCSC ncbiRefSeqSelect (one MANE/representative transcript per gene). |
| `09_chipatlas_fetch.py` | Download ChIP-Atlas 'Target Genes' tables (hg38) for every TF that is a source of a TF_target edge in the canonical network. |
| `10_chipatlas_edge_evidence.py` | ChIP-seq support for the canonical network's 770 TF_target edges, from ChIP-Atlas hg38 'Target Genes' tables (average and per-experiment MACS2 -10*log10(p) scores in TSS +/- 1 kb and +/- 10 kb windows). |
| `11_remap_fetch.py` | Query ReMap 2022 (hg38, all non-merged datasets) for every ChIP-seq peak whose summit falls within TSS +/-5 kb of each gene/TF node of the canonical network. |
| `12_remap_edge_support.py` | ReMap 2022 (hg38, all 8,103 non-merged datasets) promoter occupancy. |
| `13_remap_background_fetch.py` | Genome-wide background: ReMap 2022 promoter occupancy at 600 randomly chosen protein-coding RefSeqSelect genes (excluding the network's own nodes). |
| `14_remap_genomewide_null.py` | Genome-wide promoter null for ReMap edge support: for each network TF_target edge, replace the real target promoter by a randomly drawn promoter from 600 random RefSeqSelect genes that are NOT network nodes. |
| `15_collagen_promoter_celltypes.py` | Which factors occupy the COL1A1 / COL3A1 promoters in FIBROBLAST, STROMAL or BREAST biosamples specifically? |
| `16_collagen_promoter_specificity.py` | Which TFs are SPECIFICALLY bound at the COL1A1 / COL3A1 promoters? |
| `17_encode_counts_per_tf.py` | ENCODE portal: number of released human TF ChIP-seq experiments per network TF, read from the target.label facet of a single faceted search. |
| `18_chip_evidence_merged.py` | Merge the three ChIP resources into one per-edge evidence table for Table 1, and produce the headline coverage fractions. |
| `20_targetscan_sites_all.py` | TargetScan 8.0 sequence-level evidence for the canonical network's miRNA_target edges. |
| `21_targetscan_site_null.py` | Are the network's 4,819 miRNA-target edges enriched for TargetScan sites over random miRNA-gene pairs drawn from the same miRNA and gene universes? |
| `22_targetscan_vs_edge_rho.R` | Does TargetScan site quality predict which network miRNA-target edges are sign-concordant (rho < 0) in TCGA-BRCA? |
| `23_direct_evidence_summary.py` | Summary table of the direct molecular evidence (ChIP-seq, ReMap and TargetScan) for the network edges. |
