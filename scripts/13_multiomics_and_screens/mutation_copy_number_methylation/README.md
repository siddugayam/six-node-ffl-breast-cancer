# Mutation copy number methylation

## Scripts in this folder

| Script | What it does |
|---|---|
| `mo_01_mutation.R` | Somatic non-silent mutation frequency in TCGA-BRCA for network genes/TFs |
| `mo_02_copynumber.R` | GISTIC2 thresholded copy number for network nodes + CN-expression correlation |
| `mo_03_probeselect.py` | Selects the 450k methylation probes for the promoter analysis of the network genes and of a background gene set, and writes the miRNA loci (hg19) to cache/multiomics/. |
| `mo_04_probemap.R` | Maps Illumina 450k probes to gene promoters (IlluminaHumanMethylation450kanno.ilmn12.hg19) for the network genes and the background genes; writes cache/multiomics/promoter_probes_*.tsv. |
| `mo_05_methylation.R` | TCGA-BRCA Illumina 450k promoter methylation of network nodes |
| `mo_06_mut_lengthmatched.R` | Gene-length-matched control for the mutation frequency of network genes. |
| `mo_07_integrate.R` | Integrated multi-omic summary for the featured hubs. |
