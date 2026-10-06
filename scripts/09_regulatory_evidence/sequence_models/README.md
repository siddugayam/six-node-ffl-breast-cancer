# Sequence models

## Scripts in this folder

| Script | What it does |
|---|---|
| `20_seqreg_regions.py` | Define hg38 loci for sequence-level regulatory analysis and fetch reference sequence from the UCSC REST API. |
| `21_seqreg_alphagenome.py` | AlphaGenome (Google DeepMind, 2025) sequence-to-function predictions at the COL1A1 / COL3A1 / miR-29 loci. |
| `21b_seqreg_alphagenome_access.py` | Record, reproducibly, exactly what was attempted to obtain AlphaGenome access, what the server replied, and the single step the account holder must perform. |
| `22_seqreg_enformer.py` | Enformer (Avsec et al. 2021) sequence-to-function predictions at the COL1A1 / COL3A1 / miR-29 loci (hg38), with a genome-wide promoter background so that every statement is a calibrated percentile rather than a raw number. |
| `22b_seqreg_enformer_ism.py` | Sequence attribution at the COL1A1 and COL3A1 promoters with Enformer. |
| `22c_seqreg_ism_motifmap.py` | Map the Enformer in-silico-mutagenesis peaks onto the JASPAR/HOCOMOCO motif sites called by 28_, so that the model's most load-bearing promoter positions are named. |
| `23_seqreg_borzoi.py` | Borzoi (Linder et al., Nat Genet 2025) sequence-to-function predictions at the COL1A1/COL3A1/miR-29 loci. |
| `24_seqreg_chipatlas.py` | ChIP-Atlas (Zou et al., NAR 2022) evidence that a TF is bound within +/-1 kb and +/-5 kb of the COL1A1 and COL3A1 TSS, resolved by cell type, with particular attention to fibroblast / mesenchymal cell types. |
| `25_seqreg_remap.py` | ReMap 2022 (Hammal et al., NAR 2022) experimentally mapped TF peaks at the COL1A1 / COL3A1 / miR-29 loci (hg38), with cell-type (biotype) resolution and a random-promoter background so that "TF X binds here" is calibrated. |
| `25b_seqreg_remap_celltypes.py` | Recompute the ReMap 2022 region tables with a CURATED cell-type (biotype) whitelist. |
| `25c_seqreg_remap_background.py` | ReMap 2022 random-promoter background: for 300 random RefSeq-Select promoters, how often is each TF observed with a peak in the same TSS-1000..+500 window? |
| `25d_seqreg_remap_background_counts.py` | As 25c but storing the NUMBER of ReMap peaks per TF per background promoter, so that a TF's peak count at COL1A1/COL3A1 can be given as a percentile of its peak count at 300 random promoters (presence/absence alone is … |
| `26_seqreg_screen_ccre.py` | ENCODE SCREEN candidate cis-regulatory elements (cCREs) and UCSC CpG islands at the COL1A1/COL3A1/miR-29 loci (hg38). |
| `27_seqreg_gtex.py` | GTEx v8/v10 cis-eQTL and sQTL evidence for COL1A1, COL3A1 and the miR-29 host genes, with emphasis on breast and cultured fibroblasts. |
| `27b_seqreg_gtex_extra.py` | GTEx addendum: (a) independent cis-eQTLs for the miRNA genes themselves, (b) an explicit, reproducible record of the breast / cultured-fibroblast NEGATIVES, (c) the positions of every COL1A1/COL3A1 eQTL and fine-mapped variant … |
| `28_seqreg_motifscan.py` | FIMO-equivalent PWM scan of the COL1A1 / COL3A1 (and miRNA-locus) promoters with JASPAR 2024 CORE vertebrates non-redundant and HOCOMOCO v11 core human, calibrated against a GC-matched genome-wide promoter background. |
| `28b_seqreg_motif_calibration.py` | Higher-resolution GC-matched calibration of the focus-TF motif counts at all eight promoters. |
| `28c_seqreg_motif_bh.py` | Benjamini-Hochberg correction of the focus-motif calibration (25 motifs x 8 promoters = 200 tests), plus a note of which called sites are the same genomic site detected by several matrices (they are not independent tests). |
| `29_seqreg_cpg_architecture.py` | Promoter sequence architecture at the COL1A1/COL3A1 and miR-29 loci: GC content, CpG observed/expected, CpG-island overlap (UCSC cpgIslandExt) and ENCODE cCRE class, from sequence, independently of any deep-learning model. |
| `31_seqreg_enhancers.py` | Candidate enhancers of COL1A1 and COL3A1: ENCODE cCREs within +/-100 kb, intersected with ReMap 2022 peaks, reported by TF and by cell-type class, with fibroblast/mesenchymal peaks called out. |
| `32_seqreg_rela_enhancers.py` | The only fibroblast-context TF-binding evidence at the collagen loci in ReMap 2022 is RELA in TNF-stimulated IMR-90 lung fibroblasts (GSE43070) and nutlin-treated human dermal fibroblasts (GSE77225). |
| `33_seqreg_synthesis.py` | Assemble every sequence-level result into one tidy table (results/v3/seqreg_key_results.csv). |
| `34_seqreg_network_tf_edges.py` | Sequence-level tests of every TF → collagen edge of the network. |
| `35_seqreg_nfkb_context.py` | Every NF-kB-family ReMap 2022 peak within +/-100 kb of COL1A1 and COL3A1, annotated by cell type and by the treatment string encoded in the ReMap dataset name (TNF, LPS, IL1, etc.). |
| `36_seqreg_key_results_addendum.py` | Append the remaining headline numbers to results/v3/seqreg_key_results.csv. |
| `40_seqreg_ext_alphagenome_access.py` | Re-verification and EXTENSION of the AlphaGenome access record. |
| `41_seqreg_ext_occlusion_null.py` | EXTENSION of 22b_seqreg_enformer_ism.py. |
| `42_seqreg_ext_chipatlas_fibroblast.py` | EXTENSION of 24_seqreg_chipatlas.py. |
| `43_seqreg_ext_gtex_singletissue.py` | EXTENSION of 27_seqreg_gtex.py. |
| `44_seqreg_ext_synthesis.py` | Assemble the extension results into one tidy table. |
| `45_seqreg_ext_splicedonor_control.py` | Splice-donor control for the occlusion effects at the collagen promoters: the called ETS1 site at COL3A1 TSS+186 straddles the exon-1 splice donor. |
