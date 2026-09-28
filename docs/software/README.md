# Software requirements

These files list the software the deposited code uses, at the versions of the analysis environment.

## System used

- Ubuntu 22.04.5 LTS, Linux kernel 6.8.0
- Python 3.13.13 (Anaconda base environment)
- R 4.5.1 (2025-06-13), Bioconductor 3.21
- gcc (Ubuntu 10.5.0-1ubuntu1~22.04.3) 10.5.0 and g++ of the same version

## Python: `requirements.txt`

`pip install -r requirements.txt` installs the 24 packages the code imports, at the recorded versions.
- **torch** is the CUDA 12.4 build (2.6.0+cu124). The file points pip to the PyTorch index for it.

## R: `r_packages.tsv` and `install_r_packages.R`

`r_packages.tsv` lists 97 packages with their version and source:
- 42 from Bioconductor 3.21, including 12 array-platform design packages (pd.*) and 2 array
  annotation packages that the perturbation scripts in `analyses/perturbation_tests/` load at run time;
- 27 from CRAN at their current version, and 15 from the CRAN archive, because the recorded
  version is no longer current or the package has left CRAN;
- 8 from GitHub at the recorded commit;
- 1 from R-Forge;
- 2 installed with R itself;
- 2 installers: BiocManager and remotes;
- 2 with no recorded version: factoextra and ggupset.

The last two are used only by the original submission's `03_Network_Analysis_Hub_Enrichment.R`, which installs
them itself. They were never installed on the analysis machine.

The `used_by` column says whether a package serves the main analyses (`scripts/` and `analyses/`), the original submission code (`analyses/original_submission_code/`), or both.

`install_r_packages.R` installs every listed version:
- `Rscript install_r_packages.R` installs them;
- `--with-original` adds the two packages without a recorded version;
- `--check` installs nothing and compares the installed versions with the list.

**preprocessCore.** On the analysis machine it was built with `--disable-threading`, because oligo's RMA otherwise
failed with `pthread_create` error 22. If that happens, reinstall it:
`BiocManager::install("preprocessCore", configure.args = "--disable-threading", force = TRUE)`.

## External programs

The code calls these programs:

| program | used for | where |
|---|---|---|
| Rscript | R steps run from Python and shell scripts | many scripts in `scripts/` and `analyses/` |
| python3 | Python steps run from shell scripts | shell drivers |
| gcc | C programs: census enumerators, null models (`-O2` or `-O3 -march=native`, `-lm`) | `scripts/03_ffl_census/`, `scripts/04_motif_significance/`, `analyses/census_and_motif_nulls/`, `analyses/six_node_pattern/` |
| g++ | Rcpp compiles the jActiveModules engine `jam_engine.cpp` | `scripts/14_module_detection/jam_*.R` |
| bedtools | interval intersections of ENCODE and ATAC regions | `scripts/09_regulatory_evidence/encode_accessibility/atac_06*.py` |
| bigBedToBed (UCSC utilities) | extracting ReMap 2022 peaks from the bigBed file | `scripts/09_regulatory_evidence/chip_remap_targetscan/v2d_remap_*.py`, `scripts/09_regulatory_evidence/sequence_models/25*_seqreg_remap*.py` |
| curl | batch downloads | `v2_chip_download.sh` and `v2_dl_chunked.sh` in `scripts/09_regulatory_evidence/chip_remap_targetscan/` |
| zcat, cat | reading compressed and plain files through pipes | R scripts (`scripts/07_expression_validation/external_cohorts/N0_geo_utils.R`, `scripts/10_survival_and_clinical/multi_cohort/M12_cohort_inventory.R`), shell scripts, `scripts/09_regulatory_evidence/encode_accessibility/atac_06*.py` |
| awk, sed, cut, split, sort, seq, wc, xargs | text handling in shell scripts | shell scripts |
| nvidia-smi, df, du, ps | resource logging only | `analyses/perturbation_tests/P0/` and `P4/` |

## GPU and hardware

The Enformer and Borzoi steps need an NVIDIA GPU with CUDA and native bfloat16 support (Ampere or newer). They are
`scripts/09_regulatory_evidence/sequence_models/22_seqreg_enformer.py`, `22b_seqreg_enformer_ism.py`, `23_seqreg_borzoi.py`, `41_seqreg_ext_occlusion_null.py`,
`45_seqreg_ext_splicedonor_control.py` (all in that folder), and the Borzoi occlusion `analyses/perturbation_tests/P4/p4_borzoi_occlusion.py`:
- They call `.cuda()` directly and run under bfloat16 autocast, so they do not run on a CPU without code changes.
- They ran on an NVIDIA RTX A4000 (16 GB, driver 535.288.01).
- The Borzoi occlusion used about 3.8 GB of GPU memory at batch size 2 and took 0.29 s per sequence.
- The model weights are downloaded from Hugging Face on first use, so network access is needed.

Host memory for these steps was not recorded. The analysis machine has 125 GB of RAM and 16 CPU cores.
