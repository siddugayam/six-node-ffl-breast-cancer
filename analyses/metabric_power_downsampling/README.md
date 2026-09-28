# METABRIC downsampled to the TCGA event count (Results 3.8, Methods 2.6)

`ct_11_power_downsample_rerun.R` downsamples METABRIC to the TCGA event count 200 times and counts the prognostic hubs
at *q* < 0.05 in each subsample. It reads patient-level METABRIC data, which are not redistributed.
