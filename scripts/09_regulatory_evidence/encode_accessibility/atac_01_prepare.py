#!/usr/bin/env python3
"""Collapse TCGA-BRCA ATAC-seq (Corces 2018) technical replicates to sample level.

Input : data/atac/BRCA_log2norm.txt  (BRCA-specific peak set, log2 normalised counts,
        GDC ATAC-seq AWG publication page)
        data/atac/BRCA_peakCalls.txt (peak coords, GC, ChIPseeker annotation)
        data/atac/tcga_identifier_lookup.txt (Stanford UUID -> TCGA aliquot barcode)
Output: cache/v6/atac/brca_peaks.tsv.gz, brca_sample_matrix.npy, brca_samples.tsv
"""
import os, sys, numpy as np, pandas as pd

ROOT = "/path/to/revision"
A = os.path.join(ROOT, "data/atac")
OUT = os.path.join(ROOT, "cache/v6/atac")
os.makedirs(OUT, exist_ok=True)

print("[1] peak calls", flush=True)
pk = pd.read_csv(os.path.join(A, "BRCA_peakCalls.txt"), sep="\t")
print("   peaks:", pk.shape, flush=True)

print("[2] counts matrix (large)", flush=True)
cnt = pd.read_csv(os.path.join(A, "BRCA_log2norm.txt"), sep="\t")
print("   counts:", cnt.shape, flush=True)

meta = cnt[["seqnames", "start", "end", "name", "score"]].copy()
reps = [c for c in cnt.columns if c.startswith("BRCA_")]
print("   tech reps:", len(reps), flush=True)

# map rep -> stanford UUID (sample) -> TCGA aliquot barcode
lk = pd.read_csv(os.path.join(A, "tcga_identifier_lookup.txt"), sep="\t")
lk = lk[lk.bam_prefix.str.startswith("BRCA")].copy()
lk["rep_key"] = lk.bam_prefix.str.replace("-", "_", regex=False)
rep2uuid = dict(zip(lk.rep_key, lk.stanfordUUID))
rep2barcode = dict(zip(lk.rep_key, lk.aliquot_id))
missing = [r for r in reps if r not in rep2uuid]
print("   reps missing from lookup:", len(missing), flush=True)

X = cnt[reps].to_numpy(dtype=np.float32)
uuids = np.array([rep2uuid.get(r, "NA") for r in reps])
usamp = sorted(set(uuids[uuids != "NA"]))
S = np.zeros((X.shape[0], len(usamp)), dtype=np.float32)
nrep = []
for j, u in enumerate(usamp):
    idx = np.where(uuids == u)[0]
    S[:, j] = X[:, idx].mean(axis=1)
    nrep.append(len(idx))

samp = pd.DataFrame({"stanfordUUID": usamp, "n_techrep": nrep})
bc = {}
for r in reps:
    u = rep2uuid.get(r)
    if u:
        bc[u] = rep2barcode.get(r)
samp["aliquot_barcode"] = samp.stanfordUUID.map(bc)
samp["tcga_sample"] = samp.aliquot_barcode.str.slice(0, 15)
samp["tcga_patient"] = samp.aliquot_barcode.str.slice(0, 12)

print("[3] writing", flush=True)
meta.to_csv(os.path.join(OUT, "brca_peak_meta.tsv.gz"), sep="\t", index=False)
pk.to_csv(os.path.join(OUT, "brca_peak_calls.tsv.gz"), sep="\t", index=False)
np.save(os.path.join(OUT, "brca_sample_matrix.npy"), S)
samp.to_csv(os.path.join(OUT, "brca_samples.tsv"), sep="\t", index=False)
print("   samples:", samp.shape[0], "patients:", samp.tcga_patient.nunique(), flush=True)
print("   matrix:", S.shape, flush=True)
print("DONE", flush=True)
