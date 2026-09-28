#!/usr/bin/env python3
"""21_scrna_reference_build.py
Build a breast-tumour-specific deconvolution reference from Wu et al. 2021 (GSE176078).

Reads the 29,733 x 100,064 sparse count matrix in chunks and accumulates, per
celltype_major (9 types) and per (patient, celltype_major):
  - raw UMI sums                (for BayesPrism/InstaPrism reference phi)
  - mean CPM across cells       (for signature-matrix marker selection)
  - cell counts
Outputs go to cache/v3/deconv/ as .npz + .tsv so R can read them.
"""
import numpy as np, pandas as pd, os, sys, time

BASE = "/path/to/revision"
SC   = os.path.join(BASE, "cache/celltype/Wu_etal_2021_BRCA_scRNASeq")
OUT  = os.path.join(BASE, "cache/v3/deconv")
os.makedirs(OUT, exist_ok=True)
t0 = time.time()

genes = pd.read_csv(os.path.join(SC, "count_matrix_genes.tsv"), header=None)[0].values
bcs   = pd.read_csv(os.path.join(SC, "count_matrix_barcodes.tsv"), header=None)[0].values
meta  = pd.read_csv(os.path.join(SC, "metadata.csv"), index_col=0)
print("genes", len(genes), "barcodes", len(bcs), "meta", meta.shape, flush=True)
assert list(meta.index) == list(bcs), "metadata row order != barcode order"

ct = meta["celltype_major"].astype(str).values
pt = meta["orig.ident"].astype(str).values
ncount = meta["nCount_RNA"].values.astype(np.float64)
ct_levels = sorted(set(ct)); pt_levels = sorted(set(pt))
ct_idx = pd.Categorical(ct, categories=ct_levels).codes.astype(np.int64)
pt_idx = pd.Categorical(pt, categories=pt_levels).codes.astype(np.int64)
nct, npt, ng = len(ct_levels), len(pt_levels), len(genes)
print("celltypes", nct, ct_levels, flush=True)
print("patients", npt, flush=True)

# per-cell group index
grp_ct = ct_idx                                  # 0..nct-1
grp_pc = pt_idx * nct + ct_idx                   # 0..npt*nct-1
npc = npt * nct
inv_ncount = 1.0 / np.maximum(ncount, 1.0)

sum_ct   = np.zeros(ng * nct, dtype=np.float64)   # raw UMI sums per celltype
cpm_ct   = np.zeros(ng * nct, dtype=np.float64)   # sum over cells of count/total*1e6
sum_pc   = np.zeros(ng * npc, dtype=np.float64)   # raw UMI sums per patient x celltype
tot_read = 0

hdr = 0
with open(os.path.join(SC, "count_matrix_sparse.mtx")) as fh:
    for line in fh:
        hdr += 1
        if not line.startswith("%"):
            dims = line.split(); break
print("mtx dims line:", dims, "header lines:", hdr, flush=True)
assert int(dims[0]) == ng and int(dims[1]) == len(bcs)
NNZ = int(dims[2])

reader = pd.read_csv(os.path.join(SC, "count_matrix_sparse.mtx"), sep=r"\s+", header=None,
                     names=["g", "c", "v"], skiprows=hdr, dtype=np.int32, chunksize=25_000_000,
                     engine="c")
for k, ch in enumerate(reader):
    g = ch["g"].values.astype(np.int64) - 1
    c = ch["c"].values.astype(np.int64) - 1
    v = ch["v"].values.astype(np.float64)
    tot_read += len(v)
    np.add(sum_ct, np.bincount(g * nct + grp_ct[c], weights=v, minlength=ng * nct), out=sum_ct)
    np.add(cpm_ct, np.bincount(g * nct + grp_ct[c], weights=v * inv_ncount[c] * 1e6,
                               minlength=ng * nct), out=cpm_ct)
    np.add(sum_pc, np.bincount(g * npc + grp_pc[c], weights=v, minlength=ng * npc), out=sum_pc)
    print(f"  chunk {k}: {len(v):,} nnz  cum {tot_read:,}/{NNZ:,}  {time.time()-t0:.0f}s", flush=True)

assert tot_read == NNZ, f"read {tot_read} != declared {NNZ}"
sum_ct = sum_ct.reshape(ng, nct)
cpm_ct = cpm_ct.reshape(ng, nct)
sum_pc = sum_pc.reshape(ng, npc)
n_cells_ct = np.bincount(ct_idx, minlength=nct).astype(np.float64)
n_cells_pc = np.bincount(grp_pc, minlength=npc).astype(np.float64)
mean_cpm_ct = cpm_ct / n_cells_ct[None, :]

np.savez_compressed(os.path.join(OUT, "wu2021_reference.npz"),
                    sum_ct=sum_ct, mean_cpm_ct=mean_cpm_ct, sum_pc=sum_pc,
                    n_cells_ct=n_cells_ct, n_cells_pc=n_cells_pc)
pd.DataFrame(sum_ct, index=genes, columns=ct_levels).to_csv(os.path.join(OUT, "wu2021_sum_umi_by_celltype.tsv"), sep="\t")
pd.DataFrame(mean_cpm_ct, index=genes, columns=ct_levels).to_csv(os.path.join(OUT, "wu2021_meancpm_by_celltype.tsv"), sep="\t")
pc_names = [f"{p}|{c}" for p in pt_levels for c in ct_levels]
keep = n_cells_pc >= 20
pd.DataFrame(sum_pc[:, keep], index=genes, columns=[pc_names[i] for i in np.where(keep)[0]]).to_csv(
    os.path.join(OUT, "wu2021_sum_umi_by_patient_celltype.tsv"), sep="\t")
pd.DataFrame({"group": pc_names, "n_cells": n_cells_pc}).to_csv(os.path.join(OUT, "wu2021_patient_celltype_ncells.tsv"), sep="\t", index=False)
pd.DataFrame({"celltype": ct_levels, "n_cells": n_cells_ct,
              "total_umi": sum_ct.sum(0)}).to_csv(os.path.join(OUT, "wu2021_celltype_ncells.tsv"), sep="\t", index=False)
print("cells per type:", dict(zip(ct_levels, n_cells_ct.astype(int))), flush=True)
print("kept patient x celltype groups with >=20 cells:", int(keep.sum()), "of", npc, flush=True)
print("DONE in", round(time.time() - t0), "s", flush=True)
