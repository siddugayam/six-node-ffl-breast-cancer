#!/usr/bin/env python3
"""31_gene_lengths.py -- union-exon length per gene symbol from GENCODE v36 basic,
used to test whether the log2(norm_count) -> TPM approximation in 22_deconv_prep.R
changes any fibroblast estimate."""
import gzip, os, collections
import pandas as pd
BASE = "/path/to/revision"
GTF  = os.path.join(BASE, "cache/v3/deconv/gc.gtf.gz")
ex = collections.defaultdict(list)
n = 0
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "exon": continue
        n += 1
        attr = f[8]
        i = attr.find('gene_name "')
        if i < 0: continue
        j = attr.find('"', i + 11)
        ex[attr[i+11:j]].append((int(f[3]), int(f[4])))
print("exon records:", n, "genes:", len(ex), flush=True)
out = {}
for g, iv in ex.items():
    iv.sort(); tot = 0; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    out[g] = tot + ce - cs + 1
pd.Series(out, name="union_exon_length").rename_axis("gene").to_csv(
    os.path.join(BASE, "cache/v3/deconv/gene_lengths_gencode_v36.tsv"), sep="\t")
s = pd.Series(out)
print("median length:", int(s.median()), "COL1A1:", out.get("COL1A1"), "ACTB:", out.get("ACTB"),
      "PTPRC:", out.get("PTPRC"), flush=True)
