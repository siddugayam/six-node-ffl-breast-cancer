#!/usr/bin/env python3
"""Derive a data-driven CAF signature from GSE176078 that is deliberately
COLLAGEN-FREE and NETWORK-NODE-FREE, so scoring TCGA bulk with it and then
partialling it out of the collagen correlations cannot be circular."""
import csv, collections, math
CACHE = "/path/to/revision/cache/celltype/"
codes = {}
for line in open(CACHE + "celltype_codes.tsv"):
    a, b = line.rstrip("\n").split("\t"); codes[a] = b
# total UMI per cell type (all genes)
tot = collections.Counter()
cell_ct = {}
with open(CACHE + "cell_index_celltype.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cell_ct[int(r["cell_index"])] = r["celltype_major"]
for line in open(CACHE + "scrna_cell_totals.tsv"):
    a, b = line.split(); tot[cell_ct[int(a)]] += int(b)

cpm = collections.defaultdict(dict)
with open(CACHE + "scrna_pseudobulk_all.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); next(rd)
    for g, cc, s, n in rd:
        ct = codes[cc]
        cpm[g][ct] = int(s) / tot[ct] * 1e6
print("genes in pseudobulk:", len(cpm))
network = set(l.strip() for l in open(CACHE + "network_protein_nodes.txt"))
cts = sorted(tot)
cand = []
for g, d in cpm.items():
    caf = d.get("CAFs", 0.0)
    others = [d.get(c, 0.0) for c in cts if c != "CAFs"]
    nxt = max(others)
    if caf < 50: continue
    spec = math.log2((caf + 1) / (nxt + 1))
    if spec < 2: continue
    if g.startswith("COL"): continue
    if g in network: continue
    cand.append((spec, caf, g))
cand.sort(reverse=True)
sig = [g for _, _, g in cand[:50]]
open(CACHE + "caf_signature_scrna.txt", "w").write("\n".join(sig) + "\n")
print("CAF-specific candidates passing filters:", len(cand))
print("signature (top 50):", ", ".join(sig))
excl_col = [g for g in cpm if g.startswith("COL") and cpm[g].get("CAFs",0) >= 50
            and math.log2((cpm[g]["CAFs"]+1)/(max(cpm[g].get(c,0) for c in cts if c!="CAFs")+1)) >= 2]
print("collagen genes that WOULD have qualified and were excluded:", ", ".join(sorted(excl_col)))
