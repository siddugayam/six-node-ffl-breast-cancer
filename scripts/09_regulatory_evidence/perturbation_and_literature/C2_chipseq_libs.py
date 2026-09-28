#!/usr/bin/env python3
"""C-ii/iii) Is COL1A1 / COL3A1 a target of NFKB1, RELA, SP1 or ETS1 according to
ChIP-seq consensus (ChEA/ENCODE), curated TRRUST, ARCHS4 co-expression, GEO TF
perturbation signatures, and LINCS L1000 CRISPR-KO consensus up/down sets?"""
import os, csv, re
LIBDIR = "/path/to/scratch/enrichr_libs"
OUT = "/path/to/revision/results/multiomics"
TFS = {"NFKB1","RELA","SP1","ETS1"}
READ = ["COL1A1","COL3A1"]
LIBS = ["ChEA_2022","ENCODE_TF_ChIP-seq_2015","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
        "TRRUST_Transcription_Factors_2019","ARCHS4_TFs_Coexp",
        "TF_Perturbations_Followed_by_Expression","LINCS_L1000_CRISPR_KO_Consensus_Sigs",
        "Enrichr_Submissions_TF-Gene_Coocurrence"]

rows = []
for lib in LIBS:
    p = os.path.join(LIBDIR, lib + ".txt")
    if not os.path.exists(p):
        print("MISSING LIBRARY", lib); continue
    nsets_tf = 0
    with open(p) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            genes = {g.split(",")[0].strip().upper() for g in parts[1:] if g.strip()}
            tok = re.split(r"[ _]", name)[0].upper()
            if tok not in TFS: continue
            nsets_tf += 1
            for g in READ:
                rows.append(dict(resource=lib, tf=tok, set_name=name, set_size=len(genes),
                                 readout_gene=g, present=g in genes))
    print(f"{lib}: {nsets_tf} sets for NFKB1/RELA/SP1/ETS1")

with open(os.path.join(OUT,"chipseq_TF_collagen_membership.csv"),"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["resource","tf","set_name","set_size","readout_gene","present"])
    w.writeheader(); [w.writerow(r) for r in rows]
print("wrote chipseq_TF_collagen_membership.csv rows:", len(rows))

# summary
import collections
agg = collections.defaultdict(lambda: [0,0])
for r in rows:
    k = (r["resource"], r["tf"], r["readout_gene"])
    agg[k][0] += 1
    agg[k][1] += 1 if r["present"] else 0
print(f"{'resource':<45}{'TF':<7}{'gene':<9}{'n_sets':>7}{'n_hit':>7}")
for k in sorted(agg):
    n, h = agg[k]
    print(f"{k[0]:<45}{k[1]:<7}{k[2]:<9}{n:>7}{h:>7}")
