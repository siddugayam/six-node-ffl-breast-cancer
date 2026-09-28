#!/usr/bin/env python3
"""GRCh38 coordinates for all 30 prioritised nodes + comparators.
Protein-coding: Ensembl REST /lookup/symbol (GRCh38, current release) -- the same
reference the GWAS Catalog uses to map CHR_ID/CHR_POS.
miRNA: miRBase v23 GRCh38 GFF3."""
import re, sys, time, json
import requests, pandas as pd

DATA = "/path/to/revision/data"
RES  = "/path/to/revision/results/v6"

TFS   = ["E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3"]
GENES = ["CCND2","COL1A1","STAT5A","MYBL2","FN1","PDGFRB","MET","CXCL12","MMP14","PLAU"]
MIRS  = {"miR-21":["hsa-mir-21"], "miR-195":["hsa-mir-195"], "miR-204":["hsa-mir-204"],
         "miR-383":["hsa-mir-383"], "miR-124":["hsa-mir-124-1","hsa-mir-124-2","hsa-mir-124-3"],
         "miR-155":["hsa-mir-155"], "miR-429":["hsa-mir-429"], "miR-141":["hsa-mir-141"],
         "miR-34a":["hsa-mir-34a"], "miR-101":["hsa-mir-101-1","hsa-mir-101-2"]}
COMP  = ["COL3A1","NFKB1","RELA","SP1","ETS1","POSTN","COL1A2","DCN","LUM","EPCAM"]

S = requests.Session()
S.headers.update({"Content-Type":"application/json","Accept":"application/json"})
syms = TFS + GENES + COMP
CACHE = f"{DATA}/gwas/ensembl_symbol_lookup_cache.json"
try:
    res = json.load(open(CACHE))
except Exception:
    res = {}
for attempt in range(30):
    todo = [x for x in syms if x not in res]
    if not todo:
        break
    for sym in todo:
        try:
            r = S.get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}", timeout=60)
        except Exception as ex:
            print("retry", sym, type(ex).__name__, file=sys.stderr); time.sleep(3); continue
        if r.status_code == 200:
            try:
                j = r.json()
            except Exception:
                time.sleep(2); continue
            if isinstance(j, dict) and j.get("id"):
                res[sym] = j
                json.dump(res, open(CACHE, "w"))
        else:
            time.sleep(3)
        time.sleep(0.4)
missing = [x for x in syms if x not in res]
if missing:
    sys.exit(f"FATAL: Ensembl did not resolve {missing}")

rows = []
grp = {**{s:"TF" for s in TFS}, **{s:"gene" for s in GENES}, **{s:"comparator" for s in COMP}}
for s in syms:
    if s not in res: continue
    g = res[s]
    strand = "+" if g["strand"] == 1 else "-"
    rows.append(dict(node=s, node_type=grp[s], locus=s, chrom="chr"+str(g["seq_region_name"]),
                     start=int(g["start"]), end=int(g["end"]), strand=strand,
                     tss=int(g["start"]) if strand=="+" else int(g["end"]),
                     gencode_id=g["id"],
                     source=f"Ensembl REST lookup/symbol, {g.get('assembly_name')}"))

gff = {}
for ln in open(f"{DATA}/gwas/hsa_mirbase.gff3"):
    if ln.startswith("#"): continue
    f = ln.rstrip("\n").split("\t")
    if f[2] != "miRNA_primary_transcript": continue
    nm = re.search(r"Name=([^;]+)", f[8]).group(1)
    mi = re.search(r"ID=([^;]+)", f[8]).group(1)
    gff[nm] = (f[0], int(f[3]), int(f[4]), f[6], mi)
for node, loci in MIRS.items():
    for L in loci:
        if L not in gff:
            print(f"WARNING miRBase locus missing: {L}", file=sys.stderr); continue
        c,s,e,st,mi = gff[L]
        rows.append(dict(node=node, node_type="miRNA", locus=L, chrom=c, start=s, end=e,
                         strand=st, tss=(s if st=="+" else e), gencode_id=mi,
                         source="miRBase v23 GFF3, GRCh38"))

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/gwas_node_coordinates_grch38.csv", index=False)
print(df.to_string(index=False))
