#!/usr/bin/env python3
"""Part A: inventory of the HPA files actually obtained, with release, URL, bytes, rows."""
import os, csv, hashlib
D="/path/to/revision/data/hpa"
R="/path/to/revision/results/v6"
files=[
 ("proteinatlas.tsv","proteinatlas.tsv.zip","HPA 25.1 (Ensembl 109)","https://www.proteinatlas.org/download/proteinatlas.tsv.zip",
  "consolidated per-gene table: prognostics, reliability, subcellular, antibodies"),
 ("pathology.tsv","pathology.tsv.zip","HPA 23.0","https://v23.proteinatlas.org/download/pathology.tsv.zip",
  "per-gene per-cancer IHC staining counts + prognostic p and direction"),
 ("normal_tissue.tsv","normal_tissue.tsv.zip","HPA 23.0","https://v23.proteinatlas.org/download/normal_tissue.tsv.zip",
  "IHC staining level per normal tissue and cell type"),
 ("rna_single_cell_type.tsv","rna_single_cell_type.tsv.zip","HPA 23.0","https://v23.proteinatlas.org/download/rna_single_cell_type.tsv.zip",
  "scRNA nTPM per pan-tissue cell type"),
 ("rna_single_cell_type_tissue.tsv","rna_single_cell_type_tissue.tsv.zip","HPA 23.0","https://v23.proteinatlas.org/download/rna_single_cell_type_tissue.tsv.zip",
  "scRNA nTPM per tissue-specific cluster (includes breast)"),
 ("rna_tissue_consensus.tsv","rna_tissue_consensus.tsv.zip","HPA 23.0","https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
  "consensus tissue RNA nTPM"),
 ("subcellular_location.tsv","subcellular_location.tsv.zip","HPA 23.0","https://v23.proteinatlas.org/download/subcellular_location.tsv.zip",
  "ICC/IF subcellular location and reliability"),
]
rows=[]
for tsv,zipf,rel,url,desc in files:
    p=os.path.join(D,tsv); z=os.path.join(D,zipf)
    n=sum(1 for _ in open(p,encoding="utf-8",errors="replace"))
    rows.append(dict(file=tsv, hpa_release=rel, url=url,
                     zip_bytes=os.path.getsize(z) if os.path.exists(z) else "",
                     tsv_bytes=os.path.getsize(p), rows_incl_header=n, data_rows=n-1,
                     contents=desc))
# per-gene XML
xd=os.path.join(D,"xml"); xs=sorted(os.listdir(xd))
rows.append(dict(file=f"xml/ (per-gene entries, {len(xs)} genes)", hpa_release="HPA 25.1 (Ensembl 109)",
                 url="https://www.proteinatlas.org/<ENSG>.xml", zip_bytes="",
                 tsv_bytes=sum(os.path.getsize(os.path.join(xd,x)) for x in xs),
                 rows_incl_header=len(xs), data_rows=len(xs),
                 contents="full current entry: normal IHC by cell type, per-patient cancer IHC, prognostics, subcellular, per-antibody validation"))
with open(os.path.join(R,"hpa_file_inventory.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows: print(r["file"], r["hpa_release"], r["data_rows"], "rows")
