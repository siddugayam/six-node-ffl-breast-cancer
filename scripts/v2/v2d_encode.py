#!/usr/bin/env python3
"""ENCODE portal: number of released human TF ChIP-seq experiments per network TF,
read from the target.label facet of a single faceted search."""
import csv, json, urllib.request, collections
BASE="/path/to/revision"
U=("https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq"
   "&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"
   "&status=released&limit=0")
req=urllib.request.Request(U,headers={"Accept":"application/json","User-Agent":"revision-check"})
d=json.load(urllib.request.urlopen(req,timeout=180))
print("released human TF ChIP-seq experiments on ENCODE:",d["total"])
fac={f["field"]:f for f in d["facets"]}
terms={t["key"]:t["doc_count"] for t in fac["target.label"]["terms"]}
print("distinct targets:",len(terms))
bios={t["key"]:t["doc_count"] for t in fac["biosample_ontology.term_name"]["terms"]}
edges=list(csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t"))
tfs=sorted({e["source"] for e in edges if e["edge_type"]=="TF_target"})
rows=[]
for tf in tfs:
    rows.append(dict(TF=tf,encode_human_TF_ChIPseq_experiments=terms.get(tf,0)))
n0=sum(1 for r in rows if r["encode_human_TF_ChIPseq_experiments"]==0)
print("network TFs with >=1 ENCODE experiment: %d/%d ; with none: %d"%(len(tfs)-n0,len(tfs),n0))
for tf in ("NFKB1","RELA","SP1","ETS1","SMAD3","TWIST1","MYB","TFAP2A","STAT6","MKL1"):
    print("   %-8s %d"%(tf,terms.get(tf,0)))
with open(BASE+"/results/v2/chip_encode_counts_v2d.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=["TF","encode_human_TF_ChIPseq_experiments"])
    w.writeheader(); w.writerows(rows)
print("WROTE chip_encode_counts_v2d.csv")
# MCF-7 / breast context
U2=U+"&biosample_ontology.term_name=MCF-7"
d2=json.load(urllib.request.urlopen(urllib.request.Request(U2,headers={"Accept":"application/json"}),timeout=180))
t2={t["key"]:t["doc_count"] for t in {f["field"]:f for f in d2["facets"]}["target.label"]["terms"]}
print("ENCODE MCF-7 TF ChIP-seq experiments:",d2["total"],"; distinct targets:",len(t2))
print("  of the 4 named TFs in MCF-7:",{k:t2.get(k,0) for k in ("NFKB1","RELA","SP1","ETS1")})
json.dump({"human_total":d["total"],"mcf7_total":d2["total"],
           "mcf7_named":{k:t2.get(k,0) for k in ("NFKB1","RELA","SP1","ETS1")}},
          open(BASE+"/cache/direct/encode_meta.json","w"))
