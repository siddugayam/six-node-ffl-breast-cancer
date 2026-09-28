#!/usr/bin/env python3
"""Cellular source of the miRNA loci, approached through their pri-miRNA HOST GENES.
Droplet scRNA-seq cannot measure mature miRNAs; the host transcript is the only
cell-type-resolved read-out available for a miRNA locus in GSE176078."""
import csv, collections, re, math
CACHE="/path/to/revision/cache/celltype/"
RES="/path/to/revision/results/multiomics/"
codes={}
for l in open(CACHE+"celltype_codes.tsv"):
    a,b=l.rstrip("\n").split("\t"); codes[a]=b
cell_ct={int(r["cell_index"]):r["celltype_major"]
         for r in csv.DictReader(open(CACHE+"cell_index_celltype.tsv"),delimiter="\t")}
tot=collections.Counter(); ncell=collections.Counter(cell_ct.values())
for l in open(CACHE+"scrna_cell_totals.tsv"):
    a,b=l.split(); tot[cell_ct[int(a)]]+=int(b)
cpm=collections.defaultdict(dict); nexp=collections.defaultdict(dict)
for g,cc,s,n in csv.reader(open(CACHE+"scrna_pseudobulk_all.tsv"),delimiter="\t"):
    if g=="gene": continue
    ct=codes[cc]; cpm[g][ct]=int(s)/tot[ct]*1e6; nexp[g][ct]=int(n)
net_mir=[r[0] for r in csv.reader(open("data/canonical_nodes.tsv"),delimiter="\t") if r[1]=="miRNA"]
netset=set(m.lower() for m in net_mir)

def host_to_mirna(h):
    core=h[:-2]
    if core.startswith("MIRLET"):
        m="hsa-let-"+core[6:].lower()          # MIRLET7B -> hsa-let-7b
    else:
        m="hsa-miR-"+core[3:]
    m=re.sub(r'-(\d)$', r'-\1', m)
    return m
cts=sorted(tot); rows=[]
for h in sorted(g for g in cpm if re.match(r'^MIR(LET)?[0-9A-Za-z\-]*HG$', g)):
    cand=host_to_mirna(h)
    variants={cand, re.sub(r'\d+$','',cand), re.sub(r'-\d$','',cand)}
    matched=sorted({m for m in net_mir if m.lower() in {v.lower() for v in variants}})
    d=cpm[h]; top=max(cts,key=lambda c: d.get(c,0))
    caf=d.get("CAFs",0.0); can=d.get("Cancer Epithelial",0.0)
    rows.append(dict(host_gene=h, implied_miRNA=cand,
        in_network=("YES" if matched else "no"), matched_network_miRNA=";".join(matched),
        top_celltype=top, top_CPM=round(d.get(top,0),3),
        CAF_CPM=round(caf,3), CancerEpi_CPM=round(can,3),
        log2_CAF_over_CancerEpi=round(math.log2((caf+0.1)/(can+0.1)),3),
        pct_CAFs_expressing=round(100*nexp[h].get("CAFs",0)/ncell["CAFs"],2),
        pct_CancerEpi_expressing=round(100*nexp[h].get("Cancer Epithelial",0)/ncell["Cancer Epithelial"],2),
        **{("CPM_"+c.replace(" ","")): round(d.get(c,0),3) for c in cts}))
with open(RES+"celltype_source_mirna_hostgenes.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("host genes in atlas:",len(rows), "| matching a network miRNA:", sum(r["in_network"]=="YES" for r in rows))
for r in rows:
    if r["in_network"]=="YES":
        print("%-12s -> %-14s top=%-18s CAF=%8.2f Epi=%8.2f log2FC=%6.2f"%(
            r["host_gene"],r["matched_network_miRNA"],r["top_celltype"],r["CAF_CPM"],r["CancerEpi_CPM"],r["log2_CAF_over_CancerEpi"]))
print("\nMIR29 host genes present in atlas:", [g for g in cpm if "MIR29" in g] or "NONE")
