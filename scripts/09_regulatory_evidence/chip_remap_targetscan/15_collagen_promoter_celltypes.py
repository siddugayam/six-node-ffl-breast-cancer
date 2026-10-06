#!/usr/bin/env python3
"""Which factors occupy the COL1A1 / COL3A1 promoters in FIBROBLAST, STROMAL or
BREAST biosamples specifically?  ReMap 2022 biotype strings, TSS +/- 1 kb and 5 kb."""
import json, os, collections, csv, re
BASE="/path/to/revision"
tss=json.load(open(BASE+"/cache/direct/node_tss_hg38.json"))
FIB=re.compile(r"fibro|stromal|stroma|mesenchym|MSC|myofibro|LX2|IMR-90|WI-38|BJ$|hMSC|SGBS|"
               r"myometri|leiomyoma|HFF|dermal|ESF|smooth-muscle|HSC",re.I)
BRE=re.compile(r"breast|mammary|MCF|MDA-MB|T-47D|SK-BR|SUM|HCC19|ZR-75|BT-4|MCF10|MCF-10|"
               r"HMEC|226LDM|WHIM",re.I)
rows=[]
for gene in ("COL1A1","COL3A1"):
    chrom,t,strand,tx=tss[gene]
    pk=[]
    for line in open(BASE+"/cache/direct/remap/regions/%s.bed"%gene):
        q=line.rstrip("\n").split("\t")
        if len(q)>=11: pk.append((int(q[6]),q[3],q[9],q[10]))
    for w in (1000,5000):
        m=collections.defaultdict(lambda: collections.defaultdict(set))
        for summit,ds,tf,bio in pk:
            if abs(summit-t)<=w: m[tf][bio].add(ds)
        for tf,bios in m.items():
            fb=[b for b in bios if FIB.search(b)]; br=[b for b in bios if BRE.search(b)]
            rows.append(dict(gene=gene,window_bp=w,TF=tf,
                n_datasets_total=sum(len(v) for v in bios.values()),
                n_biotypes_total=len(bios),
                n_datasets_fibroblast_stromal=sum(len(bios[b]) for b in fb),
                fibroblast_stromal_biotypes=";".join(sorted(fb)),
                n_datasets_breast=sum(len(bios[b]) for b in br),
                breast_biotypes=";".join(sorted(br))))
with open(BASE+"/results/v2/chip_collagen_celltype_context.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader()
    w.writerows(sorted(rows,key=lambda r:(r["gene"],r["window_bp"],
        -r["n_datasets_fibroblast_stromal"],-r["n_datasets_total"])))
print("WROTE chip_collagen_celltype_context.csv rows=",len(rows))
for gene in ("COL1A1","COL3A1"):
    for w in (1000,5000):
        sub=[r for r in rows if r["gene"]==gene and r["window_bp"]==w]
        fb=[r for r in sub if r["n_datasets_fibroblast_stromal"]>0]
        br=[r for r in sub if r["n_datasets_breast"]>0]
        print("\n== %s TSS+/-%dbp: %d TFs total; %d bound in a fibroblast/stromal biosample; "
              "%d in a breast biosample =="%(gene,w,len(sub),len(fb),len(br)))
        print("   fibroblast/stromal: "+", ".join("%s(%d:%s)"%(r["TF"],r["n_datasets_fibroblast_stromal"],
              r["fibroblast_stromal_biotypes"][:28]) for r in sorted(fb,key=lambda z:-z["n_datasets_fibroblast_stromal"])[:12]))
        print("   breast           : "+", ".join("%s(%d:%s)"%(r["TF"],r["n_datasets_breast"],
              r["breast_biotypes"][:28]) for r in sorted(br,key=lambda z:-z["n_datasets_breast"])[:12]))
