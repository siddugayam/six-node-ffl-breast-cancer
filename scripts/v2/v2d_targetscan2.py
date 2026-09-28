#!/usr/bin/env python3
"""TargetScan 8.0 sequence-level evidence for the canonical network's miRNA_target
edges.  Uses ALL site classes (conserved + nonconserved) and the all_predictions
summary counts, both pre-filtered to Species/Gene Tax ID 9606 and to the 364 gene
symbols that appear as network targets.

Files downloaded this run from targetscan.org/vert_80/vert_80_data_download/:
  miR_Family_Info.txt
  Conserved_Site_Context_Scores.txt          (264,563 human site rows)
  Nonconserved_Site_Context_Scores.txt       (14,337,075 human site rows -> filtered)
  Summary_Counts.all_predictions.txt         (filtered to network genes)
Site Type numeric coding verified against Predicted_Targets_Info 'Seed match':
  1 = 7mer-A1, 2 = 7mer-m8, 3 = 8mer, negative = 3'-compensatory.
"""
import csv, json, os, collections, math
BASE="/path/to/revision"
TS=BASE+"/cache/direct/ts"; OUT=BASE+"/results/v2"
ST={"1":"7mer-A1","2":"7mer-m8","3":"8mer","-2":"3'-compensatory","-3":"3'-compensatory"}
Q=json.load(open(BASE+"/cache/direct/ts_query_sets.json"))
arm_map=Q["arm_map"]; GENES=set(Q["genes"])
ALIAS={"CTGF":"CCN2","CYR61":"CCN1","MKL1":"MRTFA"}
REV={v:k for k,v in ALIAS.items()}

fam_of={}; seed_of={}; famcons={}
for line in open(TS+"/miR_Family_Info.txt"):
    p=line.rstrip("\n").split("\t")
    if len(p)<6 or p[2]!="9606": continue
    fam_of[p[3]]=p[0]; seed_of[p[3]]=p[1]; famcons[p[3]]=p[5]

def fnum(x):
    try: return float(x)
    except: return None

sites=collections.defaultdict(list)   # (gene, mature miRNA) -> [site dicts]
def read(path, conserved, only_genes=True):
    n=0
    with open(path) as f:
        next(f)
        for line in f:
            p=line.rstrip("\n").split("\t")
            if p[3]!="9606": continue
            if only_genes and p[1] not in GENES: continue
            sites[(p[1],p[4])].append(dict(transcript=p[2],site_type=ST.get(p[5],p[5]),
                utr_start=p[6],utr_end=p[7],context=fnum(p[8]),context_pct=fnum(p[9]),
                wcontext=fnum(p[10]),wcontext_pct=fnum(p[11]),relKD=fnum(p[12]),
                conserved=conserved))
            n+=1
    print("%s -> %d site rows kept"%(os.path.basename(path),n))
    return n
nc=read(TS+"/Conserved_Site_Context_Scores.txt",1)
nn=read(TS+"/NC_sites_network.txt",0)
print("total (gene,miRNA) pairs with >=1 site:",len(sites))
print("total site rows:",nc+nn)

# NOTE: in Summary_Counts.*, the column headed "miRNA family" actually holds the
# 7-nt SEED (Seed+m8), not the family name -- verified by inspection.  Join on seed.
summ=collections.defaultdict(list)
with open(TS+"/SC_all_network.txt") as f:
    hdr=next(f).rstrip("\n").split("\t"); C={c:i for i,c in enumerate(hdr)}
    for line in f:
        p=line.rstrip("\n").split("\t")
        summ[(p[1],p[2])].append(p)
print("(gene,seed) summary rows:",sum(len(v) for v in summ.values()),
      "over",len(summ),"gene x seed pairs")
TS_GENES={g for (g,_) in sites} | {g for (g,_) in summ}
print("TargetScan human gene symbols represented among network targets:",len(TS_GENES))
def resolve(name):
    """TargetScan 8.0 uses 2018-vintage HGNC symbols (CTGF not CCN2 etc.);
    try the network symbol first, then the modern alias, then the reverse alias."""
    if name in TS_GENES: return name
    for alt in (ALIAS.get(name), REV.get(name)):
        if alt and alt in TS_GENES: return alt
    return name

def pair(gene,mir):
    s=sites.get((gene,mir),[]); fam=fam_of.get(mir); sd=seed_of.get(mir)
    best=None
    for p in summ.get((gene,sd),[]):
        t=fnum(p[C["Total context++ score"]])
        if t is None: continue
        if best is None or t<best[0]: best=(t,p)
    cs=[x["context"] for x in s if x["context"] is not None]
    pc=[x["context_pct"] for x in s if x["context_pct"] is not None]
    d=dict(gene=gene,mirna=mir,family=fam,seed=seed_of.get(mir),
        family_conservation=famcons.get(mir),n_sites=len(s),
        n_conserved_sites=sum(x["conserved"] for x in s),
        n_8mer=sum(1 for x in s if x["site_type"]=="8mer"),
        n_7mer_m8=sum(1 for x in s if x["site_type"]=="7mer-m8"),
        n_7mer_A1=sum(1 for x in s if x["site_type"]=="7mer-A1"),
        n_8mer_conserved=sum(1 for x in s if x["site_type"]=="8mer" and x["conserved"]),
        best_site_context=min(cs) if cs else None,
        best_site_context_pct=max(pc) if pc else None)
    if best:
        p=best[1]
        d.update(transcript=p[C["Transcript ID"]],
                 total_context=fnum(p[C["Total context++ score"]]),
                 cum_wcontext=fnum(p[C["Cumulative weighted context++ score"]]),
                 aggregate_PCT=fnum(p[C["Aggregate PCT"]]),
                 sc_conserved_sites=p[C["Total num conserved sites"]],
                 sc_nonconserved_sites=p[C["Total num nonconserved sites"]],
                 sc_6mer=p[C["Number of 6mer sites"]],
                 rep_miRNA=p[C["Representative miRNA"]])
    else:
        for k in ("transcript","total_context","cum_wcontext","aggregate_PCT",
                  "sc_conserved_sites","sc_nonconserved_sites","sc_6mer","rep_miRNA"): d[k]=None
    return d,s

edges=list(csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t"))
mt=[e for e in edges if e["edge_type"]=="miRNA_target"]
net={(e["source"],e["target"]) for e in mt}

NAMED=[("hsa-miR-29a","COL1A1"),("hsa-miR-29b","COL1A1"),("hsa-miR-29c","COL1A1"),
       ("hsa-miR-29a","COL3A1"),("hsa-miR-29b","COL3A1"),("hsa-miR-29c","COL3A1"),
       ("hsa-let-7b","COL3A1"),("hsa-let-7e","COL3A1"),("hsa-miR-101","EZH2"),
       ("hsa-miR-130a","VEGFA"),("hsa-miR-124","STAT3"),("hsa-let-7b","HK2")]
named_pairs=[]; named_sites=[]
print("\n=== NAMED AXES ===")
for lab,gene in NAMED:
    g=resolve(gene)
    for mir in arm_map[lab]:
        d,s=pair(g,mir); d["network_label"]=lab; d["query_gene"]=gene
        d["in_network"]=int((lab,gene) in net)
        named_pairs.append(d)
        for x in s:
            y=dict(x); y.update(network_label=lab,mirna=mir,gene=g); named_sites.append(y)
    tot=[p for p in named_pairs if p["network_label"]==lab and p["gene"]==g]
    print("%-14s -> %-7s  in_network=%d"%(lab,gene,int((lab,gene) in net)))
    for p in sorted(tot,key=lambda z:-(z["n_sites"])):
        if p["n_sites"]==0 and p["total_context"] is None: continue
        print("   %-20s fam=%-18s sites=%d (cons %d) 8mer=%d 7m8=%d 7A1=%d "
              "totalCtx++=%s aggPCT=%s bestSiteCtx++=%s (pct %s)"%(
              p["mirna"],p["family"],p["n_sites"],p["n_conserved_sites"],p["n_8mer"],
              p["n_7mer_m8"],p["n_7mer_A1"],p["total_context"],p["aggregate_PCT"],
              p["best_site_context"],p["best_site_context_pct"]))
    if not any(p["n_sites"] for p in tot):
        print("   NO TargetScan site of any class for any arm")

# ---- network-wide ----
rows=[]
for e in mt:
    lab=e["source"]; g=resolve(e["target"])
    best=None
    for mir in arm_map[lab]:
        d,_=pair(g,mir)
        k=(d["total_context"] if d["total_context"] is not None else 0.0,-d["n_sites"])
        if best is None or k<best[0]: best=(k,d,mir)
    d=best[1]
    rows.append(dict(source=lab,target=e["target"],ts_gene=g,best_arm=best[2],
        family=d["family"],n_sites=d["n_sites"],n_conserved_sites=d["n_conserved_sites"],
        n_8mer=d["n_8mer"],n_7mer_m8=d["n_7mer_m8"],n_7mer_A1=d["n_7mer_A1"],
        n_8mer_conserved=d["n_8mer_conserved"],
        total_context=d["total_context"],cum_wcontext=d["cum_wcontext"],
        aggregate_PCT=d["aggregate_PCT"],best_site_context=d["best_site_context"],
        best_site_context_pct=d["best_site_context_pct"],
        gene_in_targetscan=int(g in TS_GENES) ))
print("\n=== NETWORK-WIDE ===")
print("miRNA_target edges:",len(rows))
print("edges whose target gene is in the TargetScan human UTR set:",sum(r["gene_in_targetscan"] for r in rows))
print("genes NOT in TargetScan:",sorted({r["target"] for r in rows if not r["gene_in_targetscan"]}))
print("edges with >=1 site of any class:",sum(1 for r in rows if r["n_sites"]>0))
print("edges with >=1 CONSERVED site:",sum(1 for r in rows if r["n_conserved_sites"]>0))
print("edges with >=1 8mer:",sum(1 for r in rows if r["n_8mer"]>0))
print("edges with a total context++ score:",sum(1 for r in rows if r["total_context"] is not None))

with open(OUT+"/targetscan_sites.csv","w",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["record_type","network_miRNA","targetscan_miRNA","gene","in_network","transcript",
        "site_type","utr_start","utr_end","conserved_site","context_pp_score","context_pp_percentile",
        "weighted_context_pp","weighted_context_pp_percentile","predicted_relative_KD",
        "miR_family","seed_m8","family_conservation","n_sites","n_conserved_sites","n_8mer",
        "n_7mer_m8","n_7mer_A1","n_8mer_conserved","total_context_pp","cum_weighted_context_pp",
        "aggregate_PCT","best_site_context_pp","best_site_context_pp_percentile","best_arm",
        "gene_in_targetscan"])
    for s in named_sites:
        w.writerow(["named_axis_site",s["network_label"],s["mirna"],s["gene"],"",s["transcript"],
            s["site_type"],s["utr_start"],s["utr_end"],s["conserved"],s["context"],s["context_pct"],
            s["wcontext"],s["wcontext_pct"],s["relKD"],fam_of.get(s["mirna"]),seed_of.get(s["mirna"]),
            famcons.get(s["mirna"]),"","","","","","","","","","","","",""])
    for d in named_pairs:
        w.writerow(["named_axis_pair",d["network_label"],d["mirna"],d["gene"],d["in_network"],
            d["transcript"],"","","","","","","","","",d["family"],d["seed"],d["family_conservation"],
            d["n_sites"],d["n_conserved_sites"],d["n_8mer"],d["n_7mer_m8"],d["n_7mer_A1"],
            d["n_8mer_conserved"],d["total_context"],d["cum_wcontext"],d["aggregate_PCT"],
            d["best_site_context"],d["best_site_context_pct"],"",""])
    for r in rows:
        w.writerow(["network_edge",r["source"],r["best_arm"],r["ts_gene"],1,"","","","","","","","","","",
            r["family"],"","",r["n_sites"],r["n_conserved_sites"],r["n_8mer"],r["n_7mer_m8"],
            r["n_7mer_A1"],r["n_8mer_conserved"],r["total_context"],r["cum_wcontext"],
            r["aggregate_PCT"],r["best_site_context"],r["best_site_context_pct"],r["best_arm"],
            r["gene_in_targetscan"]])
print("WROTE",OUT+"/targetscan_sites.csv")
json.dump(rows,open(BASE+"/cache/direct/ts_network_rows.json","w"))
