#!/usr/bin/env python3
"""TargetScan 8.0 site-level evidence for the canonical network's miRNA_target edges.

Inputs (downloaded this run from targetscan.org/vert_80/vert_80_data_download/):
  miR_Family_Info.txt                                   family / seed / MiRBase ID
  Conserved_Site_Context_Scores.txt                     conserved sites, context++ per site
  Predicted_Targets_Context_Scores.default_predictions.txt  all default-prediction sites
  Summary_Counts.default_predictions.txt                per transcript x family aggregates
All restricted to Gene Tax ID / Species ID 9606.
Site Type numeric coding verified empirically against Predicted_Targets_Info
'Seed match' text: 1 = 7mer-A1, 2 = 7mer-m8, 3 = 8mer.
"""
import csv, os, re, sys, json, collections
BASE="/path/to/revision"
TS=BASE+"/cache/direct/ts"
OUT=BASE+"/results/v2"
SITETYPE={"1":"7mer-A1","2":"7mer-m8","3":"8mer","-2":"3'-compensatory","-3":"3'-compensatory"}

# ---------------- network ----------------
edges=list(csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t"))
mt=[e for e in edges if e["edge_type"]=="miRNA_target"]
print("miRNA_target edges:",len(mt))
net_mirs=sorted({e["source"] for e in mt})
net_tgts=sorted({e["target"] for e in mt})
print("distinct network miRNAs:",len(net_mirs),"distinct targets:",len(net_tgts))

ALIAS={"CTGF":"CCN2","CYR61":"CCN1","MKL1":"MRTFA","MLLT4":"AFDN","WISP1":"CCN4",
       "NOV":"CCN3","PTGS2":"PTGS2","KIAA1524":"CIP2A","C11orf30":"EMSY",
       "SEPP1":"SELENOP","ADAMTS5":"ADAMTS5","TMSB4X":"TMSB4X","H19":"H19"}

# ---------------- miRNA family info ----------------
fam_of={}; seed_of={}; cons_of={}
with open(TS+"/miR_Family_Info.txt") as f:
    h=next(f)
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[2]!="9606": continue
        fam_of[p[3]]=p[0]; seed_of[p[3]]=p[1]; cons_of[p[3]]=p[5]
print("human mature miRNAs in TargetScan:",len(fam_of))

def arms(label):
    """TargetScan MiRBase IDs belonging to a network miRNA label (all arms)."""
    pat=re.compile(r"^%s(-[12])?(-[35]p)?$"%re.escape(label))
    return sorted([m for m in fam_of if pat.match(m)])

mir_arms={m:arms(m) for m in net_mirs}
unmapped=[m for m,a in mir_arms.items() if not a]
print("network miRNAs mapped to >=1 TargetScan arm:",sum(1 for a in mir_arms.values() if a))
print("network miRNAs with NO TargetScan arm:",len(unmapped))
print("  unmapped:",unmapped)

# ---------------- site tables ----------------
def read_sites(path, tag):
    rows=[]
    with open(path) as f:
        next(f)
        for line in f:
            p=line.rstrip("\n").split("\t")
            if p[3]!="9606": continue
            rows.append((p[1],p[2],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12]))
    print("%s: human site rows = %d"%(tag,len(rows)))
    return rows
cons=read_sites(TS+"/Conserved_Site_Context_Scores.txt","Conserved_Site_Context_Scores")
allp=read_sites(TS+"/Predicted_Targets_Context_Scores.default_predictions.txt","Predicted_Targets_Context_Scores.default")

cons_key={(r[0],r[1],r[2],r[4],r[5]) for r in cons}   # gene,tx,miR,start,end
print("distinct conserved site keys:",len(cons_key))
n_in=sum(1 for r in allp if (r[0],r[1],r[2],r[4],r[5]) in cons_key)
print("default-prediction sites that are conserved:",n_in,"of",len(allp))

# index all sites by (gene symbol, miRNA)
site_idx=collections.defaultdict(list)
for r in allp:
    gene,tx,mir,st,us,ue,cs,csp,wcs,wcsp,kd=r
    site_idx[(gene,mir)].append(dict(transcript=tx,site_type=SITETYPE.get(st,st),
        site_type_code=st,utr_start=us,utr_end=ue,
        context_score=cs,context_pctile=csp,weighted_context=wcs,weighted_context_pctile=wcsp,
        rel_KD=kd,conserved=int((gene,tx,mir,us,ue) in cons_key)))
# conserved-only sites that are absent from default predictions (rare) get added
for r in cons:
    gene,tx,mir,st,us,ue,cs,csp,wcs,wcsp,kd=r
    k=(gene,mir)
    if not any(d["transcript"]==tx and d["utr_start"]==us and d["utr_end"]==ue for d in site_idx[k]):
        site_idx[k].append(dict(transcript=tx,site_type=SITETYPE.get(st,st),site_type_code=st,
            utr_start=us,utr_end=ue,context_score=cs,context_pctile=csp,weighted_context=wcs,
            weighted_context_pctile=wcsp,rel_KD=kd,conserved=1))
print("(gene,miRNA) pairs with >=1 site:",len(site_idx))

# ---------------- summary counts (per transcript x family) ----------------
summ=collections.defaultdict(list)
with open(TS+"/Summary_Counts.default_predictions.txt") as f:
    hdr=next(f).rstrip("\n").split("\t")
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[3]!="9606": continue
        summ[(p[1],p[2])].append(p)      # (gene symbol, miRNA family)
print("(gene,family) summary rows:",sum(len(v) for v in summ.values()))
COLS={c:i for i,c in enumerate(hdr)}

def fnum(x):
    try: return float(x)
    except: return None

def pair_summary(gene, mir):
    """Aggregate TargetScan evidence for one (gene symbol, mature miRNA)."""
    fam=fam_of.get(mir)
    sites=site_idx.get((gene,mir),[])
    best=None
    if fam and (gene,fam) in summ:
        for p in summ[(gene,fam)]:
            tot=fnum(p[COLS["Total context++ score"]])
            if tot is None: continue
            if best is None or tot < best[0]: best=(tot,p)
    d=dict(gene=gene, mirna=mir, family=fam, seed=seed_of.get(mir),
           family_conservation=cons_of.get(mir),
           n_sites=len(sites),
           n_conserved_sites=sum(s["conserved"] for s in sites),
           n_8mer=sum(1 for s in sites if s["site_type"]=="8mer"),
           n_7mer_m8=sum(1 for s in sites if s["site_type"]=="7mer-m8"),
           n_7mer_A1=sum(1 for s in sites if s["site_type"]=="7mer-A1"),
           n_8mer_conserved=sum(1 for s in sites if s["site_type"]=="8mer" and s["conserved"]),
           n_7mer_m8_conserved=sum(1 for s in sites if s["site_type"]=="7mer-m8" and s["conserved"]),
           n_7mer_A1_conserved=sum(1 for s in sites if s["site_type"]=="7mer-A1" and s["conserved"]))
    cscores=[fnum(s["context_score"]) for s in sites if fnum(s["context_score"]) is not None]
    d["best_site_context"]=min(cscores) if cscores else None
    pcts=[fnum(s["context_pctile"]) for s in sites if fnum(s["context_pctile"]) is not None]
    d["best_site_context_pctile"]=max(pcts) if pcts else None
    if best:
        p=best[1]
        d["transcript"]=p[COLS["Transcript ID"]]
        d["total_context_pp"]=fnum(p[COLS["Total context++ score"]])
        d["cum_weighted_context_pp"]=fnum(p[COLS["Cumulative weighted context++ score"]])
        d["aggregate_PCT"]=fnum(p[COLS["Aggregate PCT"]])
        d["summ_total_conserved_sites"]=p[COLS["Total num conserved sites"]]
        d["summ_total_nonconserved_sites"]=p[COLS["Total num nonconserved sites"]]
        d["summ_6mer_sites"]=p[COLS["Number of 6mer sites"]]
        d["representative_miRNA"]=p[COLS["Representative miRNA"]]
    else:
        for k in ("transcript","total_context_pp","cum_weighted_context_pp","aggregate_PCT",
                  "summ_total_conserved_sites","summ_total_nonconserved_sites",
                  "summ_6mer_sites","representative_miRNA"): d[k]=None
    return d, sites

# ---------------- named axes ----------------
NAMED=[("hsa-miR-29a","COL1A1"),("hsa-miR-29b","COL1A1"),("hsa-miR-29c","COL1A1"),
       ("hsa-miR-29a","COL3A1"),("hsa-miR-29b","COL3A1"),("hsa-miR-29c","COL3A1"),
       ("hsa-let-7b","COL3A1"),("hsa-let-7e","COL3A1"),
       ("hsa-miR-101","EZH2"),("hsa-miR-130a","VEGFA"),
       ("hsa-miR-124","STAT3"),("hsa-let-7b","HK2")]
net_edge_set={(e["source"],e["target"]) for e in mt}
named_rows=[]; named_sites=[]
for lab,gene in NAMED:
    g=ALIAS.get(gene,gene)
    for mir in (mir_arms.get(lab) or [lab]):
        d,sites=pair_summary(g,mir)
        d["network_label"]=lab; d["query_gene"]=gene
        d["in_network"]=int((lab,gene) in net_edge_set)
        named_rows.append(d)
        for s in sites:
            s2=dict(s); s2.update(network_label=lab,mirna=mir,gene=g)
            named_sites.append(s2)
print("named-axis (label,arm,gene) combinations:",len(named_rows))

# ---------------- network-wide ----------------
rows=[]
missing_gene=set(); no_arm=set()
ts_genes={g for (g,m) in site_idx}
for e in mt:
    lab=e["source"]; gene=ALIAS.get(e["target"],e["target"])
    arms_l=mir_arms.get(lab,[])
    if not arms_l: no_arm.add(lab)
    if gene not in ts_genes: missing_gene.add(e["target"])
    bestd=None
    for mir in arms_l:
        d,_=pair_summary(gene,mir)
        key=(d["total_context_pp"] if d["total_context_pp"] is not None else 0.0,
             -(d["n_sites"]))
        if bestd is None: bestd=(key,d,mir)
        elif key<bestd[0]: bestd=(key,d,mir)
    if bestd is None:
        rows.append(dict(source=lab,target=e["target"],ts_gene=gene,best_arm=None,
                         mapped_to_targetscan=0,n_sites=0,n_conserved_sites=0,
                         n_8mer=0,n_7mer_m8=0,n_7mer_A1=0,total_context_pp=None,
                         cum_weighted_context_pp=None,aggregate_PCT=None,
                         best_site_context=None,best_site_context_pctile=None,family=None))
    else:
        d=bestd[1]
        rows.append(dict(source=lab,target=e["target"],ts_gene=gene,best_arm=bestd[2],
            mapped_to_targetscan=1,n_sites=d["n_sites"],n_conserved_sites=d["n_conserved_sites"],
            n_8mer=d["n_8mer"],n_7mer_m8=d["n_7mer_m8"],n_7mer_A1=d["n_7mer_A1"],
            total_context_pp=d["total_context_pp"],cum_weighted_context_pp=d["cum_weighted_context_pp"],
            aggregate_PCT=d["aggregate_PCT"],best_site_context=d["best_site_context"],
            best_site_context_pctile=d["best_site_context_pctile"],family=d["family"]))
print("network edges scored:",len(rows))
print("edges whose miRNA has no TargetScan arm:",sum(1 for r in rows if not r["mapped_to_targetscan"]))
print("network miRNA labels with no TargetScan arm:",len(no_arm),sorted(no_arm))
print("network target symbols absent from TargetScan UTR set:",len(missing_gene),sorted(missing_gene)[:30])
print("edges with >=1 predicted site:",sum(1 for r in rows if r["n_sites"]>0))
print("edges with >=1 CONSERVED site:",sum(1 for r in rows if r["n_conserved_sites"]>0))

os.makedirs(OUT,exist_ok=True)
# per-site table for the named axes + per-edge network table, in one file
with open(OUT+"/targetscan_sites.csv","w",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["record_type","network_miRNA","targetscan_miRNA","gene","in_network","transcript",
                "site_type","utr_start","utr_end","conserved_site","context_pp_score",
                "context_pp_percentile","weighted_context_pp","weighted_context_pp_percentile",
                "predicted_relative_KD","miR_family","seed_m8","family_conservation",
                "n_sites","n_conserved_sites","n_8mer","n_7mer_m8","n_7mer_A1",
                "total_context_pp","cum_weighted_context_pp","aggregate_PCT",
                "best_site_context_pp","best_site_context_pp_percentile","best_arm"])
    for s in named_sites:
        w.writerow(["named_axis_site",s["network_label"],s["mirna"],s["gene"],"",s["transcript"],
                    s["site_type"],s["utr_start"],s["utr_end"],s["conserved"],s["context_score"],
                    s["context_pctile"],s["weighted_context"],s["weighted_context_pctile"],
                    s["rel_KD"],fam_of.get(s["mirna"]),seed_of.get(s["mirna"]),
                    cons_of.get(s["mirna"]),"","","","","","","","","","",""])
    for d in named_rows:
        w.writerow(["named_axis_pair",d["network_label"],d["mirna"],d["gene"],d["in_network"],
                    d["transcript"],"","","","", "","","","","",d["family"],d["seed"],
                    d["family_conservation"],d["n_sites"],d["n_conserved_sites"],d["n_8mer"],
                    d["n_7mer_m8"],d["n_7mer_A1"],d["total_context_pp"],
                    d["cum_weighted_context_pp"],d["aggregate_PCT"],d["best_site_context"],
                    d["best_site_context_pctile"],""])
    for r in rows:
        w.writerow(["network_edge",r["source"],r["best_arm"],r["ts_gene"],1,"","","","","",
                    "","","","","",r["family"],"","",r["n_sites"],r["n_conserved_sites"],
                    r["n_8mer"],r["n_7mer_m8"],r["n_7mer_A1"],r["total_context_pp"],
                    r["cum_weighted_context_pp"],r["aggregate_PCT"],r["best_site_context"],
                    r["best_site_context_pctile"],r["best_arm"]])
print("WROTE",OUT+"/targetscan_sites.csv")
json.dump({"n_edges":len(rows),
           "n_named_pairs":len(named_rows),
           "n_named_sites":len(named_sites)},
          open(BASE+"/cache/direct/ts_run_meta.json","w"))
