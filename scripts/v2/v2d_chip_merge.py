#!/usr/bin/env python3
"""Merge the three ChIP resources into one per-edge evidence table for Table 1,
and produce the headline coverage fractions."""
import csv, json, collections
BASE="/path/to/revision"; OUT=BASE+"/results/v2"
edges=[e for e in csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t")
       if e["edge_type"]=="TF_target"]
ca={ }
for r in csv.DictReader(open(OUT+"/chip_evidence.csv")):
    ca[(r["TF"],r["target"],r["window_kb"])]=r
rm={(r["TF"],r["target"]):r for r in csv.DictReader(open(OUT+"/remap_edge_support.csv"))}
enc={r["TF"]:int(r["encode_human_TF_ChIPseq_experiments"])
     for r in csv.DictReader(open(OUT+"/chip_encode_counts_v2d.csv"))}
inv={r["TF"]:r for r in csv.DictReader(open(OUT+"/chip_atlas_tf_inventory_v2d.csv"))}
REMAP_CAT=set(json.load(open(BASE+"/cache/direct/remap_tf_catalogue.json")))
rows=[]
for e in edges:
    tf,tg=e["source"],e["target"]
    a1=ca.get((tf,tg,"1")); a10=ca.get((tf,tg,"10")); r=rm.get((tf,tg))
    d=dict(TF=tf,target=tg,
        chipatlas_table_available=int(a1 is not None),
        chipatlas_n_experiments=(a1["n_experiments"] if a1 else ""),
        chipatlas_bound_1kb=(a1["n_experiments_bound"] if a1 else ""),
        chipatlas_frac_1kb=(a1["frac_bound"] if a1 else ""),
        chipatlas_meanMACS2_1kb=(a1["mean_MACS2_all_exp"] if a1 else ""),
        chipatlas_bound_10kb=(a10["n_experiments_bound"] if a10 else ""),
        chipatlas_support_ge2exp_1kb=(a1["chip_support_ge2"] if a1 else ""),
        chipatlas_cells_1kb=(a1["cell_types_bound"] if a1 else ""),
        remap_in_catalogue=int(tf in REMAP_CAT),
        remap_datasets_1kb=(r["remap_datasets_1000"] if r else ""),
        remap_datasets_5kb=(r["remap_datasets_5000"] if r else ""),
        remap_support_ge2_1kb=((r["remap_support_ge2_1000"] if r else "") if tf in REMAP_CAT else ""),
        remap_biotypes_1kb=(r["remap_biotypes_1000"] if r else ""),
        encode_n_experiments=enc.get(tf,0))
    d["any_chip_support_ge2"]=int((d["chipatlas_support_ge2exp_1kb"]=="1") or
                                   (d["remap_support_ge2_1kb"]=="1"))
    d["both_resources_support_ge2"]=int((d["chipatlas_support_ge2exp_1kb"]=="1") and
                                        (d["remap_support_ge2_1kb"]=="1"))
    rows.append(d)
with open(OUT+"/chip_edge_evidence_merged.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
n=len(rows)
ca_t=[r for r in rows if r["chipatlas_table_available"]==1]
rm_t=[r for r in rows if r["remap_in_catalogue"]==1]
both=[r for r in rows if r["chipatlas_table_available"]==1 and r["remap_in_catalogue"]==1]
def f(sub,k): return sum(1 for r in sub if r[k]=="1" or r[k]==1)
print("TF_target edges:",n)
print("ChIP-Atlas testable:",len(ca_t)," supported (>=2 exp, +/-1kb):",f(ca_t,"chipatlas_support_ge2exp_1kb"),
      "= %.1f%%"%(100*f(ca_t,"chipatlas_support_ge2exp_1kb")/len(ca_t)))
print("ReMap testable    :",len(rm_t)," supported (>=2 datasets, +/-1kb):",f(rm_t,"remap_support_ge2_1kb"),
      "= %.1f%%"%(100*f(rm_t,"remap_support_ge2_1kb")/len(rm_t)))
print("edges testable in BOTH:",len(both))
print("  supported by BOTH:",f(both,"both_resources_support_ge2"),
      "= %.1f%% of the 770 published edges"%(100*f(both,"both_resources_support_ge2")/n))
print("  supported by EITHER:",f(both,"any_chip_support_ge2"),
      "= %.1f%% of testable"%(100*f(both,"any_chip_support_ge2")/len(both)))
print("  supported by EITHER, as a fraction of all 770 published edges: %.1f%%"%(
      100*sum(1 for r in rows if r["any_chip_support_ge2"]==1)/n))
# agreement between resources
agree=sum(1 for r in both if (r["chipatlas_support_ge2exp_1kb"]=="1")==(r["remap_support_ge2_1kb"]=="1"))
print("  resource agreement on the binary call: %d/%d = %.1f%%"%(agree,len(both),100*agree/len(both)))
notest=[r for r in rows if r["chipatlas_table_available"]==0 and r["remap_in_catalogue"]==0]
print("edges with NO ChIP data in either resource:",len(notest),
      "(TFs: %s)"%sorted({r["TF"] for r in notest}))
print("network TFs never assayed by ENCODE:",sum(1 for t in {r['TF'] for r in rows} if enc.get(t,0)==0),
      "of",len({r['TF'] for r in rows}))
