#!/usr/bin/env python3
"""Consolidated verification table for the direct-molecular-evidence task."""
import csv, json, os
BASE="/path/to/revision"; OUT=BASE+"/results/v2"
rows=[]
def add(sec,item,value,source):
    rows.append(dict(section=sec,item=item,value=value,source=source))
inv={r["TF"]:r for r in csv.DictReader(open(OUT+"/chip_atlas_tf_inventory_v2d.csv"))}
enc={r["TF"]:r["encode_human_TF_ChIPseq_experiments"] for r in csv.DictReader(open(OUT+"/chip_encode_counts_v2d.csv"))}
spec={}
for r in csv.DictReader(open(OUT+"/chip_collagen_promoter_specificity.csv")):
    spec[(r["gene"],r["window_bp"],r["TF"])]=r
ca={}
for r in csv.DictReader(open(OUT+"/chip_evidence.csv")): ca[(r["TF"],r["target"],r["window_kb"])]=r
for tf in ("NFKB1","RELA","SP1","ETS1"):
    i=inv.get(tf,{})
    add("A_named_TF",tf+"_chipatlas_hg38_experiments",i.get("n_experiments_table",""),"ChIP-Atlas target table header")
    add("A_named_TF",tf+"_encode_human_experiments",enc.get(tf,"0"),"ENCODE portal facet")
    for g in ("COL1A1","COL3A1"):
        a=ca.get((tf,g,"1"))
        add("A_named_TF","%s_%s_chipatlas_bound_1kb"%(tf,g),
            ("%s/%s experiments"%(a["n_experiments_bound"],a["n_experiments"]) if a else "edge not in network"),
            "ChIP-Atlas +/-1kb")
        s=spec.get((g,"1000",tf))
        add("A_named_TF","%s_%s_remap_datasets_1kb"%(tf,g),
            (s["n_datasets_at_promoter"] if s else "NA"),"ReMap 2022 +/-1kb")
        add("A_named_TF","%s_%s_remap_enrichment_vs_own_background"%(tf,g),
            (s["enrichment_vs_background"] if s else "NA"),"ReMap vs 600 random promoters")
summ={ (r["window_kb"]):r for r in csv.DictReader(open(OUT+"/chip_edge_support_summary_v2d.csv"))}
for w in ("1","10"):
    s=summ[w]
    add("A_edge_support","chipatlas_%skb_testable_edges"%w,s["n_edges_testable"],"ChIP-Atlas")
    add("A_edge_support","chipatlas_%skb_supported_ge2exp"%w,"%s (%s)"%(s["n_supported_ge2_exp"],s["frac_supported_ge2"]),"ChIP-Atlas")
    add("A_edge_support","chipatlas_%skb_null_mean_ge2"%w,"%s +/- %s (%s reps)"%(s["null_mean_ge2"],s["null_sd_ge2"],s["null_reps"]),"random gene from ChIP-Atlas universe")
    add("A_edge_support","chipatlas_%skb_fold"%w,s["fold_enrichment_ge2"],"observed/null")
for r in csv.DictReader(open(OUT+"/chip_edge_null_genomewide.csv")):
    add("A_edge_support","remap_%sbp_ge%s_dataset"%(r["window_bp"],r["min_datasets"]),
        "obs %s/%s (%s) vs null %s (%s), Z=%s fold=%s p=%s, %s reps"%(
        r["observed"],r["n_edges"],r["obs_frac"],r["null_mean"],r["null_frac"],r["Z"],r["fold"],r["p_emp"],r["reps"]),
        "ReMap vs 600 random non-network promoters")
for r in csv.DictReader(open(OUT+"/targetscan_site_null.csv")):
    add("B_site_null","targetscan_%s"%r["metric"],
        "obs %s (%s) vs null %s (%s), Z=%s fold=%s p=%s, %s reps"%(
        r["observed"],r["obs_frac"],r["null_mean"],r["null_frac"],r["Z"],r["fold"],r["p_emp"],r["reps"]),
        "TargetScan 8.0, target reshuffled within network genes")
for r in csv.DictReader(open(OUT+"/targetscan_concordance_stats_v2d.csv")):
    add("B_site_quality_vs_concordance",r["metric"],r["value"],"TargetScan 8.0 x TCGA-BRCA n=1066")
ts={}
for r in csv.DictReader(open(OUT+"/targetscan_sites.csv")):
    if r["record_type"]=="named_axis_pair":
        ts.setdefault((r["network_miRNA"],r["gene"]),[]).append(r)
for (m,g),v in sorted(ts.items()):
    best=sorted(v,key=lambda x:float(x["total_context_pp"] or 0))[0]
    add("B_named_axis","%s -> %s"%(m,g),
        "arm=%s sites=%s (conserved %s) 8mer=%s 7mer-m8=%s 7mer-A1=%s totalCtx++=%s aggPCT=%s in_network=%s"%(
        best["targetscan_miRNA"],best["n_sites"],best["n_conserved_sites"],best["n_8mer"],
        best["n_7mer_m8"],best["n_7mer_A1"],best["total_context_pp"],best["aggregate_PCT"],best["in_network"]),
        "TargetScan 8.0")
with open(OUT+"/verify_direct_evidence.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=["section","item","value","source"]); w.writeheader(); w.writerows(rows)
print("WROTE verify_direct_evidence.csv rows=",len(rows))
for r in rows:
    if r["section"] in ("A_named_TF",): print("  %-55s %s"%(r["item"],r["value"]))
