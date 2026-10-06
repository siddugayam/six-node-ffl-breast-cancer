#!/usr/bin/env python3
"""ChIP-seq support for the canonical network's 770 TF_target edges, from
ChIP-Atlas hg38 'Target Genes' tables (average and per-experiment MACS2
-10*log10(p) scores in TSS +/- 1 kb and +/- 10 kb windows)."""
import csv, os, sys, json, collections, random, math
BASE="/path/to/revision"
TDIR=BASE+"/cache/direct/chipatlas/target"
OUT=BASE+"/results/v2"

edges=list(csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t"))
tfe=[e for e in edges if e["edge_type"]=="TF_target"]
print("TF_target edges:",len(tfe))
tf_targets=collections.defaultdict(set)
for e in tfe: tf_targets[e["source"]].add(e["target"])
tfs=sorted(tf_targets); print("distinct TFs:",len(tfs))
ALIAS={"CTGF":"CCN2","CYR61":"CCN1","MKL1":"MRTFA"}

# ---- experiment counts from ChIP-Atlas metadata (hg38, "TFs and others") ----
meta_n=collections.Counter(); meta_cells=collections.defaultdict(set)
with open(BASE+"/cache/direct/chipatlas/experimentList.tab",encoding="utf-8",errors="replace") as f:
    for line in f:
        p=line.rstrip("\n").split("\t")
        if len(p)<6: continue
        if p[1]=="hg38" and p[2]=="TFs and others":
            meta_n[p[3]]+=1; meta_cells[p[3]].add(p[5])
print("hg38 TF antigens in ChIP-Atlas metadata:",len(meta_n))

def load(tf,dist):
    path=os.path.join(TDIR,"%s.%s.tsv"%(tf,dist))
    if not os.path.exists(path) or os.path.getsize(path)<100: return None
    with open(path,encoding="utf-8",errors="replace") as f:
        hdr=f.readline().rstrip("\n").split("\t")
        exps=[]
        for c in hdr[2:]:
            srx,_,cell=c.partition("|")
            exps.append((srx,cell))
        genes={}
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)<2: continue
            vals=[]
            for v in p[2:]:
                try: vals.append(float(v))
                except: vals.append(0.0)
            genes[p[0]]=vals
    return dict(exps=exps,genes=genes)

rows=[]; inv=[]; reverse=[]
uni1=set(); uni10=set()
per_tf_bound={}
BOUND={}   # (tf,dist) -> (set bound in >=1 exp, set bound in >=2 exp)
for tf in tfs:
    for dist in ("1","10"):
        d=load(tf,dist)
        if d is None:
            continue
        (uni1 if dist=="1" else uni10).update(d["genes"])
        b1={g for g,v in d["genes"].items() if any(x>0 for x in v)}
        b2={g for g,v in d["genes"].items() if sum(1 for x in v if x>0)>=2}
        BOUND[(tf,dist)]=(b1,b2)
        nb1=len(b1); nb2=len(b2)
        per_tf_bound[(tf,dist)]=(nb1,nb2,len(d["genes"]),len(d["exps"]))
        if dist=="1":
            inv.append(dict(TF=tf,n_experiments_table=len(d["exps"]),
                n_experiments_metadata=meta_n.get(tf,0),
                n_genes_listed=len(d["genes"]),
                n_genes_bound_ge1_exp=nb1,n_genes_bound_ge2_exp=nb2,
                n_cell_types=len({c for _,c in d["exps"]}),
                n_network_targets=len(tf_targets[tf])))
        # reverse query for the collagens
        for col in ("COL1A1","COL3A1"):
            if col in d["genes"]:
                v=d["genes"][col]
                nb=sum(1 for x in v if x>0)
                cells=sorted({d["exps"][i][1] for i,x in enumerate(v) if x>0})
                reverse.append(dict(gene=col,TF=tf,window_kb=dist,
                    n_experiments=len(d["exps"]),n_experiments_bound=nb,
                    frac_bound=round(nb/len(d["exps"]),4) if d["exps"] else None,
                    mean_MACS2_all_exp=round(sum(v)/len(v),3) if v else None,
                    max_MACS2=max(v) if v else None,
                    cell_types_bound=";".join(cells[:25]),
                    n_cell_types_bound=len(cells),
                    is_network_TF_of_gene=int(col in tf_targets[tf])))
        # edge-level
        for tgt in sorted(tf_targets[tf]):
            g=tgt if tgt in d["genes"] else ALIAS.get(tgt,tgt)
            v=d["genes"].get(g)
            if v is None:
                rows.append(dict(TF=tf,target=tgt,window_kb=dist,in_table=0,
                    n_experiments=len(d["exps"]),n_experiments_bound=0,frac_bound=0.0,
                    mean_MACS2_all_exp=0.0,max_MACS2=0.0,n_cell_types_bound=0,
                    cell_types_bound="",chip_support_ge1=0,chip_support_ge2=0))
            else:
                nb=sum(1 for x in v if x>0)
                cells=sorted({d["exps"][i][1] for i,x in enumerate(v) if x>0})
                rows.append(dict(TF=tf,target=tgt,window_kb=dist,in_table=1,
                    n_experiments=len(d["exps"]),n_experiments_bound=nb,
                    frac_bound=round(nb/len(d["exps"]),4) if d["exps"] else None,
                    mean_MACS2_all_exp=round(sum(v)/len(v),3),max_MACS2=max(v),
                    n_cell_types_bound=len(cells),cell_types_bound=";".join(cells[:25]),
                    chip_support_ge1=int(nb>=1),chip_support_ge2=int(nb>=2)))
print("ChIP-Atlas gene universe: 1kb=%d  10kb=%d"%(len(uni1),len(uni10)))
tf_with_1=len({t for (t,d) in per_tf_bound if d=="1"})
tf_with_10=len({t for (t,d) in per_tf_bound if d=="10"})
print("TFs with a +/-1kb table: %d/%d ; +/-10kb table: %d/%d"%(tf_with_1,len(tfs),tf_with_10,len(tfs)))
missing=sorted(set(tfs)-{t for (t,d) in per_tf_bound if d=="1"})
print("TFs with NO +/-1kb table:",missing)

# ---------------- edge support summary + matched null ----------------
random.seed(20260909)
summary=[]
for dist in ("1","10"):
    sub=[r for r in rows if r["window_kb"]==dist]
    have=[r for r in sub if (r["TF"],dist) in per_tf_bound]
    n=len(have)
    s1=sum(r["chip_support_ge1"] for r in have); s2=sum(r["chip_support_ge2"] for r in have)
    uni=len(uni1 if dist=="1" else uni10)
    exp1=sum(per_tf_bound[(r["TF"],dist)][0]/uni for r in have)
    exp2=sum(per_tf_bound[(r["TF"],dist)][1]/uni for r in have)
    # empirical null: for each edge draw a random gene from the universe, 200 reps
    U=sorted(uni1 if dist=="1" else uni10)
    Uset={}
    reps=200
    null1=[];null2=[]
    bysets={tf:BOUND[(tf,dist)] for tf in {r["TF"] for r in have}}
    for _ in range(reps):
        a=b=0
        for r in have:
            g=U[random.randrange(len(U))]
            s=bysets[r["TF"]]
            if g in s[0]: a+=1
            if g in s[1]: b+=1
        null1.append(a); null2.append(b)
    m1=sum(null1)/len(null1); m2=sum(null2)/len(null2)
    sd1=math.sqrt(sum((x-m1)**2 for x in null1)/(len(null1)-1))
    sd2=math.sqrt(sum((x-m2)**2 for x in null2)/(len(null2)-1))
    p1=(sum(1 for x in null1 if x>=s1)+1)/(reps+1)
    p2=(sum(1 for x in null2 if x>=s2)+1)/(reps+1)
    summary.append(dict(window_kb=dist,n_edges_testable=n,n_edges_total=len(tfe),
        n_supported_ge1_exp=s1,frac_supported_ge1=round(s1/n,4),
        n_supported_ge2_exp=s2,frac_supported_ge2=round(s2/n,4),
        analytic_expected_ge1=round(exp1,1),analytic_expected_ge2=round(exp2,1),
        null_mean_ge1=round(m1,1),null_sd_ge1=round(sd1,2),null_p_ge1=round(p1,4),
        null_mean_ge2=round(m2,1),null_sd_ge2=round(sd2,2),null_p_ge2=round(p2,4),
        null_reps=reps,gene_universe=uni,
        fold_enrichment_ge2=round((s2/n)/(m2/n),3) if m2>0 else None))
    print("window %skb: %d/%d testable edges; >=1 exp %d (%.1f%%); >=2 exp %d (%.1f%%); "
          "null(%d reps) >=2 = %.1f +/- %.2f ; p=%.4f ; fold=%.2f"%(
        dist,n,len(tfe),s1,100*s1/n,s2,100*s2/n,reps,m2,sd2,p2,(s2/n)/(m2/n) if m2>0 else float('nan')))

# ---------------- write ----------------
with open(OUT+"/chip_evidence.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=["TF","target","window_kb","in_table","n_experiments",
        "n_experiments_bound","frac_bound","mean_MACS2_all_exp","max_MACS2",
        "n_cell_types_bound","cell_types_bound","chip_support_ge1","chip_support_ge2"])
    w.writeheader(); w.writerows(rows)
print("WROTE chip_evidence.csv rows=",len(rows))
with open(OUT+"/chip_atlas_tf_inventory_v2d.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(inv[0].keys())); w.writeheader(); w.writerows(inv)
with open(OUT+"/chip_collagen_reverse_chipatlas.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(reverse[0].keys())); w.writeheader()
    w.writerows(sorted(reverse,key=lambda r:(r["gene"],r["window_kb"],-r["n_experiments_bound"])))
with open(OUT+"/chip_edge_support_summary_v2d.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
print("WROTE inventory / reverse / summary")
