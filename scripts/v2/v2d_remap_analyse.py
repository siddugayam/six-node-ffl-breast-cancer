#!/usr/bin/env python3
"""ReMap 2022 (hg38, all 8,103 non-merged datasets) promoter occupancy.
Region-restricted reads of reMap2022.bb for TSS +/-5 kb of every gene/TF node.
Peak 'summit' = thickStart. Promoter binding = summit within TSS +/- w."""
import csv, json, os, collections, math, random
BASE="/path/to/revision"
RD=BASE+"/cache/direct/remap/regions"
OUT=BASE+"/results/v2"
tss=json.load(open(BASE+"/cache/direct/node_tss_hg38.json"))
edges=list(csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t"))
tfe=[e for e in edges if e["edge_type"]=="TF_target"]
ALIAS={"CTGF":"CCN2","CYR61":"CCN1","MKL1":"MRTFA"}

def load(gene):
    p=os.path.join(RD,gene+".bed")
    if not os.path.exists(p): return None
    out=[]
    with open(p) as f:
        for line in f:
            q=line.rstrip("\n").split("\t")
            if len(q)<11: continue
            out.append((int(q[1]),int(q[2]),q[3],int(q[6]),q[9],q[10]))  # s,e,dataset,summit,TF,biotype
    return out

# ---- global inventory of ReMap TFs seen in these promoter windows ----
allTF=set(); ndatasets=0; total=0
prom={}   # gene -> {w: {TF: [datasets, biotypes]}}
for g,(chrom,t,strand,tx) in tss.items():
    pk=load(g)
    if pk is None: continue
    total+=len(pk)
    d={}
    for w in (1000,5000):
        m=collections.defaultdict(lambda: [set(),set()])
        for s,e,ds,summit,tf,bt in pk:
            if abs(summit-t)<=w:
                m[tf][0].add(ds); m[tf][1].add(bt)
        d[w]={tf:(sorted(v[0]),sorted(v[1])) for tf,v in m.items()}
    prom[g]=d
    allTF.update(prom[g][5000])
print("genes with ReMap windows:",len(prom))
print("total ReMap peak records read (TSS+/-5kb):",total)
print("distinct ReMap TFs seen in any network promoter (+/-5kb):",len(allTF))

# ---- reverse query: TFs at COL1A1 / COL3A1 promoters ----
rev=[]
for gene in ("COL1A1","COL3A1","EZH2","VEGFA","FN1"):
    if gene not in prom: continue
    for w in (1000,5000):
        for tf,(dss,bts) in sorted(prom[gene][w].items(),key=lambda x:-len(x[1][0])):
            rev.append(dict(gene=gene,window_bp=w,TF=tf,n_datasets=len(dss),
                n_biotypes=len(bts),biotypes=";".join(bts[:20]),
                datasets=";".join(dss[:12])))
with open(OUT+"/remap_collagen_promoter_reverse.csv","w",newline="") as fh:
    w_=csv.DictWriter(fh,fieldnames=list(rev[0].keys())); w_.writeheader(); w_.writerows(rev)
print("WROTE remap_collagen_promoter_reverse.csv rows=",len(rev))
for gene in ("COL1A1","COL3A1"):
    for w in (1000,5000):
        n=len(prom[gene][w]); tot=sum(len(v[0]) for v in prom[gene][w].values())
        top=sorted(prom[gene][w].items(),key=lambda x:-len(x[1][0]))[:12]
        print("%s +/-%dbp : %d distinct TFs, %d datasets. Top: %s"%(gene,w,n,tot,
              ", ".join("%s(%d)"%(t,len(v[0])) for t,v in top)))
for tf in ("NFKB1","RELA","SP1","ETS1","SMAD3","TWIST1"):
    for gene in ("COL1A1","COL3A1"):
        for w in (1000,5000):
            v=prom.get(gene,{}).get(w,{}).get(tf)
            print("  %-7s -> %-7s +/-%4dbp : %s"%(tf,gene,w,
                  ("%d datasets in %s"%(len(v[0]),";".join(v[1][:6]))) if v else "NOT BOUND"))

# ---- edge-level ReMap support for the 770 TF_target edges ----
rows=[]
for e in tfe:
    tf=e["source"]; tgt=e["target"]
    g=tgt if tgt in prom else ALIAS.get(tgt,tgt)
    r=dict(TF=tf,target=tgt,tss_found=int(g in prom))
    for w in (1000,5000):
        v=prom.get(g,{}).get(w,{}).get(tf)
        r["remap_datasets_%d"%w]=len(v[0]) if v else 0
        r["remap_biotypes_%d"%w]=";".join(v[1][:15]) if v else ""
        r["remap_support_%d"%w]=int(bool(v))
        r["remap_support_ge2_%d"%w]=int(bool(v) and len(v[0])>=2)
    rows.append(r)
# restrict to TFs that exist in the ReMap catalogue at all
remap_catalog=allTF
rows_test=[r for r in rows if r["TF"] in remap_catalog and r["tss_found"]]
print("TF_target edges:",len(rows),"; edges whose TF appears anywhere in ReMap promoter data:",len(rows_test))
tf_in=sorted({r["TF"] for r in rows}); tf_ok=sorted({r["TF"] for r in rows_test})
print("network TFs:",len(tf_in),"present in ReMap:",len(tf_ok))
print("network TFs ABSENT from ReMap:",sorted(set(tf_in)-set(tf_ok)))
for w in (1000,5000):
    s1=sum(r["remap_support_%d"%w] for r in rows_test)
    s2=sum(r["remap_support_ge2_%d"%w] for r in rows_test)
    print("ReMap +/-%dbp : %d/%d edges bound (>=1 dataset) = %.1f%% ; >=2 datasets %d = %.1f%%"%(
        w,s1,len(rows_test),100*s1/len(rows_test),s2,100*s2/len(rows_test)))

# ---- matched null: shuffle target assignment within the set of network promoters ----
random.seed(20260909)
genes=[g for g in prom]
reps=1000
for w in (1000,5000):
    obs=sum(r["remap_support_ge2_%d"%w] for r in rows_test)
    null=[]
    for _ in range(reps):
        c=0
        for r in rows_test:
            g=genes[random.randrange(len(genes))]
            v=prom[g][w].get(r["TF"])
            if v and len(v[0])>=2: c+=1
        null.append(c)
    m=sum(null)/reps; sd=math.sqrt(sum((x-m)**2 for x in null)/(reps-1))
    p=(sum(1 for x in null if x>=obs)+1)/(reps+1)
    print("ReMap null +/-%dbp (%d reps, targets reshuffled over the %d network promoters): "
          "obs %d vs %.1f +/- %.2f, Z=%.2f, fold=%.3f, p=%.4f"%(
          w,reps,len(genes),obs,m,sd,(obs-m)/sd if sd else float('nan'),obs/m if m else float('nan'),p))

with open(OUT+"/remap_edge_support.csv","w",newline="") as fh:
    w_=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w_.writeheader(); w_.writerows(rows)
print("WROTE remap_edge_support.csv rows=",len(rows))
