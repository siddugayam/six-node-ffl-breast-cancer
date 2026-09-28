#!/usr/bin/env python3
"""Are the network's 4,819 miRNA-target edges enriched for TargetScan sites over
random miRNA-gene pairs drawn from the same miRNA and gene universes?
This is the sequence-level mirror of the ChIP-seq null for the TF arm."""
import csv, json, os, random, math, collections
BASE="/path/to/revision"; TS=BASE+"/cache/direct/ts"
Q=json.load(open(BASE+"/cache/direct/ts_query_sets.json")); arm=Q["arm_map"]
seed_of={}
for line in open(TS+"/miR_Family_Info.txt"):
    p=line.rstrip("\n").split("\t")
    if len(p)<6 or p[2]!="9606": continue
    seed_of[p[3]]=p[1]
# pair-level site presence, from the network-gene-filtered site tables
site=collections.defaultdict(lambda:[0,0,0])   # (gene,mir) -> [n_sites, n_conserved, n_8mer]
def rd(path,cons):
    n=0
    with open(path) as f:
        next(f)
        for line in f:
            p=line.rstrip("\n").split("\t")
            if p[3]!="9606": continue
            v=site[(p[1],p[4])]; v[0]+=1; v[1]+=cons; v[2]+= (1 if p[5]=="3" else 0); n+=1
    return n
GEN=set(open(BASE+"/cache/direct/ts_genes.txt").read().split())
def rd_f(path,cons,filt):
    n=0
    with open(path) as f:
        next(f)
        for line in f:
            p=line.rstrip("\n").split("\t")
            if p[3]!="9606" or (filt and p[1] not in GEN): continue
            v=site[(p[1],p[4])]; v[0]+=1; v[1]+=cons; v[2]+= (1 if p[5]=="3" else 0); n+=1
    return n
a=rd_f(TS+"/Conserved_Site_Context_Scores.txt",1,True)
b=rd_f(TS+"/NC_sites_network.txt",0,False)
print("site rows loaded:",a+b,"over",len(site),"(gene,miRNA) pairs")
E=[e for e in csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t")
   if e["edge_type"]=="miRNA_target"]
mirs=sorted({e["source"] for e in E}); genes=sorted({e["target"] for e in E})
ALIAS={"CTGF":"CCN2","CYR61":"CCN1","MKL1":"MRTFA"}
TSG={g for (g,_) in site}
def res(g): return g if g in TSG else (ALIAS.get(g) if ALIAS.get(g) in TSG else g)
def score(lab,g):
    g=res(g); tot=[0,0,0]
    for m in arm[lab]:
        v=site.get((g,m))
        if v: tot=[tot[i]+v[i] for i in range(3)]
    return tot
obs=[score(e["source"],e["target"]) for e in E]
o_site=sum(1 for t in obs if t[0]>0); o_cons=sum(1 for t in obs if t[1]>0); o_8=sum(1 for t in obs if t[2]>0)
n=len(E)
print("observed: %d/%d edges with >=1 site (%.1f%%), %d conserved (%.1f%%), %d 8mer (%.1f%%)"%(
    o_site,n,100*o_site/n,o_cons,100*o_cons/n,o_8,100*o_8/n))
random.seed(20260909); reps=1000
cache={}
def sc(lab,g):
    k=(lab,g)
    if k not in cache: cache[k]=score(lab,g)
    return cache[k]
ns=[];nc=[];n8=[]
for _ in range(reps):
    a1=b1=c1=0
    for e in E:
        g=genes[random.randrange(len(genes))]
        t=sc(e["source"],g)
        if t[0]>0: a1+=1
        if t[1]>0: b1+=1
        if t[2]>0: c1+=1
    ns.append(a1);nc.append(b1);n8.append(c1)
out=[]
for nm,o,nl in (("any_site",o_site,ns),("conserved_site",o_cons,nc),("8mer_site",o_8,n8)):
    m=sum(nl)/reps; sd=math.sqrt(sum((x-m)**2 for x in nl)/(reps-1))
    p=(sum(1 for x in nl if x>=o)+1)/(reps+1)
    print("%-15s obs %d (%.1f%%) vs null %.1f +/- %.1f (%.1f%%) Z=%+.2f fold=%.3f p=%.4f [%d reps]"%(
        nm,o,100*o/n,m,sd,100*m/n,(o-m)/sd,o/m,p,reps))
    out.append(dict(metric=nm,observed=o,obs_frac=round(o/n,4),null_mean=round(m,1),
        null_sd=round(sd,2),null_frac=round(m/n,4),Z=round((o-m)/sd,2),fold=round(o/m,3),
        p_emp=round(p,4),reps=reps,n_edges=n,
        null_type="target gene reshuffled within the 363 network target genes, miRNA fixed"))
with open(BASE+"/results/v2/targetscan_site_null.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
print("WROTE targetscan_site_null.csv")
