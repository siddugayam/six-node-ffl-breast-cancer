#!/usr/bin/env python3
"""Genome-wide promoter null for ReMap edge support: for each network TF_target
edge, replace the real target promoter by a randomly drawn promoter from 600
random RefSeqSelect genes that are NOT network nodes."""
import csv, json, os, collections, random, math
BASE="/path/to/revision"
def load_dir(d):
    out={}
    for fn in os.listdir(d):
        if not fn.endswith(".bed"): continue
        g=fn[:-4]; pk=[]
        for line in open(os.path.join(d,fn)):
            q=line.rstrip("\n").split("\t")
            if len(q)<11: continue
            pk.append((int(q[6]),q[3],q[9]))   # summit, dataset, TF
        out[g]=pk
    return out
netpk=load_dir(BASE+"/cache/direct/remap/regions")
bgpk =load_dir(BASE+"/cache/direct/remap/bg")
nettss=json.load(open(BASE+"/cache/direct/node_tss_hg38.json"))
bgtss =json.load(open(BASE+"/cache/direct/bg_tss.json"))
print("network promoters:",len(netpk)," background promoters:",len(bgpk))
print("background peak records:",sum(len(v) for v in bgpk.values()))

def index(pk,tss,getT):
    idx={}
    for g,peaks in pk.items():
        t=getT(g); d={}
        for w in (1000,5000):
            m=collections.defaultdict(set)
            for summit,ds,tf in peaks:
                if abs(summit-t)<=w: m[tf].add(ds)
            d[w]={tf:len(s) for tf,s in m.items()}
        idx[g]=d
    return idx
NET=index(netpk,nettss,lambda g:nettss[g][1])
BG =index(bgpk ,bgtss ,lambda g:bgtss[g][1])
print("median distinct TFs per network promoter (+/-1kb):",
      sorted(len(NET[g][1000]) for g in NET)[len(NET)//2])
print("median distinct TFs per background promoter (+/-1kb):",
      sorted(len(BG[g][1000]) for g in BG)[len(BG)//2])

edges=list(csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t"))
tfe=[e for e in edges if e["edge_type"]=="TF_target"]
ALIAS={"CTGF":"CCN2","CYR61":"CCN1","MKL1":"MRTFA"}
catalog={tf for g in NET for tf in NET[g][5000]} | {tf for g in BG for tf in BG[g][5000]}
test=[e for e in tfe if e["source"] in catalog and (e["target"] in NET or ALIAS.get(e["target"]) in NET)]
print("edges testable (TF in ReMap catalogue, target promoter fetched):",len(test),"of",len(tfe))
bgn=sorted(BG); random.seed(20260909); reps=1000
out=[]
for w in (1000,5000):
    for minds in (1,2):
        obs=0
        for e in test:
            g=e["target"] if e["target"] in NET else ALIAS[e["target"]]
            if NET[g][w].get(e["source"],0)>=minds: obs+=1
        null=[]
        for _ in range(reps):
            c=0
            for e in test:
                g=bgn[random.randrange(len(bgn))]
                if BG[g][w].get(e["source"],0)>=minds: c+=1
            null.append(c)
        m=sum(null)/reps; sd=math.sqrt(sum((x-m)**2 for x in null)/(reps-1))
        p=(sum(1 for x in null if x>=obs)+1)/(reps+1)
        z=(obs-m)/sd if sd else float("nan")
        out.append(dict(window_bp=w,min_datasets=minds,n_edges=len(test),observed=obs,
            obs_frac=round(obs/len(test),4),null_mean=round(m,1),null_sd=round(sd,2),
            null_frac=round(m/len(test),4),Z=round(z,2),
            fold=round(obs/m,3) if m else None,p_emp=round(p,4),reps=reps,
            null_type="genome-wide random promoter (600 non-network RefSeqSelect genes)"))
        print("+/-%dbp, >=%d dataset(s): obs %d/%d (%.1f%%) vs null %.1f +/- %.2f (%.1f%%), "
              "Z=%.2f fold=%.3f p=%.4f  [%d reps]"%(w,minds,obs,len(test),100*obs/len(test),
              m,sd,100*m/len(test),z,obs/m if m else float('nan'),p,reps))
with open(BASE+"/results/v2/chip_edge_null_genomewide.csv","w",newline="") as fh:
    wr=csv.DictWriter(fh,fieldnames=list(out[0].keys())); wr.writeheader(); wr.writerows(out)
print("WROTE chip_edge_null_genomewide.csv")
# promiscuity of the collagen promoters vs background
for g in ("COL1A1","COL3A1"):
    for w in (1000,5000):
        n=len(NET[g][w]); dist=sorted(len(BG[x][w]) for x in BG)
        rank=sum(1 for x in dist if x<n)/len(dist)
        print("%s +/-%dbp: %d distinct ReMap TFs = %.1fth percentile of random promoters "
              "(background median %d)"%(g,w,n,100*rank,dist[len(dist)//2]))
