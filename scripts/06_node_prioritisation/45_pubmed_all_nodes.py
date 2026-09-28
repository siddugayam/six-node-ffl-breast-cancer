#!/usr/bin/env python3
"""PubMed record counts for all 587 network nodes, to test whether network centrality
tracks how heavily a node has been studied (the core argument of PRIORITISATION.md)."""
import json, time, urllib.parse, urllib.request, csv, os, re
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL="&tool=ffl_revision&email=your.email@example.org"
ROOT="/path/to/revision"
def get(u,tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(u,timeout=60) as r: return json.loads(r.read().decode())
        except Exception:
            if i==tries-1: raise
            time.sleep(2+2*i)
def count(term):
    u=EUT+"esearch.fcgi?db=pubmed&retmode=json&retmax=0&term=%s%s"%(urllib.parse.quote(term),TOOL)
    time.sleep(0.36); return int(get(u)["esearchresult"]["count"])
nodes=[]
with open(os.path.join(ROOT,"data/canonical_nodes.tsv")) as f:
    rd=csv.DictReader(f,delimiter="\t")
    for r in rd: nodes.append((r["name"],r["type"]))
BC='"Breast Neoplasms"[Mesh]'
def term(name,typ):
    if typ=="miRNA":
        s=name.replace("hsa-","")
        return f'("{s}"[tiab] OR "micro{s}"[tiab])'
    return f'"{name}"[tiab]'
out=os.path.join(ROOT,"results/v5/network_literature_volume.csv")
done=set()
if os.path.exists(out):
    with open(out) as f:
        for r in csv.DictReader(f): done.add(r["name"])
mode="a" if done else "w"
with open(out,mode,newline="") as f:
    w=csv.DictWriter(f,fieldnames=["name","type","pubmed_total","pubmed_breast"])
    if not done: w.writeheader()
    for k,(n,t) in enumerate(nodes,1):
        if n in done: continue
        tm=term(n,t)
        try:
            tot=count(tm); bc=count(tm+" AND "+BC)
        except Exception as e:
            print("FAIL",n,e,flush=True); continue
        w.writerow(dict(name=n,type=t,pubmed_total=tot,pubmed_breast=bc)); f.flush()
        if k%25==0: print(f"{k}/{len(nodes)} {n} {tot}/{bc}",flush=True)
print("done")
