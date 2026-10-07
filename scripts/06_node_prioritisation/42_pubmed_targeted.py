#!/usr/bin/env python3
"""Targeted PubMed lookups for specific statements; prints PMID/year/journal/title for selection."""
import json, time, urllib.parse, urllib.request, sys
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL="&tool=ffl_revision&email=your.email@example.org"
def get(u,tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(u,timeout=45) as r: return json.loads(r.read().decode())
        except Exception:
            if i==tries-1: raise
            time.sleep(2+2*i)
def look(term,n=5,sort="relevance"):
    u=EUT+"esearch.fcgi?db=pubmed&retmode=json&retmax=%d&sort=%s&term=%s%s"%(n,sort,urllib.parse.quote(term),TOOL)
    time.sleep(0.4); ids=get(u)["esearchresult"]["idlist"]
    if not ids: return []
    u2=EUT+"esummary.fcgi?db=pubmed&retmode=json&id=%s%s"%(",".join(ids),TOOL)
    time.sleep(0.4); res=get(u2)["result"]
    out=[]
    for i in ids:
        s=res.get(i)
        if isinstance(s,dict):
            out.append((i,s.get("pubdate","")[:4],s.get("source",""),s.get("title","")))
    return out
QUERIES=json.load(open(sys.argv[1]))
for label,term in QUERIES:
    print(f"##### {label}\n    q: {term}")
    for pmid,yr,jr,ti in look(term):
        print(f"    {pmid} {yr} {jr[:26]:26s} | {ti[:120]}")
    print()
