#!/usr/bin/env python3
"""Verify every PMID cited in node_literature_text.py against live NCBI esummary records.
Writes results/v5/node_literature_verified.csv. Any PMID that cannot be retrieved, or whose
record is empty, is reported as FAILED and must be removed before the compendium is built."""
import json, re, sys, time, urllib.parse, urllib.request, csv, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import node_literature_text as L

EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "&tool=ffl_revision&email=your.email@example.org"
OUT = "/path/to/revision/results/v5/node_literature_verified.csv"

def get(u, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(u, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:
            if i == tries - 1: raise
            time.sleep(2 + 2 * i)

pmids = sorted({p for v in L.LIT.values() for p in re.findall(r"PMID:(\d+)", v)})
print("PMIDs to verify:", len(pmids))
rows, failed = [], []
for i in range(0, len(pmids), 20):
    chunk = pmids[i:i+20]
    res = get(EUT + "esummary.fcgi?db=pubmed&retmode=json&id=%s%s" % (",".join(chunk), TOOL))["result"]
    time.sleep(0.4)
    for p in chunk:
        s = res.get(p)
        if not isinstance(s, dict) or not s.get("title"):
            failed.append(p); continue
        auth = s.get("authors") or []
        rows.append(dict(pmid=p,
                         first_author=(auth[0]["name"] if auth else ""),
                         last_author=s.get("lastauthor", ""),
                         n_authors=len(auth),
                         year=s.get("pubdate", "")[:4],
                         journal=s.get("source", ""),
                         title=s.get("title", "").rstrip("."),
                         volume=s.get("volume", ""), pages=s.get("pages", ""),
                         doi=next((x["value"] for x in s.get("articleids", [])
                                   if x.get("idtype") == "doi"), ""),
                         pubtype="; ".join(s.get("pubtype", []))))
        print(f"  OK {p} {rows[-1]['year']} {rows[-1]['journal'][:24]:24s} {rows[-1]['title'][:70]}")
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("verified:", len(rows), "FAILED:", failed)
if failed:
    sys.exit("PMID verification FAILED for: " + ", ".join(failed))
