#!/usr/bin/env python3
"""Triage the GEO search hits: keep human breast miRNA series with n>=40,
then probe each for outcome fields."""
import json, os, re, sys, time, urllib.parse, urllib.request, csv

CACHE = "/path/to/revision/cache/v6/cohorts"
OUT = "/path/to/revision/results/v6"
recs = json.load(open(os.path.join(CACHE, "geo_raw_records.json")))
print("loaded", len(recs))

MIRNA_RE = re.compile(r"mirna|microrna|mir-|small rna|smallrna|ncrna|non-coding", re.I)
BREAST_RE = re.compile(r"breast|mammary", re.I)
OUTCOME_RE = re.compile(r"surviv|prognos|outcome|relapse|recurrenc|follow.?up|disease.free|"
                        r"metastasis.free|pathologic(al)? complete response|\bpCR\b|"
                        r"response to|chemotherapy|neoadjuvant|event.free|distant", re.I)

rows = []
for acc, r in recs.items():
    txt = (r["title"] or "") + " || " + (r["summary"] or "")
    n = r.get("n_samples", 0)
    try: n = int(n)
    except Exception: n = 0
    gdst = r.get("gdstype","")
    is_mirna = bool(MIRNA_RE.search(txt)) or ("non-coding" in gdst.lower())
    is_breast = bool(BREAST_RE.search(txt))
    rows.append(dict(accession=acc, n_samples=n, gdstype=gdst, gpl=r.get("gpl",""),
        pdat=r.get("pdat",""), title=(r["title"] or "")[:250],
        mirna_kw=is_mirna, breast_kw=is_breast,
        outcome_kw=bool(OUTCOME_RE.search(txt)),
        queries=";".join(r["queries"]), summary=(r["summary"] or "")[:600]))

rows.sort(key=lambda x: -x["n_samples"])
cand = [r for r in rows if r["mirna_kw"] and r["breast_kw"] and r["n_samples"] >= 40]
print("breast+mirna+n>=40:", len(cand))
print("of those with outcome keyword:", sum(1 for r in cand if r["outcome_kw"]))
with open(os.path.join(CACHE, "triage_all.json"), "w") as f: json.dump(rows, f)
with open(os.path.join(CACHE, "triage_candidates.json"), "w") as f: json.dump(cand, f)
for r in cand[:70]:
    print("%-12s n=%-5d out=%-5s %s | GPL %s" % (r["accession"], r["n_samples"],
          r["outcome_kw"], r["title"][:110], r["gpl"][:40]))
