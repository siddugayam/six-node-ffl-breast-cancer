#!/usr/bin/env python3
"""Systematic GEO (E-utilities) + ArrayExpress/BioStudies search for breast cancer
miRNA datasets. Writes results/v6/cohorts_geo_search.csv"""
import json, os, re, sys, time, urllib.parse, urllib.request

OUT = "/path/to/revision/results/v6"
CACHE = "/path/to/revision/cache/v6/cohorts"
os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = "your.email@example.org"
TOOL = "mirna_ffl_revision"

def get(url, tries=5, sleep=0.4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mirna-ffl/1.0 (%s)" % EMAIL})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            sys.stderr.write("retry %d %s : %s\n" % (i, url[:120], e))
            time.sleep(2 + 2*i)
    return None

def esearch(db, term, retmax=2000):
    u = EU + "esearch.fcgi?" + urllib.parse.urlencode(
        dict(db=db, term=term, retmax=retmax, retmode="json", tool=TOOL, email=EMAIL))
    t = get(u)
    if t is None: return []
    j = json.loads(t)
    return j.get("esearchresult", {}).get("idlist", [])

def esummary(db, ids, chunk=200):
    out = []
    for i in range(0, len(ids), chunk):
        sub = ids[i:i+chunk]
        u = EU + "esummary.fcgi?" + urllib.parse.urlencode(
            dict(db=db, id=",".join(sub), retmode="json", tool=TOOL, email=EMAIL))
        t = get(u)
        if t is None: continue
        try: j = json.loads(t)
        except Exception: continue
        res = j.get("result", {})
        for k in res.get("uids", []):
            out.append(res[k])
        time.sleep(0.4)
    return out

# ---- systematic query set (not a guessed accession list) -----------------
QUERIES = {
 # broad: breast + mirna, GSE series only, human
 "Q1_breast_mirna_series":
   '("breast"[All Fields]) AND (microRNA[All Fields] OR miRNA[All Fields] OR "non-coding RNA profiling by array"[DataSet Type] OR "non-coding RNA profiling by high throughput sequencing"[DataSet Type]) AND "Homo sapiens"[Organism] AND GSE[Entry Type]',
 # dataset-type restricted, breast
 "Q2_breast_ncRNA_datasettype":
   '("breast"[All Fields]) AND ("non-coding RNA profiling by array"[DataSet Type] OR "non-coding RNA profiling by high throughput sequencing"[DataSet Type]) AND "Homo sapiens"[Organism] AND GSE[Entry Type]',
 # breast + mirna + survival/outcome words
 "Q3_breast_mirna_outcome":
   '("breast"[All Fields]) AND (microRNA[All Fields] OR miRNA[All Fields]) AND (survival[All Fields] OR prognosis[All Fields] OR outcome[All Fields] OR relapse[All Fields] OR recurrence[All Fields] OR "follow-up"[All Fields] OR "disease-free"[All Fields] OR response[All Fields]) AND "Homo sapiens"[Organism] AND GSE[Entry Type]',
 # serum / plasma / circulating (Part C)
 "Q4_breast_mirna_circulating":
   '("breast"[All Fields]) AND (microRNA[All Fields] OR miRNA[All Fields]) AND (serum[All Fields] OR plasma[All Fields] OR circulating[All Fields] OR "liquid biopsy"[All Fields] OR exosom*[All Fields]) AND "Homo sapiens"[Organism] AND GSE[Entry Type]',
 # neoadjuvant / chemotherapy response + mirna
 "Q5_breast_mirna_neoadjuvant":
   '("breast"[All Fields]) AND (microRNA[All Fields] OR miRNA[All Fields]) AND (neoadjuvant[All Fields] OR chemotherapy[All Fields] OR "pathological complete response"[All Fields] OR pCR[All Fields]) AND "Homo sapiens"[Organism] AND GSE[Entry Type]',
}
# platform-anchored search: find miRNA GPLs then find breast GSEs on them
GPL_QUERY = '(miRNA[All Fields] OR microRNA[All Fields]) AND "Homo sapiens"[Organism] AND GPL[Entry Type]'

records = {}
qhits = {}
for name, q in QUERIES.items():
    ids = esearch("gds", q, retmax=3000)
    qhits[name] = len(ids)
    sys.stderr.write("%s -> %d uids\n" % (name, len(ids)))
    time.sleep(0.5)
    for s in esummary("gds", ids):
        acc = s.get("accession", "")
        if not acc.startswith("GSE"): continue
        r = records.setdefault(acc, {"accession": acc, "title": s.get("title",""),
             "summary": (s.get("summary","") or "")[:2000],
             "n_samples": s.get("n_samples", ""), "gdstype": s.get("gdstype",""),
             "taxon": s.get("taxon",""), "pdat": s.get("pdat",""),
             "gpl": s.get("gpl",""), "queries": set(), "suppfile": s.get("suppfile","")})
        r["queries"].add(name)

sys.stderr.write("total unique GSE from queries: %d\n" % len(records))
with open(os.path.join(CACHE, "geo_raw_records.json"), "w") as f:
    json.dump({k: {kk: (sorted(vv) if isinstance(vv, set) else vv) for kk, vv in v.items()}
               for k, v in records.items()}, f)
with open(os.path.join(CACHE, "query_hits.json"), "w") as f:
    json.dump(qhits, f, indent=1)
print(json.dumps(qhits, indent=1))
