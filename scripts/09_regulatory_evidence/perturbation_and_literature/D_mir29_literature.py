#!/usr/bin/env python3
"""D) Literature evidence (PubMed) for miR-29 family perturbation -> collagen genes."""
import json, urllib.parse, urllib.request, time, csv, os, re
OUT = "/path/to/revision/results/multiomics"
SCR = "/path/to/scratch"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get(url, params, tries=3):
    u = url + "?" + urllib.parse.urlencode(params)
    for i in range(tries):
        try:
            with urllib.request.urlopen(u, timeout=90) as r: return r.read().decode("utf-8","replace")
        except Exception as e:
            if i == tries-1: raise
            time.sleep(3)

QUERIES = {
 "miR29_COL1A1_direct": '("miR-29"[tiab] OR "microRNA-29"[tiab] OR "miRNA-29"[tiab]) AND ("COL1A1"[tiab] OR "collagen type I alpha 1"[tiab] OR "alpha1(I) collagen"[tiab])',
 "miR29_COL3A1_direct": '("miR-29"[tiab] OR "microRNA-29"[tiab] OR "miRNA-29"[tiab]) AND ("COL3A1"[tiab] OR "collagen type III"[tiab] OR "collagen III"[tiab])',
 "miR29_collagen_perturbation": '("miR-29"[tiab] OR "microRNA-29"[tiab]) AND (collagen[tiab]) AND (overexpression[tiab] OR mimic[tiab] OR inhibitor[tiab] OR antagomir[tiab] OR "knockdown"[tiab] OR transfect*[tiab])',
 "miR29_breast_cancer_collagen": '("miR-29"[tiab] OR "microRNA-29"[tiab]) AND (breast[tiab]) AND (collagen[tiab] OR "extracellular matrix"[tiab])',
 "miR29_luciferase_collagen": '("miR-29"[tiab] OR "microRNA-29"[tiab]) AND collagen[tiab] AND (luciferase[tiab] OR "3\'UTR"[tiab] OR "3 UTR"[tiab])',
}

allrows = []
for tag, q in QUERIES.items():
    r = json.loads(get(E+"esearch.fcgi", {"db":"pubmed","term":q,"retmax":80,"retmode":"json"}))
    ids = r["esearchresult"]["idlist"]
    print(f"{tag}: {r['esearchresult']['count']} hits, retrieved {len(ids)}", flush=True)
    if not ids: continue
    s = json.loads(get(E+"esummary.fcgi", {"db":"pubmed","id":",".join(ids),"retmode":"json"}))["result"]
    ab = get(E+"efetch.fcgi", {"db":"pubmed","id":",".join(ids),"retmode":"text","rettype":"abstract"})
    with open(os.path.join(SCR, f"pubmed_{tag}.txt"),"w") as fh: fh.write(ab)
    for u in s["uids"]:
        it = s[u]
        allrows.append(dict(query_tag=tag, pmid=u, title=it.get("title","").strip(),
                            journal=it.get("source"), year=(it.get("pubdate") or "")[:4],
                            doi=next((a["value"] for a in it.get("articleids",[]) if a["idtype"]=="doi"), "")))
    time.sleep(0.4)

# dedupe, keep all query tags
byid = {}
for r in allrows:
    if r["pmid"] in byid: byid[r["pmid"]]["query_tag"] += ";" + r["query_tag"]
    else: byid[r["pmid"]] = r
rows = sorted(byid.values(), key=lambda r: -int(r["year"] or 0))
with open(os.path.join(SCR,"mir29_pubmed_all.csv"),"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["pmid","year","journal","title","doi","query_tag"])
    w.writeheader()
    for r in rows: w.writerow({k:r[k] for k in w.fieldnames})
print("unique PMIDs:", len(rows), flush=True)
for r in rows[:60]:
    print(f"  {r['pmid']} ({r['year']}) {r['journal']}: {r['title'][:150]}")
