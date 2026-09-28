#!/usr/bin/env python3
"""Query DGIdb v5 GraphQL for the 20 protein-coding prioritised nodes.
Writes per-interaction and per-gene summary tables to results/v5/."""
import json, time, sys, urllib.request
import pandas as pd

GENES = ['E2F1','EZH2','GATA3','BRCA1','JUN','EGR2','ESR1','SREBF1','DNMT1','E2F3',
         'CCND2','COL1A1','STAT5A','MYBL2','FN1','PDGFRB','MET','CXCL12','MMP14','PLAU']
URL = "https://dgidb.org/api/graphql"
Q = """{genes(names:[%s]){nodes{name longName interactions{
  drug{name conceptId approved antiNeoplastic immunotherapy}
  interactionScore
  interactionTypes{type directionality}
  sources{sourceDbName}}}}}"""

rows = []
found = set()
for g in GENES:
    q = Q % json.dumps(g)
    req = urllib.request.Request(URL, data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            r = json.load(urllib.request.urlopen(req, timeout=90)); break
        except Exception as e:
            sys.stderr.write(f"{g} attempt {attempt}: {e}\n"); time.sleep(5)
    else:
        sys.stderr.write(f"{g}: FAILED\n"); continue
    nodes = r.get("data", {}).get("genes", {}).get("nodes", [])
    for n in nodes:
        if n["name"] != g:
            continue
        found.add(g)
        for it in n["interactions"]:
            rows.append(dict(gene=g, long_name=n.get("longName"),
                             drug=it["drug"]["name"], concept_id=it["drug"]["conceptId"],
                             approved=it["drug"]["approved"],
                             anti_neoplastic=it["drug"]["antiNeoplastic"],
                             immunotherapy=it["drug"]["immunotherapy"],
                             interaction_score=it["interactionScore"],
                             interaction_types="|".join(sorted({t["type"] for t in it["interactionTypes"] if t["type"]})),
                             directionality="|".join(sorted({t["directionality"] for t in it["interactionTypes"] if t["directionality"]})),
                             n_sources=len(it["sources"]),
                             sources="|".join(sorted({s["sourceDbName"] for s in it["sources"]}))))
    time.sleep(1)

d = pd.DataFrame(rows)
d.to_csv("results/v5/node_druggability_dgidb_interactions.csv", index=False)

summ = []
for g in GENES:
    s = d[d.gene == g]
    if len(s) == 0:
        summ.append(dict(gene=g, in_dgidb=(g in found), n_drugs=0, n_approved=0,
                         n_antineoplastic=0, n_approved_antineoplastic=0,
                         interaction_types="", top_approved_antineoplastic="",
                         top_approved=""))
        continue
    ap = s[s.approved == True]
    an = s[s.anti_neoplastic == True]
    apan = s[(s.approved == True) & (s.anti_neoplastic == True)]
    def top(df, k=6):
        return "; ".join(df.sort_values("interaction_score", ascending=False)
                           .drug.drop_duplicates().head(k).tolist())
    summ.append(dict(gene=g, in_dgidb=True, n_drugs=s.drug.nunique(),
                     n_approved=ap.drug.nunique(), n_antineoplastic=an.drug.nunique(),
                     n_approved_antineoplastic=apan.drug.nunique(),
                     interaction_types="|".join(sorted({t for x in s.interaction_types for t in x.split("|") if t})),
                     top_approved_antineoplastic=top(apan), top_approved=top(ap)))
pd.DataFrame(summ).to_csv("results/v5/node_druggability_dgidb_summary.csv", index=False)
print(pd.DataFrame(summ).to_string(index=False))
