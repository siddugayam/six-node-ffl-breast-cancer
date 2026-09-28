#!/usr/bin/env python3
"""Resolve breast-cancer EFO trait hierarchy via OLS4 and export the term list."""
import json, time, urllib.parse, sys
import requests

OUT = "/path/to/revision/data/gwas"
S = requests.Session()
S.headers.update({"User-Agent": "miRNA-FFL-revision/1.0"})

def enc(iri):
    return urllib.parse.quote(urllib.parse.quote(iri, safe=""), safe="")

def term(iri):
    r = S.get(f"https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms/{enc(iri)}", timeout=60)
    return r.json() if r.status_code == 200 else None

def descendants(iri):
    out, page = [], 0
    while True:
        u = (f"https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms/{enc(iri)}"
             f"/descendants?size=200&page={page}")
        r = S.get(u, timeout=90)
        if r.status_code != 200:
            break
        j = r.json()
        terms = j.get("_embedded", {}).get("terms", [])
        out.extend(terms)
        pg = j.get("page", {})
        if page >= pg.get("totalPages", 1) - 1:
            break
        page += 1
        time.sleep(0.2)
    return out

roots = {
    "EFO_0000305": "http://www.ebi.ac.uk/efo/EFO_0000305",
    "MONDO_0004989": "http://purl.obolibrary.org/obo/MONDO_0004989",
}
rows = []
for sf, iri in roots.items():
    t = term(iri)
    if t is None:
        print(f"{sf}: NOT RESOLVED", file=sys.stderr); continue
    print(f"root {sf}: label={t.get('label')} obsolete={t.get('is_obsolete')}")
    rows.append({"short_form": t.get("short_form"), "label": t.get("label"),
                 "iri": t.get("iri"), "obsolete": t.get("is_obsolete"),
                 "relation": "root", "root": sf})
    if t.get("is_obsolete"):
        continue
    for d in descendants(iri):
        rows.append({"short_form": d.get("short_form"), "label": d.get("label"),
                     "iri": d.get("iri"), "obsolete": d.get("is_obsolete"),
                     "relation": "descendant", "root": sf})

import pandas as pd
df = pd.DataFrame(rows).drop_duplicates(subset=["short_form"])
df.to_csv(f"{OUT}/efo_breast_carcinoma_terms.csv", index=False)
print("terms:", len(df))
print(df["label"].head(40).tolist())
