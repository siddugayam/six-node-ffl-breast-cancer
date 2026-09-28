#!/usr/bin/env python3
"""C-i) SigCom LINCS: does knocking down / over-expressing NFKB1, RELA, SP1, ETS1
change COL1A1 / COL3A1 in the L1000 CD-coefficient signatures?"""
import json, os, sys, gzip, io, urllib.parse, time
from concurrent.futures import ThreadPoolExecutor
import urllib.request

MDAPI = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
SCRATCH = "/path/to/scratch/lincs_sigs"
OUT = "/path/to/revision/results/multiomics"
os.makedirs(SCRATCH, exist_ok=True); os.makedirs(OUT, exist_ok=True)
TFS = ["NFKB1","RELA","SP1","ETS1"]
READOUT = ["COL1A1","COL3A1"]

def get(url, params=None, tries=3):
    if params: url = url + "?" + urllib.parse.urlencode(params)
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: return r.read()
        except Exception as e:
            if i == tries-1: raise
            time.sleep(3)

meta = {}
for tf in TFS:
    f = json.dumps({"where":{"meta.pert_name":tf},"limit":1000,"fields":["id","meta"]})
    d = json.loads(get(MDAPI+"/signatures", {"filter": f}))
    meta[tf] = d
    print(f"{tf}: {len(d)} signatures in SigCom LINCS", flush=True)

records = []
for tf, sigs in meta.items():
    for s in sigs:
        m = s["meta"]
        records.append(dict(tf=tf, sig_id=s["id"], local_id=m.get("local_id"),
                            pert_type=m.get("pert_type"), cell_line=m.get("cell_line"),
                            pert_time=m.get("pert_time"), pert_dose=m.get("pert_dose"),
                            url=m.get("persistent_id"), filename=m.get("filename")))
print("total signature records:", len(records), flush=True)

def fetch(r):
    if not r["url"]: return None
    fn = os.path.join(SCRATCH, r["local_id"].replace("/","_") + (".tsv.gz" if r["url"].endswith(".gz") else ".tsv"))
    if not os.path.exists(fn) or os.path.getsize(fn) == 0:
        try:
            raw = get(r["url"])
        except Exception as e:
            print("FAIL", r["local_id"], e, flush=True); return None
        with open(fn,"wb") as fh: fh.write(raw)
    return fn

with ThreadPoolExecutor(max_workers=12) as ex:
    paths = list(ex.map(fetch, records))

import csv
rows = []
n_ok = 0
for r, p in zip(records, paths):
    if p is None: continue
    op = gzip.open if p.endswith(".gz") else open
    try:
        with op(p, "rt") as fh:
            rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
            vals = {}
            for line in rd:
                if len(line) < 2: continue
                try: vals[line[0]] = float(line[1])
            	# noqa
                except ValueError: pass
    except Exception as e:
        print("PARSE FAIL", p, e, flush=True); continue
    if not vals: continue
    n_ok += 1
    ordered = sorted(vals.values())
    n = len(ordered)
    import bisect
    out = dict(r); out.pop("url"); out["n_genes"] = n
    out["value_col"] = hdr[1] if len(hdr)>1 else "?"
    for g in READOUT + ["NFKB1","RELA","SP1","ETS1"]:
        v = vals.get(g)
        out[g+"_value"] = v
        out[g+"_pctile"] = (bisect.bisect_left(ordered, v)/n) if v is not None else None
    rows.append(out)

print("signatures parsed:", n_ok, flush=True)
keys = sorted({k for r in rows for k in r})
order = ["tf","pert_type","local_id","cell_line","pert_time","pert_dose","filename","n_genes","value_col"] + \
        [k for k in keys if k.endswith("_value") or k.endswith("_pctile")] + ["sig_id"]
order = [k for k in order if k in keys] + [k for k in keys if k not in order]
with open(os.path.join(OUT,"lincs_TF_signatures_collagen.csv"),"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=order); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote lincs_TF_signatures_collagen.csv rows:", len(rows), flush=True)
