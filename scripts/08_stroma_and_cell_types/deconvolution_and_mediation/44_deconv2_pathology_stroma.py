#!/usr/bin/env python3
"""44_deconv2_pathology_stroma.py

Retrieve the PATHOLOGIST-SCORED slide annotations for every TCGA-BRCA case from the
GDC API (percent_stromal_cells, percent_tumor_cells, percent_lymphocyte_infiltration,
percent_necrosis, percent_normal_cells).

Why this matters for the deconvolution question: every other mediator in this analysis
is derived from the same bulk RNA matrix as COL1A1 itself, so one can argue that
the mediation is a within-transcriptome tautology. The slide scores are read off
haematoxylin-eosin sections by a pathologist and share no measurement channel with the
RNA-seq. Together with the ABSOLUTE DNA-based purity estimate they give two mediators
that cannot be circular with COL1A1 mRNA.

Output: data/deconv/tcga_brca_slide_pathology.tsv  (one row per slide)
        results/v3/deconv_pathology_stroma_summary.csv
"""
import json, sys, urllib.request, urllib.parse, os
import csv

BASE = "/path/to/revision"
LOG = os.path.join(BASE, "logs/v3/44_deconv2_pathology.log")
_lf = open(LOG, "w")
def logf(*a):
    m = " ".join(str(x) for x in a)
    print(m); _lf.write(m + "\n"); _lf.flush()

filt = {"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-BRCA"]}}
params = {
    "filters": json.dumps(filt),
    "expand": "samples,samples.portions,samples.portions.slides",
    "fields": "submitter_id",
    "size": "2000",
    "format": "JSON",
}
url = "https://api.gdc.cancer.gov/cases?" + urllib.parse.urlencode(params)
logf("querying GDC cases endpoint for TCGA-BRCA slide annotations ...")
with urllib.request.urlopen(url, timeout=300) as r:
    d = json.load(r)
hits = d["data"]["hits"]
logf("cases returned:", len(hits), "of total", d["data"]["pagination"]["total"])

PCT = ["percent_stromal_cells", "percent_tumor_cells", "percent_tumor_nuclei",
       "percent_normal_cells", "percent_necrosis", "percent_lymphocyte_infiltration",
       "percent_monocyte_infiltration", "percent_neutrophil_infiltration"]

rows = []
for c in hits:
    case = c.get("submitter_id")
    for s in c.get("samples", []):
        samp = s.get("submitter_id")
        stype = s.get("sample_type")
        for p in s.get("portions", []):
            for sl in p.get("slides", []):
                rec = {"case": case, "sample": samp, "sample_type": stype,
                       "slide": sl.get("submitter_id")}
                for k in PCT:
                    rec[k] = sl.get(k)
                rows.append(rec)
logf("slide records:", len(rows))

out = os.path.join(BASE, "data/deconv/tcga_brca_slide_pathology.tsv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["case", "sample", "sample_type", "slide"] + PCT,
                       delimiter="\t")
    w.writeheader()
    for r_ in rows:
        w.writerow(r_)
logf("WROTE", out)

# quick summary
n_with_stroma = sum(1 for r_ in rows if r_["percent_stromal_cells"] is not None)
logf("slides with a percent_stromal_cells value:", n_with_stroma)
kinds = {}
for r_ in rows:
    sl = r_["slide"] or ""
    tag = sl.split("-")[-1][:2] if "-" in sl else "?"
    kinds[tag] = kinds.get(tag, 0) + 1
logf("slide-type suffixes:", sorted(kinds.items(), key=lambda x: -x[1])[:10])

# how many of our 1,097 tumours get a BS (frozen, sequenced portion) slide score
ss = [l.strip() for l in open(os.path.join(BASE, "cache/v3/deconv/samples_primary_tumour.txt"))]
sset = set(ss)
by_sample = {}
for r_ in rows:
    sl = r_["slide"] or ""
    if r_["percent_stromal_cells"] is None:
        continue
    short = (r_["sample"] or "")[:15]
    if short not in sset:
        continue
    tag = "BS" if "-BS" in sl else ("TS" if "-TS" in sl else ("DX" if "-DX" in sl else "other"))
    by_sample.setdefault(short, []).append((tag, r_["percent_stromal_cells"]))
logf("of our", len(ss), "primary tumours,", len(by_sample), "have >=1 scored slide")
nbs = sum(1 for v in by_sample.values() if any(t == "BS" for t, _ in v))
logf("  ... of which", nbs, "have a BS (frozen/sequenced portion) slide")

with open(os.path.join(BASE, "results/v3/deconv_pathology_stroma_summary.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["metric", "value"])
    w.writerow(["cases_returned", len(hits)])
    w.writerow(["slide_records", len(rows)])
    w.writerow(["slides_with_percent_stromal_cells", n_with_stroma])
    w.writerow(["our_primary_tumours", len(ss)])
    w.writerow(["our_tumours_with_any_scored_slide", len(by_sample)])
    w.writerow(["our_tumours_with_BS_slide_scored", nbs])
logf("DONE 44")
