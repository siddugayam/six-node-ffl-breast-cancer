#!/usr/bin/env python3
"""Third resource: ENCODE portal REST API - how many human TF ChIP-seq experiments exist for
NFKB1 / RELA / SP1 / ETS1, and in which biosamples (breast highlighted)."""
import json, urllib.request, urllib.parse, csv, os
OUT = "/path/to/revision/results/v2"
H = {"Accept": "application/json", "User-Agent": "revision-analysis"}
rows = []
for tf in ["NFKB1", "RELA", "SP1", "ETS1", "EGR1", "TFAP2A", "SMAD3", "GATA2", "TWIST1"]:
    q = ("https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq"
         f"&target.label={urllib.parse.quote(tf)}&replicates.library.biosample.donor.organism.scientific_name="
         "Homo+sapiens&status=released&limit=all&format=json")
    req = urllib.request.Request(q, headers=H)
    try:
        d = json.load(urllib.request.urlopen(req, timeout=120))
    except urllib.error.HTTPError as e:
        # the ENCODE portal answers a search with zero hits with HTTP 404
        if e.code == 404:
            rows.append(dict(TF=tf, n_experiments=0, error="404 = zero released experiments",
                             biosamples="", breast_biosamples=""))
            print(f"{tf}: 0 released human TF ChIP-seq experiments on ENCODE (HTTP 404 = empty search)")
            continue
        raise
    n = d.get("total", 0)
    bios = sorted({g.get("biosample_ontology", {}).get("term_name", "?") for g in d.get("@graph", [])})
    br = [b for b in bios if any(k in b.lower() for k in ("breast", "mcf", "mammary", "t47", "sk-br", "mda-mb"))]
    rows.append(dict(TF=tf, n_experiments=n, error="", biosamples=";".join(bios),
                     breast_biosamples=";".join(br)))
    print(f"{tf}: {n} released human TF ChIP-seq experiments on ENCODE; "
          f"{len(bios)} biosamples; breast: {br if br else 'NONE'}")
with open(os.path.join(OUT, "chip_encode_experiment_counts.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote chip_encode_experiment_counts.csv")
