#!/usr/bin/env python3
"""Fetch ENCODE hg38 peak calls for breast carcinoma lines, breast epithelium and
fibroblasts, to ask whether the COL1A1/COL3A1 promoters are accessible in
carcinoma cells vs fibroblasts."""
import json, os, urllib.request, urllib.parse, gzip, sys
B = "https://www.encodeproject.org"
CACHE = "/path/to/revision/cache/v6/atac/encode"
os.makedirs(CACHE, exist_ok=True)

def q(url):
    r = urllib.request.Request(url, headers={'Accept': 'application/json'})
    return json.load(urllib.request.urlopen(r, timeout=180))

def experiments(term, assay):
    url = (B + "/search/?type=Experiment&status=released&format=json&limit=100"
           "&field=accession&field=assay_title&field=biosample_ontology.term_name"
           "&field=biosample_summary&assay_title=" + urllib.parse.quote(assay) +
           "&biosample_ontology.term_name=" + urllib.parse.quote(term))
    try: d = q(url)
    except Exception: return []
    return d.get('@graph', [])

WANT = [("MCF-7", "ATAC-seq", 1), ("MCF-7", "DNase-seq", 1), ("T47D", "DNase-seq", 1),
        ("MCF 10A", "DNase-seq", 1), ("breast epithelium", "ATAC-seq", 3),
        ("IMR-90", "ATAC-seq", 1), ("fibroblast of mammary gland", "DNase-seq", 1),
        ("fibroblast of dermis", "DNase-seq", 1), ("fibroblast of lung", "DNase-seq", 2)]

rows = []
for term, assay, k in WANT:
    exps = experiments(term, assay)[:k]
    for e in exps:
        acc = e['accession']
        det = q(f"{B}/experiments/{acc}/?format=json")
        best = None
        for f in det.get('files', []):
            if (f.get('file_format') == 'bed' and f.get('assembly') == 'GRCh38'
                    and f.get('status') == 'released'
                    and f.get('file_format_type') in ('narrowPeak', 'broadPeak')):
                ot = f.get('output_type', '')
                rank = {'conservative IDR thresholded peaks': 0, 'replicated peaks': 1,
                        'IDR thresholded peaks': 2, 'pseudoreplicated peaks': 3,
                        'peaks': 4, 'hotspots': 5}.get(ot, 9)
                pd_ = 0 if f.get('preferred_default') else 1
                key = (pd_, rank, f['accession'])
                if best is None or key < best[0]: best = (key, f)
        if best is None:
            print("no hg38 peak file:", acc, term, assay); continue
        f = best[1]
        dest = os.path.join(CACHE, f"{f['accession']}.bed.gz")
        if not os.path.exists(dest):
            urllib.request.urlretrieve(B + f['href'], dest)
        n = sum(1 for _ in gzip.open(dest, 'rt'))
        rows.append(dict(biosample=term, assay=assay, experiment=acc, file=f['accession'],
                         output_type=f.get('output_type'), preferred_default=bool(f.get('preferred_default')),
                         n_peaks=n, summary=det.get('biosample_summary', '')[:80]))
        print(f"{term:32s} {assay:10s} {acc} {f['accession']} {f.get('output_type'):38s} n={n}")

json.dump(rows, open(os.path.join(CACHE, "manifest.json"), 'w'), indent=1)
print("total files", len(rows))
