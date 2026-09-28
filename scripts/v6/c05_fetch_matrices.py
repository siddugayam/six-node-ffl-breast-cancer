#!/usr/bin/env python3
"""Download GEO series matrices + platform annotations for the usable miRNA cohorts."""
import os, sys, time, urllib.request, urllib.error
D="/path/to/revision/cache/v6/cohorts/matrix"
P="/path/to/revision/cache/v6/cohorts/gpl"
os.makedirs(D,exist_ok=True); os.makedirs(P,exist_ok=True)

def dl(url, path, tries=4):
    if os.path.exists(path) and os.path.getsize(path)>1000:
        print("cached", os.path.basename(path), os.path.getsize(path)); return True
    for i in range(tries):
        try:
            req=urllib.request.Request(url, headers={"User-Agent":"mirna-ffl/1.0"})
            with urllib.request.urlopen(req, timeout=600) as r, open(path,"wb") as f:
                while True:
                    b=r.read(1<<20)
                    if not b: break
                    f.write(b)
            print("OK", url, os.path.getsize(path), flush=True); return True
        except urllib.error.HTTPError as e:
            print("HTTP %s %s"%(e.code,url), flush=True); return False
        except Exception as e:
            print("retry", e, flush=True); time.sleep(4+4*i)
    return False

def series_matrix(gse, gpl=None):
    stub = gse[:-3]+"nnn"
    fn = "%s_series_matrix.txt.gz"%gse if gpl is None else "%s-%s_series_matrix.txt.gz"%(gse,gpl)
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s"%(stub,gse,fn)
    return dl(url, os.path.join(D,fn))

def gpl_annot(gpl):
    stub = gpl[:-3]+"nnn" if len(gpl)>6 else "GPLnnn"
    for fn,url in [("%s.annot.gz"%gpl,"https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz"%(stub,gpl,gpl)),
                   ("%s_family.soft.gz"%gpl,"https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/soft/%s_family.soft.gz"%(stub,gpl,gpl))]:
        if dl(url, os.path.join(P,fn)): return True
    return False

TARGETS = [("GSE22216",None), ("GSE59829",None), ("GSE78870",None),
           ("GSE37405","GPL13703"), ("GSE37405","GPL14149"), ("GSE37405","GPL15462"),
           ("GSE73002",None), ("GSE57897",None), ("GSE103161",None), ("GSE126125",None),
           ("GSE40267",None), ("GSE26666","GPL8227"), ("GSE97811",None), ("GSE28321",None)]
for gse,gpl in TARGETS:
    series_matrix(gse,gpl)
for gpl in ["GPL8178","GPL8179","GPL8227","GPL20662","GPL13703","GPL14149","GPL15462",
            "GPL18941","GPL18722","GPL23960","GPL10850","GPL21263","GPL5106"]:
    gpl_annot(gpl)
