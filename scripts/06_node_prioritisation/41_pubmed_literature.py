#!/usr/bin/env python3
"""Literature volume + candidate references for the 30 prioritised nodes.
Uses NCBI E-utilities. Every PMID returned here is subsequently verified by esummary,
so no citation is ever asserted without a retrieved record."""
import json, time, urllib.parse, urllib.request, os, sys, csv

EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = "/path/to/revision/results/v5"
TOOL = "&tool=ffl_revision&email=your.email@example.org"

def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i == tries - 1: raise
            time.sleep(2 + 2*i)

def esearch(term, retmax=0, sort="relevance"):
    u = (EUT + "esearch.fcgi?db=pubmed&retmode=json&retmax=%d&sort=%s&term=%s%s"
         % (retmax, sort, urllib.parse.quote(term), TOOL))
    time.sleep(0.4)
    return get(u)["esearchresult"]

def esummary(pmids):
    if not pmids: return {}
    u = (EUT + "esummary.fcgi?db=pubmed&retmode=json&id=%s%s" % (",".join(pmids), TOOL))
    time.sleep(0.4)
    return get(u)["result"]

NODES = {
 # name : (pubmed search term for the entity, human-readable)
 "CCND2": '("CCND2"[tiab] OR "cyclin D2"[tiab])',
 "COL1A1": '("COL1A1"[tiab] OR "collagen type I alpha 1"[tiab])',
 "STAT5A": '("STAT5A"[tiab] OR "STAT5"[tiab])',
 "MYBL2": '("MYBL2"[tiab] OR "B-Myb"[tiab])',
 "FN1": '("FN1"[tiab] OR "fibronectin"[tiab])',
 "PDGFRB": '("PDGFRB"[tiab] OR "PDGFR-beta"[tiab] OR "PDGFR beta"[tiab])',
 "MET": '("MET proto-oncogene"[tiab] OR "c-Met"[tiab] OR "HGF receptor"[tiab])',
 "CXCL12": '("CXCL12"[tiab] OR "SDF-1"[tiab] OR "stromal cell-derived factor 1"[tiab])',
 "MMP14": '("MMP14"[tiab] OR "MT1-MMP"[tiab] OR "membrane type 1 matrix metalloproteinase"[tiab])',
 "PLAU": '("PLAU"[tiab] OR "urokinase plasminogen activator"[tiab] OR "uPA"[tiab])',
 "E2F1": '"E2F1"[tiab]',
 "EZH2": '("EZH2"[tiab] OR "enhancer of zeste homolog 2"[tiab])',
 "GATA3": '"GATA3"[tiab]',
 "BRCA1": '"BRCA1"[tiab]',
 "JUN": '("c-Jun"[tiab] OR "JUN proto-oncogene"[tiab] OR "AP-1"[tiab])',
 "EGR2": '("EGR2"[tiab] OR "early growth response 2"[tiab] OR "Krox20"[tiab])',
 "ESR1": '("ESR1"[tiab] OR "estrogen receptor alpha"[tiab] OR "oestrogen receptor alpha"[tiab])',
 "SREBF1": '("SREBF1"[tiab] OR "SREBP-1"[tiab] OR "SREBP1"[tiab])',
 "DNMT1": '("DNMT1"[tiab] OR "DNA methyltransferase 1"[tiab])',
 "E2F3": '"E2F3"[tiab]',
 "hsa-miR-21": '("miR-21"[tiab] OR "microRNA-21"[tiab] OR "miRNA-21"[tiab])',
 "hsa-miR-195": '("miR-195"[tiab] OR "microRNA-195"[tiab])',
 "hsa-miR-204": '("miR-204"[tiab] OR "microRNA-204"[tiab])',
 "hsa-miR-383": '("miR-383"[tiab] OR "microRNA-383"[tiab])',
 "hsa-miR-124": '("miR-124"[tiab] OR "microRNA-124"[tiab])',
 "hsa-miR-155": '("miR-155"[tiab] OR "microRNA-155"[tiab])',
 "hsa-miR-429": '("miR-429"[tiab] OR "microRNA-429"[tiab])',
 "hsa-miR-141": '("miR-141"[tiab] OR "microRNA-141"[tiab])',
 "hsa-miR-34a": '("miR-34a"[tiab] OR "microRNA-34a"[tiab])',
 "hsa-miR-101": '("miR-101"[tiab] OR "microRNA-101"[tiab])',
}
BC = '"Breast Neoplasms"[Mesh]'

rows = []
cands = {}
for name, term in NODES.items():
    tot = int(esearch(term)["count"])
    bc  = esearch(term + " AND " + BC)
    nbc = int(bc["count"])
    # candidate references: most relevant, and most relevant reviews
    rel = esearch(term + " AND " + BC, retmax=12, sort="relevance")["idlist"]
    rev = esearch(term + " AND " + BC + ' AND (review[pt] OR "systematic review"[pt])',
                  retmax=6, sort="relevance")["idlist"]
    ids = list(dict.fromkeys(rel + rev))
    summ = esummary(ids)
    c = []
    for pid in ids:
        s = summ.get(pid)
        if not isinstance(s, dict): continue
        c.append(dict(pmid=pid, year=s.get("pubdate", "")[:4], journal=s.get("source", ""),
                      title=s.get("title", ""),
                      pubtype=";".join(s.get("pubtype", [])),
                      first=(s.get("authors") or [{}])[0].get("name", "")))
    cands[name] = c
    rows.append(dict(node=name, pubmed_total=tot, pubmed_breast=nbc,
                     frac_breast=(nbc / tot if tot else 0)))
    print(f"{name:14s} total={tot:7d}  breast={nbc:6d}  candidates={len(c)}", flush=True)

with open(os.path.join(OUT, "node_literature_volume.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["node", "pubmed_total", "pubmed_breast", "frac_breast"])
    w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT, "node_literature_candidates.json"), "w") as f:
    json.dump(cands, f, indent=1)
print("done")
