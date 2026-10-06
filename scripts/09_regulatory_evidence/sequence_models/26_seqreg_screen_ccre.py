#!/usr/bin/env python3
"""26_seqreg_screen_ccre.py -- ENCODE SCREEN candidate cis-regulatory elements (cCREs)
and UCSC CpG islands at the COL1A1/COL3A1/miR-29 loci (hg38).
Source: UCSC REST API tracks encodeCcreCombined (ENCODE Registry of cCREs) and cpgIslandExt.
"""
import json, os, time, urllib.request
import pandas as pd

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"
API="https://api.genome.ucsc.edu"

def get(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=240) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i==tries-1: raise
            time.sleep(4*(i+1))

reg = pd.read_csv(f"{RES}/seqreg_regions.csv")
HALF = 100_000
rows_c, rows_cpg = [], []
for _,r in reg.iterrows():
    s=max(0,int(r.anchor)-HALF); e=int(r.anchor)+HALF
    for track,store in [("encodeCcreCombined",rows_c),("cpgIslandExt",rows_cpg)]:
        d=get(f"{API}/getData/track?genome=hg38;track={track};chrom={r.chrom};start={s};end={e}")
        for it in d.get(track,[]):
            base=dict(region=r.region, anchor=int(r.anchor), strand=r.strand,
                      chrom=it['chrom'], start=it['chromStart'], end=it['chromEnd'])
            mid=(it['chromStart']+it['chromEnd'])//2
            base['dist_to_anchor']=mid-int(r.anchor)
            if r.strand=='-': base['dist_to_anchor']=-base['dist_to_anchor']
            if track=="encodeCcreCombined":
                base.update(accession=it['name'], ccre_class=it['ccre'],
                            encode_label=it['encodeLabel'], ucsc_label=it['ucscLabel'],
                            zscore=it['zScore'], score=it['score'])
            else:
                base.update(name=it['name'], length=it['length'], cpgNum=it['cpgNum'],
                            gcNum=it['gcNum'], perCpg=it['perCpg'], perGc=it['perGc'],
                            obsExp=it['obsExp'])
            store.append(base)
    print(r.region, "cCRE so far", len(rows_c), "CpG", len(rows_cpg), flush=True)

pd.DataFrame(rows_c).to_csv(f"{RES}/seqreg_encode_ccre.csv", index=False)
pd.DataFrame(rows_cpg).to_csv(f"{RES}/seqreg_cpg_islands.csv", index=False)
print("cCREs", len(rows_c), "CpG islands", len(rows_cpg))
