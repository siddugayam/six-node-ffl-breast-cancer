#!/usr/bin/env python3
"""27_seqreg_gtex.py -- GTEx v8/v10 cis-eQTL and sQTL evidence for COL1A1, COL3A1 and the
miR-29 host genes, with emphasis on breast and cultured fibroblasts.
Endpoints used (GTEx portal API v2): /reference/gene, /association/independentEqtl,
/association/fineMapping, /association/singleTissueSqtl, /association/singleTissueEqtlByLocation.
"""
import json, os, time, urllib.request, urllib.parse
import pandas as pd

ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"
BASE="https://gtexportal.org/api/v2"

def get(path, **params):
    q=urllib.parse.urlencode(params, doseq=True)
    url=f"{BASE}/{path}?{q}"
    for i in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i==4: raise RuntimeError(f"{url}: {e}")
            time.sleep(4*(i+1))

def paged(path, **params):
    out=[]; page=0
    while True:
        d=get(path, page=page, itemsPerPage=250, **params)
        out += d.get('data', [])
        pi=d.get('paging_info',{})
        if page+1 >= pi.get('numberOfPages',1) or not d.get('data'): break
        page+=1
    return out

GENES=["COL1A1","COL3A1","MIR29B2CHG","MIR29A","MIR29B1","MIR29C","LINC-PINT",
       "NFKB1","RELA","SP1","ETS1"]
g=get("reference/gene", geneId=GENES)
gdf=pd.DataFrame(g['data'])
gdf.to_csv(f"{RES}/seqreg_gtex_genes.csv", index=False)
print(gdf[['geneSymbol','gencodeId','chromosome','tss','strand']].to_string(index=False))
gmap={r.geneSymbol:r.gencodeId for r in gdf.itertuples()}

FOCUS=["COL1A1","COL3A1","MIR29B2CHG","LINC-PINT"]
FOCUS=[x for x in FOCUS if x in gmap]

rows=[]
for s in FOCUS:
    d=paged("association/independentEqtl", gencodeId=gmap[s])
    for r in d: r['querySymbol']=s
    rows+=d
    print("independentEqtl", s, len(d), flush=True)
ind=pd.DataFrame(rows)
ind.to_csv(f"{RES}/seqreg_gtex_independent_eqtl.csv", index=False)

fm=[]
for s in FOCUS:
    try:
        d=paged("association/fineMapping", gencodeId=gmap[s])
    except Exception as e:
        print("fineMapping fail", s, e); d=[]
    for r in d: r['querySymbol']=s
    fm+=d
    print("fineMapping", s, len(d), flush=True)
pd.DataFrame(fm).to_csv(f"{RES}/seqreg_gtex_finemapping.csv", index=False)

sq=[]
for s in FOCUS:
    for t in ["Breast_Mammary_Tissue","Cells_Cultured_fibroblasts","Artery_Tibial",
              "Skin_Sun_Exposed_Lower_leg","Adipose_Subcutaneous"]:
        try:
            d=paged("association/singleTissueSqtl", gencodeId=gmap[s], tissueSiteDetailId=t)
        except Exception as e:
            print("sqtl fail", s, t, e); d=[]
        for r in d: r['querySymbol']=s
        sq+=d
        print("sQTL", s, t, len(d), flush=True)
pd.DataFrame(sq).to_csv(f"{RES}/seqreg_gtex_sqtl.csv", index=False)

# eQTLs by location: promoter windows (+/-5 kb of TSS) in fibroblast/breast tissues
reg=pd.read_csv(f"{RES}/seqreg_regions.csv")
loc=[]
for _,r in reg[reg.region.isin(["COL1A1","COL3A1","MIR29B2CHG"])].iterrows():
    for t in ["Breast_Mammary_Tissue","Cells_Cultured_fibroblasts"]:
        try:
            d=get("association/singleTissueEqtlByLocation", tissueSiteDetailId=t,
                  chromosome=r.chrom, start=int(r.anchor)-5000, end=int(r.anchor)+5000)
        except Exception as e:
            print("byLoc fail", r.region, t, e); continue
        for x in d.get('singleTissueEqtl',[]):
            x['queryRegion']=r.region; x['queryTissue']=t; loc.append(x)
        print("byLocation", r.region, t, len(d.get('singleTissueEqtl',[])), flush=True)
pd.DataFrame(loc).to_csv(f"{RES}/seqreg_gtex_eqtl_by_location.csv", index=False)
print("done")
