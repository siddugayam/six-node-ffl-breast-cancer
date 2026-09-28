#!/usr/bin/env python3
"""27b -- GTEx addendum: (a) independent cis-eQTLs for the miRNA genes themselves,
(b) an explicit, reproducible record of the breast / cultured-fibroblast NEGATIVES,
(c) the positions of every COL1A1/COL3A1 eQTL and fine-mapped variant relative to the TSS,
    ready for intersection with the motif scan and ReMap peaks.
"""
import json, time, urllib.request, urllib.parse
import numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"
BASE="https://gtexportal.org/api/v2"
def get(path,**p):
    url=f"{BASE}/{path}?"+urllib.parse.urlencode(p,doseq=True)
    for i in range(5):
        try:
            with urllib.request.urlopen(url,timeout=180) as r: return json.loads(r.read().decode())
        except Exception as e:
            if i==4: raise RuntimeError(f"{url}: {e}")
            time.sleep(4*(i+1))
def paged(path,**p):
    out=[];pg=0
    while True:
        d=get(path,page=pg,itemsPerPage=250,**p); out+=d.get('data',[])
        pi=d.get('paging_info',{})
        if pg+1>=pi.get('numberOfPages',1) or not d.get('data'): break
        pg+=1
    return out

g=pd.read_csv(f"{RES}/seqreg_gtex_genes.csv")
gmap=dict(zip(g.geneSymbol,g.gencodeId))
rows=[]
for s in ["MIR29A","MIR29B1","MIR29C"]:
    if s not in gmap: rows.append(dict(querySymbol=s,status="not in GTEx v8 gencode v26")); continue
    d=paged("association/independentEqtl", gencodeId=gmap[s])
    print(s, len(d))
    if not d: rows.append(dict(querySymbol=s, gencodeId=gmap[s], status="no independent cis-eQTL in any GTEx v8 tissue"))
    for r in d: r['querySymbol']=s; rows.append(r)
for s in ["MIR130AHG","MIR29B2CHG"]:
    rows.append(dict(querySymbol=s, status="gene not present in the GTEx v8 reference (gencode v26); no eQTL can be queried"))
pd.DataFrame(rows).to_csv(f"{RES}/seqreg_gtex_mirna_eqtl.csv", index=False)

ind=pd.read_csv(f"{RES}/seqreg_gtex_independent_eqtl.csv")
neg=[]
for s in ["COL1A1","COL3A1"]:
    have=set(ind[ind.querySymbol==s].tissueSiteDetailId)
    for t in ["Breast_Mammary_Tissue","Cells_Cultured_fibroblasts"]:
        neg.append(dict(gene=s, tissue=t, has_independent_cis_eQTL=(t in have),
                        n_tissues_with_signal=len(have),
                        tissues_with_signal=";".join(sorted(have))))
pd.DataFrame(neg).to_csv(f"{RES}/seqreg_gtex_breast_fibroblast_negatives.csv", index=False)
print(pd.DataFrame(neg)[["gene","tissue","has_independent_cis_eQTL","n_tissues_with_signal"]].to_string(index=False))

fm=pd.read_csv(f"{RES}/seqreg_gtex_finemapping.csv")
reg=pd.read_csv(f"{RES}/seqreg_regions.csv").set_index("region")
var=[]
for src,df in [("independentEqtl",ind),("fineMapping",fm)]:
    for _,r in df.iterrows():
        s=r.querySymbol
        if s not in ("COL1A1","COL3A1"): continue
        vid=r.variantId; ch,pos,ref,alt,_=vid.split("_")
        a=int(reg.loc[s,"anchor"]); st=reg.loc[s,"strand"]
        d=(int(pos)-a) if st=='+' else (a-int(pos))
        var.append(dict(source=src, gene=s, variantId=vid, chrom=ch, pos=int(pos), ref=ref, alt=alt,
                        tissue=r.get("tissueSiteDetailId"), pip=r.get("pip"), pValue=r.get("pValue"),
                        nes=r.get("nes"), dist_to_TSS=d))
v=pd.DataFrame(var).drop_duplicates(["source","gene","variantId","tissue"])
v.to_csv(f"{RES}/seqreg_gtex_collagen_variants.csv", index=False)
print("\ncollagen eQTL/fine-mapped variants:", len(v),
      " within 2 kb of TSS:", int((v.dist_to_TSS.abs()<=2000).sum()),
      " within 10 kb:", int((v.dist_to_TSS.abs()<=10000).sum()))
print(v.nsmallest(12,'dist_to_TSS',keep='all').head(12).to_string(index=False))
