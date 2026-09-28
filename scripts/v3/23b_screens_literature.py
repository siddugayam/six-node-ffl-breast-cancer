#!/usr/bin/env python3
"""PART 1D, second half) Literature search for published CRISPR screens in BREAST models
with migration / invasion / metastasis / ECM read-outs, since BioGRID ORCS 2.0.18 has none.
E-utilities, with the authors' institutional e-mail."""
import urllib.parse, urllib.request, json, time, re, sys
import pandas as pd
EMAIL="your.email@example.org"; TOOL="mirnaFFL_revision"
OUT="/path/to/revision/results/v3"
B="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def get(url):
    for _ in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r: return r.read().decode('utf-8','replace')
        except Exception as e:
            print("retry:", e); time.sleep(3)
    raise RuntimeError("failed "+url)

queries = {
 'breast_crispr_migration_invasion':
   '(CRISPR[tiab] AND (screen[tiab] OR screening[tiab] OR library[tiab])) AND (breast[tiab] OR mammary[tiab]) AND (migration[tiab] OR invasion[tiab] OR invasive[tiab] OR metastasis[tiab] OR metastatic[tiab])',
 'breast_crispr_ecm_collagen':
   '(CRISPR[tiab] AND (screen[tiab] OR screening[tiab])) AND (breast[tiab] OR mammary[tiab]) AND ("extracellular matrix"[tiab] OR collagen[tiab] OR fibroblast[tiab] OR stroma*[tiab])',
 'breast_crispr_metabolism':
   '(CRISPR[tiab] AND (screen[tiab] OR screening[tiab])) AND (breast[tiab] OR mammary[tiab]) AND (metabolism[tiab] OR metabolic[tiab])',
 'caf_crispr_screen':
   '(CRISPR[tiab] AND (screen[tiab] OR screening[tiab])) AND ("cancer-associated fibroblast"[tiab] OR "cancer associated fibroblasts"[tiab] OR CAF[tiab])',
}
rows=[]
for name,q in queries.items():
    u=f"{B}/esearch.fcgi?db=pubmed&retmax=200&sort=relevance&term={urllib.parse.quote(q)}&email={EMAIL}&tool={TOOL}&retmode=json"
    d=json.loads(get(u)); ids=d['esearchresult']['idlist']
    print(f"{name}: {d['esearchresult']['count']} records, fetching {len(ids)}")
    if not ids: continue
    for i in range(0,len(ids),100):
        chunk=ids[i:i+100]
        u2=f"{B}/esummary.fcgi?db=pubmed&id={','.join(chunk)}&email={EMAIL}&tool={TOOL}&retmode=json"
        s=json.loads(get(u2))['result']
        for pid in chunk:
            r=s.get(pid,{})
            rows.append(dict(query=name, pmid=pid, year=str(r.get('pubdate',''))[:4],
                             journal=r.get('source',''), title=r.get('title','')))
        time.sleep(0.4)
P=pd.DataFrame(rows).drop_duplicates(['query','pmid'])
P.to_csv(f"{OUT}/screens_literature_crispr_breast.csv", index=False)
print(f"\nwrote {len(P)} records to screens_literature_crispr_breast.csv")

# fetch abstracts for the migration/invasion set and keyword-triage
mig=P[P['query']=='breast_crispr_migration_invasion'].copy()
ids=mig.pmid.tolist()
abst={}
for i in range(0,len(ids),50):
    u=f"{B}/efetch.fcgi?db=pubmed&id={','.join(ids[i:i+50])}&rettype=abstract&retmode=text&email={EMAIL}&tool={TOOL}"
    txt=get(u); time.sleep(0.4)
    for block in txt.split("\n\n\n"):
        m=re.search(r'PMID:\s*(\d+)', block)
        if m: abst[m.group(1)]=block
print(f"abstracts retrieved for {len(abst)} of {len(ids)} migration/invasion records")
mig['abstract']=mig.pmid.map(lambda p: abst.get(p,''))
kw=re.compile(r'genome[- ]?wide|pooled|sgRNA librar|GeCKO|Brunello|transwell|matrigel|boyden|metastasis screen', re.I)
mig['looks_like_a_pooled_screen']=mig.abstract.str.contains(kw)
mig.to_csv(f"{OUT}/screens_literature_migration_abstracts.csv", index=False)
print("\ncandidate pooled CRISPR screens in breast with a migration/invasion/metastasis read-out:")
for _,r in mig[mig.looks_like_a_pooled_screen].iterrows():
    print(f"  PMID {r.pmid} ({r.year}, {r.journal}): {r.title[:120]}")
print(f"\n{int(mig.looks_like_a_pooled_screen.sum())} of {len(mig)} retrieved records match pooled-screen keywords")
