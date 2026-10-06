#!/usr/bin/env python3
"""Fetch NCI PDMR (cBioPortal study pancan_pdmr_2025) RNA-seq V2 RSEM expression
for breast-cancer samples, for a reference gene panel + 40 target genes.
Writes a long TSV: sampleId, entrezGeneId, value."""
import json, urllib.request, time, sys, os

OUT = "/path/to/revision/cache/v6/pdx"
os.makedirs(OUT, exist_ok=True)
PROF = "pancan_pdmr_2025_rna_seq_v2_mrna"
URL = f"https://www.cbioportal.org/api/molecular-profiles/{PROF}/molecular-data/fetch?projection=ID"

TARGETS = {"COL1A1":1277,"COL3A1":1281,"FN1":2335,"PDGFRB":5159,"CXCL12":6387,"POSTN":10631,
 "E2F1":1869,"EZH2":2146,"MYBL2":4605,"CCND2":894,"GATA3":2625,"BRCA1":672,"JUN":3725,"EGR2":1959,
 "ESR1":2099,"SREBF1":6720,"DNMT1":1786,"E2F3":1871,"STAT5A":6776,"MET":4233,"MMP14":4323,"PLAU":5328,
 "NFKB1":4790,"ETS1":2113,"EPCAM":4072,"KRT8":3856,"KRT18":3875,"PTPRC":5788,"PECAM1":5175,"ACTA2":59,
 "DCN":1634,"LUM":4060,"VIM":7431,"MKI67":4288,"COL1A2":1278,"THY1":7070,"FAP":2191,"PDGFRA":5156,
 "ERBB2":2064,"KRT19":3880}

ref = json.load(open('/path/to/scratch/ref_panel_entrez.json'))
br  = json.load(open('/path/to/scratch/pdmr_breast_meta.json'))
samples = sorted(br.keys())

genes = dict(ref); genes.update(TARGETS)
json.dump(genes, open(os.path.join(OUT,'gene_symbol_entrez.json'),'w'))
ids = sorted(set(genes.values()))
print(f"{len(ids)} genes x {len(samples)} samples", flush=True)

path = os.path.join(OUT,'pdmr_breast_expr_long.tsv')
seen = 0
with open(path,'w') as fh:
    fh.write("sampleId\tentrezGeneId\tvalue\n")
    for i in range(0, len(ids), 150):
        ch = ids[i:i+150]
        body = json.dumps({"entrezGeneIds":ch, "sampleIds":samples}).encode()
        for attempt in range(4):
            try:
                req = urllib.request.Request(URL, data=body,
                        headers={"Content-Type":"application/json"}, method="POST")
                d = json.load(urllib.request.urlopen(req, timeout=600))
                break
            except Exception as ex:
                print("retry", i, ex, flush=True); time.sleep(10)
        else:
            print("FAILED chunk", i, flush=True); continue
        for r in d:
            v = r.get('value')
            if v is None: continue
            fh.write(f"{r['sampleId']}\t{r['entrezGeneId']}\t{v}\n"); seen += 1
        print(f"chunk {i} done, cumulative {seen}", flush=True)
        time.sleep(0.5)
print("TOTAL", seen)
