#!/usr/bin/env python3
"""20_seqreg_regions.py -- define hg38 loci for sequence-level regulatory analysis
and fetch reference sequence from the UCSC REST API.

All coordinates were retrieved live from the UCSC API (ncbiRefSeqCurated, hg38)
by this script; nothing is hard-coded from memory except the region names.
Outputs: results/v3/seqreg_regions.csv, cached FASTA in cache/seqreg/seq/
"""
import json, os, sys, urllib.request, urllib.parse, time
import pandas as pd

ROOT = "/path/to/revision"
RES  = f"{ROOT}/results/v3"
CACHE= f"{ROOT}/cache/seqreg"
os.makedirs(f"{CACHE}/seq", exist_ok=True)

API = "https://api.genome.ucsc.edu"

def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i == tries-1: raise
            time.sleep(3*(i+1))

def refseq(chrom, start, end):
    u = f"{API}/getData/track?genome=hg38;track=ncbiRefSeqCurated;chrom={chrom};start={start};end={end}"
    return get(u).get('ncbiRefSeqCurated', [])

# --- 1. pull the annotation for the five loci -------------------------------
WIN = [("chr17",50170000,50215000),("chr2",188960000,189025000),
       ("chr7",130800000,130950000),("chr1",207740000,207890000),
       ("chr11",57600000,57680000)]
ann = []
for c,s,e in WIN:
    for r in refseq(c,s,e):
        ann.append(dict(chrom=r['chrom'], name=r['name'], symbol=r['name2'],
                        strand=r['strand'], txStart=r['txStart'], txEnd=r['txEnd'],
                        exonCount=r['exonCount']))
ann = pd.DataFrame(ann)
ann.to_csv(f"{RES}/seqreg_refseq_annotation.csv", index=False)

def tss(sym, acc=None):
    d = ann[ann.symbol==sym] if acc is None else ann[ann.name==acc]
    if len(d)==0: raise SystemExit(f"no annotation for {sym}/{acc}")
    r = d.iloc[0]
    return r.chrom, (int(r.txEnd) if r.strand=='-' else int(r.txStart)), r.strand, int(r.txStart), int(r.txEnd)

rows = []
def add(name, chrom, pos, strand, kind, note, span=None):
    rows.append(dict(region=name, chrom=chrom, anchor=pos, strand=strand, kind=kind,
                     span_start=(span[0] if span else pos), span_end=(span[1] if span else pos),
                     note=note))

for sym, acc, label, kind in [
    ("COL1A1","NM_000088.4","COL1A1","protein_coding_TSS"),
    ("COL3A1","NM_000090.4","COL3A1","protein_coding_TSS"),
    ("MIR29A","NR_029503.1","MIR29A","miRNA_hairpin"),
    ("MIR29B1","NR_029517.1","MIR29B1","miRNA_hairpin"),
    ("MIR29C","NR_029832.1","MIR29C","miRNA_hairpin"),
    ("MIR29B2","NR_029518.1","MIR29B2","miRNA_hairpin"),
    ("MIR130A","NR_029673.1","MIR130A","miRNA_hairpin"),
    ("MIR29B2CHG","NR_135298.1","MIR29B2CHG","host_gene_TSS"),
    ("MIR130AHG","NR_186232.1","MIR130AHG","host_gene_TSS"),
    ("LINC-PINT","NR_110473.1","LINC_PINT_prox","host_gene_TSS"),
]:
    c,p,st,a,b = tss(sym, acc)
    add(label, c, p, st, kind, f"{acc} {sym}", span=(a,b))

reg = pd.DataFrame(rows)
reg.to_csv(f"{RES}/seqreg_regions.csv", index=False)
print(reg.to_string(index=False))

# --- 2. fetch sequence windows ---------------------------------------------
# Enformer needs 196,608 bp; Borzoi 524,288 bp. Fetch 600 kb centred on each anchor.
HALF = 300_000
def fetch_seq(chrom, start, end):
    key = f"{CACHE}/seq/{chrom}_{start}_{end}.txt"
    if os.path.exists(key) and os.path.getsize(key) > (end-start)*0.9:
        return open(key).read().strip()
    u = f"{API}/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    d = get(u)
    s = d['dna'].upper()
    assert len(s) == end-start, (len(s), end-start)
    open(key,'w').write(s)
    return s

meta=[]
for _,r in reg.iterrows():
    st = max(0, int(r.anchor)-HALF); en = int(r.anchor)+HALF
    s = fetch_seq(r.chrom, st, en)
    gc = (s.count('G')+s.count('C'))/max(1,len(s)-s.count('N'))
    meta.append(dict(region=r.region, chrom=r.chrom, win_start=st, win_end=en,
                     length=len(s), n_frac=s.count('N')/len(s), gc=round(gc,4)))
    print(r.region, r.chrom, st, en, len(s), round(gc,4))
pd.DataFrame(meta).to_csv(f"{RES}/seqreg_sequence_windows.csv", index=False)
print("done")
