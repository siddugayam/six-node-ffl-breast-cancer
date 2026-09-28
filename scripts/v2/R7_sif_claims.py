#!/usr/bin/env python3
"""Re-verify the deposited-SIF claims made in results/00_DATA_AUDIT.md sections A and B."""
import os, csv, collections, re
SIF='/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/SIF_files'
ATT='/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/node_attributes'
nets=['3-miR','3-TF','3-Comp','4-TF','5-TF','6-TF']
raw={}; allpairs=set(); alllabels=set()
for n in nets:
    L=[l.rstrip('\n').split('\t') for l in open(f'{SIF}/{n}.sif') if l.strip()]
    pairs=[(a,c) for a,b,c in L]
    raw[n]=(L,pairs)
    allpairs|=set(pairs); alllabels|={x for p in pairs for x in p}
    print(f"{n}.sif rows={len(L)} unique_pairs={len(set(pairs))} nodes={len({x for p in pairs for x in p})} selfloops={[p for p in set(pairs) if p[0]==p[1]]}")
print("union raw labels across SIFs:",len(alllabels))
# attribute-file types
typ=collections.defaultdict(set)
for n in nets:
    for r in csv.DictReader(open(f'{ATT}/{n}.csv')):
        typ[r['name']].add(r.get('Type',''))
print("labels in attribute files:",len(typ), "labels only in SIF:",len(alllabels-set(typ)),"only in attrs:",len(set(typ)-alllabels))
conflict=[k for k,v in typ.items() if len(v)>1]
print("labels with conflicting Type across files:",len(conflict),sorted(conflict))
# miRNA case duplication
mirs={l for l in alllabels if re.match(r'^hsa-(mir|miR|let)-',l)}
print("raw miRNA labels:",len(mirs))
low={l for l in mirs if l.startswith('hsa-mir-')}
up={l for l in mirs if not l.startswith('hsa-mir-')}
both={l[len('hsa-mir-'):] for l in low} & {l.split('-',2)[2] for l in up if l.startswith('hsa-miR-')}
print("names present in BOTH cases (raw label level):",len(both))
# upper-case node ever a target?
uptarget={t for (s,t) in allpairs if t.startswith('hsa-miR-')}
lowtarget={t for (s,t) in allpairs if t.startswith('hsa-mir-')}
print("distinct hsa-miR- (upper) labels appearing as a TARGET:",len(uptarget))
print("distinct hsa-mir- (lower) labels appearing as a TARGET:",len(lowtarget))
# miR-130a targets
def targets(lbl): return {t for (s,t) in allpairs if s==lbl}
t130=targets('hsa-miR-130a')|targets('hsa-mir-130a')
print("miR-130a pooled distinct targets:",len(t130),
      "| VEGFA in?",'VEGFA' in t130,"| MMP2 in?",'MMP2' in t130,"| COL1A1 in?",'COL1A1' in t130)
print("  regulators of VEGFA:",len({s for (s,t) in allpairs if t=='VEGFA'}),
      " MMP2:",len({s for (s,t) in allpairs if t=='MMP2'}),
      " COL1A1:",len({s for (s,t) in allpairs if t=='COL1A1'}))
print("HK2 present anywhere?", 'HK2' in alllabels)
# section 3.5 miRNAs
named=['miR-124','miR-133a','miR-133b','miR-143','miR-218','miR-3619','miR-412','miR-4776',
       'miR-4800','miR-516a','miR-516b','miR-7162','miR-92a-1']
canon={r['name'] for r in csv.DictReader(open('/path/to/revision/data/canonical_nodes.tsv'),delimiter='\t')}
col1a1_regs_raw={s for (s,t) in allpairs if t=='COL1A1'}
absent=[];present_noedge=[];present_edge=[]
for m in named:
    labs={l for l in alllabels if l.lower().endswith(m.lower()) and re.match(r'^hsa-(mir|miR)-',l)}
    incanon = any(l.replace('hsa-mir-','hsa-miR-') in canon for l in labs) or ('hsa-miR-'+m[4:] in canon)
    if not labs and not incanon: absent.append(m)
    elif labs & col1a1_regs_raw: present_edge.append(m)
    else: present_noedge.append((m,sorted(labs)))
print("\nsection-3.5 miRNAs, ABSENT from every SIF label:",len(absent),absent)
print("PRESENT with a COL1A1 edge:",len(present_edge),present_edge)
print("PRESENT but NO COL1A1 edge:",len(present_noedge),present_noedge)
# VEGFA degree in 3-Comp
L,pairs=raw['3-Comp']
deg=collections.Counter()
for s,t in set(pairs): deg[s]+=1; deg[t]+=1
print("\n3-Comp unique-pair degree: VEGFA",deg['VEGFA'],"TP53",deg['TP53'],"MYC",deg['MYC'],
      "RELA",deg['RELA'],"SP1",deg['SP1'],"CCND2",deg['CCND2'],
      "hsa-miR-21",deg['hsa-miR-21'],"hsa-mir-21",deg['hsa-mir-21'],
      "hsa-miR-34a",deg['hsa-miR-34a'],"hsa-mir-34a",deg['hsa-mir-34a'])
L2,p2=raw['3-TF']
deg2=collections.Counter()
for s,t in set(p2): deg2[s]+=1; deg2[t]+=1
print("3-TF unique-pair degree:   VEGFA",deg2['VEGFA'],"TP53",deg2['TP53'],"MYC",deg2['MYC'],
      "RELA",deg2['RELA'],"SP1",deg2['SP1'],"CCND2",deg2['CCND2'])
for n in ('3-TF','3-Comp'):
    d={r['name']:r.get('Degree') for r in csv.DictReader(open(f'{ATT}/{n}.csv'))}
    print(f"{n}.csv deposited Degree column: VEGFA={d.get('VEGFA')} TP53={d.get('TP53')} MYC={d.get('MYC')} RELA={d.get('RELA')} SP1={d.get('SP1')} CCND2={d.get('CCND2')} miR-130a={d.get('hsa-miR-130a')} mir-21={d.get('hsa-mir-21')} miR-124={d.get('hsa-miR-124')} mir-34a={d.get('hsa-mir-34a')}")
