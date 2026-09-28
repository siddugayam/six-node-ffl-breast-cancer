#!/usr/bin/env python3
"""08c -- why the Composite-FFL result depends entirely on the null model.

The motif analysis (script 08) shows Composite-FFL at Z=+28.7 under a null that destroys
reciprocity (NM1) but only Z=+2.8 under a null that preserves it (NM2).  Everything therefore
turns on ONE number: the fraction of TF->miRNA edges that carry a reciprocal miRNA->TF edge.
This script measures that fraction in the deposited network and in the same node set rebuilt
from the complete primary databases, and tests the difference.
Writes results/motif_reciprocity_audit.csv."""
import csv, collections
from scipy.stats import fisher_exact
REV = '/path/to/revision'
D = f'{REV}/data/control'
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
nodes = {r['name']: r['type'] for r in rd(f'{REV}/data/canonical_nodes.tsv')}
tfs = {v for v, t in nodes.items() if t == 'TF'}
gns = {v for v, t in nodes.items() if t == 'Gene'}
mis = {v for v, t in nodes.items() if t == 'miRNA'}
prot = tfs | gns
TRR = collections.defaultdict(set)
for r in rd(f'{D}/trrust_pairs.tsv'): TRR[r['TF']].add(r['target'])
TMR = collections.defaultdict(set)
for r in rd(f'{D}/transmir_pairs.tsv'): TMR[r['TF']].add(r['mirna'])
VAL = collections.defaultdict(set)
for r in rd(f'{D}/validated_mirna_target_pairs.tsv'): VAL[r['mirna']].add(r['gene'])

REB = set()
for t in sorted(tfs):
    for g in sorted(TRR.get(t, ())):
        if g in prot and g != t: REB.add((t, g))
    for mi in sorted(TMR.get(t, ())):
        if mi in mis: REB.add((t, mi))
for mi in sorted(mis):
    for g in sorted(VAL.get(mi, ())):
        if g in prot: REB.add((mi, g))
PUB = set()
for r in rd(f'{REV}/data/canonical_edges.tsv'):
    if r['edge_type'] in ('miRNA_target', 'TF_miRNA', 'TF_target'):
        PUB.add((r['source'], r['target']))

rows = []
def audit(name, E):
    a = [(x, y) for x, y in E if nodes[x] == 'TF' and nodes[y] == 'miRNA']
    b = [(x, y) for x, y in E if nodes[x] == 'miRNA' and nodes[y] == 'TF']
    ra = sum(1 for x, y in a if (y, x) in E)
    rb = sum(1 for x, y in b if (y, x) in E)
    rows.append(dict(network=name, n_edges=len(E), n_TF_to_miRNA=len(a),
                     n_TF_to_miRNA_reciprocated=ra,
                     pct_TF_to_miRNA_reciprocated=round(100 * ra / len(a), 2) if a else 0,
                     n_miRNA_to_TF=len(b), n_miRNA_to_TF_reciprocated=rb,
                     pct_miRNA_to_TF_reciprocated=round(100 * rb / len(b), 2) if b else 0,
                     n_mutual_pairs=ra))
    print(f'{name:<16} |E|={len(E):<7} TF->miR {len(a):<6} recip {ra:<5} '
          f'({100*ra/len(a) if a else 0:.1f}%)  miR->TF {len(b):<6} recip {rb:<5} '
          f'({100*rb/len(b) if b else 0:.1f}%)')
audit('BRCA_published', PUB)
audit('BRCA_rebuilt', REB)

# Is the published miRNA->TF layer a reciprocity-enriched subset of the available one?
pub_mt = [(x, y) for x, y in PUB if nodes[x] == 'miRNA' and nodes[y] == 'TF']
reb_mt = [(x, y) for x, y in REB if nodes[x] == 'miRNA' and nodes[y] == 'TF']
kept = [e for e in reb_mt if e in PUB]
drop = [e for e in reb_mt if e not in PUB]
kr = sum(1 for x, y in kept if (y, x) in REB)
dr = sum(1 for x, y in drop if (y, x) in REB)
odds, p = fisher_exact([[kr, len(kept) - kr], [dr, len(drop) - dr]], alternative='greater')
print(f'\nOf the {len(reb_mt)} VALIDATED miRNA->TF edges available inside the BRCA node set, '
      f'the published network keeps {len(kept)} and drops {len(drop)}.')
print(f'   kept   : {kr}/{len(kept)} ({100*kr/len(kept):.1f}%) have a reciprocal TransmiR TF->miRNA edge')
print(f'   dropped: {dr}/{len(drop)} ({100*dr/len(drop):.1f}%) have one')
print(f'   Fisher exact one-sided odds ratio = {odds:.1f}, p = {p:.3g}')
rows.append(dict(network='published_vs_available_miRNA_to_TF_selection',
                 n_edges=len(reb_mt), n_TF_to_miRNA=len(kept),
                 n_TF_to_miRNA_reciprocated=kr,
                 pct_TF_to_miRNA_reciprocated=round(100 * kr / len(kept), 2),
                 n_miRNA_to_TF=len(drop), n_miRNA_to_TF_reciprocated=dr,
                 pct_miRNA_to_TF_reciprocated=round(100 * dr / len(drop), 2),
                 n_mutual_pairs=f'fisher_OR={odds:.2f};p={p:.3g}'))
with open(f'{REV}/results/motif_reciprocity_audit.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow(r)
print('\nwrote results/motif_reciprocity_audit.csv')
