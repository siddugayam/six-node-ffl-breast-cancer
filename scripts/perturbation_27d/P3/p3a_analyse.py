#!/usr/bin/env python3
"""P3a: validate the analysed network's miRNA->target edges by evidence tier with miRNA gain/loss datasets (SETTINGS.md, P3).
Datasets: the P3a contrasts (P3a/de), the human P1b contrasts (P1b/de, tag primary) and the human P2 contrasts (P2/de,
cls starting 'human'); if more than 40, the 40 with the most samples (ties: GSE number) are kept and the rest listed.
Per dataset: members = the perturbed network miRNA node(s); effect = sign-aligned repression (-log2FC gain, +log2FC loss)
over expressed genes; for each tier (data/edge_evidence_tier.tsv: strong, weak, predicted_only) the members' network targets
in that tier; background = expressed genes with no network edge from any member and no TargetScan 8.0 site of any type
for the members' seeds; shift = median effect of the tier's targets minus the background median; SE from 1,000 bootstrap
resamples (per-dataset seed, _seed.py).  Pooling: REML per tier; decision per tier: SUPPORTED if the pooled shift > 0
with 95 % CI excluding 0, NOT SUPPORTED otherwise.  Predicted minus strong: per-dataset difference, SE = sqrt of the
summed variances, pooled by REML.
Writes p3a_per_dataset.tsv, p3a_pooled.tsv, p3a_datasets_used.tsv and p3a_summary.txt (numbered lines)."""
import os, sys, glob, json, re, io, zipfile, subprocess
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(INB); sys.path.insert(0, INB)
import _dl
from _seed import rng_for
D = f'{INB}/P3'
T = pd.read_csv(f'{REV}/data/edge_evidence_tier.tsv', sep='\t')
T = T[T.edge_type == 'miRNA_target']
TIER = {t: g.groupby('source').target.apply(set).to_dict() for t, g in T.groupby('tier')}
NET = T.groupby('source').target.apply(set).to_dict()
TS = pd.read_csv(os.path.join(_dl.CACHE, 'common/targetscan/ts80_human_family_sites.tsv.gz'), sep='\t')
ANY = {s: set(g[(g.cons_total + g.noncons_total) > 0].gene) for s, g in TS.groupby('family')}
FI = pd.read_csv(io.BytesIO(zipfile.ZipFile(os.path.join(_dl.CACHE, 'common/targetscan/miR_Family_Info.txt.zip')).read('miR_Family_Info.txt')), sep='\t')
FI = FI[FI['Species ID'] == 9606]
def seeds(node):
    rx = re.compile(r'^hsa-' + re.escape(node.replace('hsa-', '')) + r'(-[35]p)?$')
    return sorted(set(FI[FI['MiRBase ID'].str.match(rx)]['Seed+m8']))
MIR29 = 'hsa-miR-29a;hsa-miR-29b;hsa-miR-29c'
def members_of(m):
    if m.get('mirna', '').startswith('hsa-'): return m['mirna'].split(';')
    x = m.get('mirna', '') or m.get('tf', '')
    if x.startswith('miR-29'):                                   # P1b notation: miR-29, miR-29a, miR-29b, miR-29c, miR-29a/b/c
        tail = x.replace('miR-29', '')
        return MIR29.split(';') if tail in ('', 'a/b/c') else [f'hsa-miR-29{c}' for c in tail.split('/') if c]
    if x.startswith('miR-101'): return ['hsa-miR-101']
    return []
cands = []
for mp in sorted(glob.glob(f'{INB}/P3/de/*.meta.json') + glob.glob(f'{INB}/P1b/de/*.meta.json') + glob.glob(f'{INB}/P2/de/*.meta.json')):
    m = json.load(open(mp))
    if m.get('tag', 'primary') != 'primary' or m.get('cls') == 'mouse': continue
    mem = [x for x in members_of(m) if x in NET]
    if not mem: continue
    n = (m.get('n_pert') or 0) + (m.get('n_ctrl') or 0)
    cands.append(dict(meta=mp, gse=m['gse'], cell=m['cell'], members=';'.join(mem), direction=m.get('direction', ''), n_samples=int(n), source=mp.split('/')[-3]))
C = pd.DataFrame(cands).drop_duplicates(['gse', 'cell', 'members'])
C['gse_num'] = C.gse.str[3:].astype(int)
C = C.sort_values(['n_samples', 'gse_num'], ascending=[False, True])
C['used'] = [i < 40 for i in range(len(C))]
C.drop(columns=['meta']).to_csv(f'{D}/p3a_datasets_used.tsv', sep='\t', index=False)
rows = []
for c in C[C.used].itertuples():
    m = json.load(open(c.meta)); t = pd.read_csv(c.meta.replace('.meta.json', '.tsv'), sep='\t').set_index('gene')
    sgn = -1.0 if c.direction == 'gain' else 1.0
    eff = sgn * t[t.expressed.astype(str) == 'True'].logFC.dropna()
    mem = c.members.split(';')
    anyedge = set().union(*[NET.get(x, set()) for x in mem]); anysite = set().union(*[ANY.get(s, set()) for x in mem for s in seeds(x)])
    bg = eff[~eff.index.isin(anyedge) & ~eff.index.isin(anysite)]
    for tier in ('strong', 'weak', 'predicted_only'):
        tg = set().union(*[TIER.get(tier, {}).get(x, set()) for x in mem])
        a = eff[eff.index.isin(tg)]
        r = dict(gse=c.gse, cell=c.cell, members=c.members, direction=c.direction, tier=tier, n_targets=len(a), n_background=len(bg))
        if len(a) >= 3 and len(bg) >= 3:
            rng = rng_for(f'{c.gse}|{c.cell}|{tier}')
            r['shift'] = float(a.median() - bg.median())
            r['se'] = float(np.std([np.median(rng.choice(a.values, len(a))) - np.median(rng.choice(bg.values, len(bg))) for _ in range(1000)], ddof=1))
        rows.append(r)
R = pd.DataFrame(rows); R.to_csv(f'{D}/p3a_per_dataset.tsv', sep='\t', index=False)
pin = [dict(group=x.tier, yi=x.shift, sei=x.se) for x in R.dropna(subset=['shift']).itertuples()]
w = R.pivot_table(index=['gse', 'cell'], columns='tier', values=['shift', 'se'])
for (g, cl), x in w.iterrows():
    if np.isfinite(x.get(('shift', 'predicted_only'), np.nan)) and np.isfinite(x.get(('shift', 'strong'), np.nan)):
        pin.append(dict(group='predicted minus strong', yi=x[('shift', 'predicted_only')] - x[('shift', 'strong')], sei=float(np.hypot(x[('se', 'predicted_only')], x[('se', 'strong')]))))
pd.DataFrame(pin).to_csv(f'{D}/p3a_pool_in.tsv', sep='\t', index=False)
subprocess.run(['Rscript', f'{INB}/_rma.R', f'{D}/p3a_pool_in.tsv', f'{D}/p3a_pooled.tsv'], check=True)
PT = pd.read_csv(f'{D}/p3a_pooled.tsv', sep='\t')
L = []; P = L.append
P(f'candidate datasets (human, perturbing network miRNAs): {len(C)}; used (40 with the most samples): {int(C.used.sum())}; listed but not used: {int((~C.used).sum())}')
for tier in ('strong', 'weak', 'predicted_only', 'predicted minus strong'):
    x = PT[PT.group == tier]
    if not len(x) or x.iloc[0].k == 0: P(f'{tier}: no dataset with >= 3 targets'); continue
    x = x.iloc[0]
    dec = '' if tier == 'predicted minus strong' else (' -> SUPPORTED' if x.ci_lb > 0 else ' -> NOT SUPPORTED')
    P(f'{tier}: k = {int(x.k)}; pooled shift {x.est:+.4f} (95 % CI {x.ci_lb:+.4f} to {x.ci_ub:+.4f}), p {x.p:.3g}; tau2 {x.tau2:.3g}{dec}')
open(f'{D}/p3a_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
