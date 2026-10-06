#!/usr/bin/env python3
"""P2: do co-transcribed miRNAs regulate their shared targets convergently?  (SETTINGS.md P2 and change 8.)
Per dataset (P2/de/*.tsv, one contrast perturbing >= 2 members of a cluster):
  effect  = sign-aligned repression (-log2FC for gain, +log2FC for loss), expressed genes only;
  classes = number of perturbed members targeting the gene: 0, exactly 1, >= 2, under
            (a) network edges (data/canonical_edges.tsv, miRNA_target, member node -> gene), and
            (b) TargetScan 8.0 (>= 1 conserved 8mer or 7mer-m8 site for a seed of the member; seeds = the member's arms
                whose TargetScan family is conserved (Family Conservation >= 1), or all arms if none);
  six-node class = G1 of every BHAT6 composite instance (analyses/six_node_pattern S1/obs, nolegacy) whose M1 and M2 are both
                   perturbed members.
Primary: median(>= 2 members) - median(exactly 1), two-sided Wilcoxon, SE from 1,000 bootstrap resamples (seed 20250908),
pooled by REML over the human datasets of the five primary clusters.  Decision: SUPPORTED if the pooled difference is
> 0 with a 95 % CI excluding 0 under both (a) and (b); OPPOSITE if either CI excludes 0 below 0; NOT SUPPORTED otherwise.
Secondary: (i) linear model effect ~ [>= 2 members] + total sites of the members' families + log10(3' UTR length),
among genes targeted by >= 1 member (TargetScan class definition; total sites = conserved + non-conserved, UTR length of
the TargetScan representative human transcript); the class coefficient pooled by REML; (ii) six-node class vs 1-member
targets (network edges).  Human secondary-cluster datasets are reported separately.
Writes p2_datasets.tsv (one row per dataset and test), p2_pooled.tsv and p2_summary.txt (numbered lines)."""
import os, sys, glob, json, zipfile, io, re, subprocess
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(os.path.dirname(INB)); sys.path.insert(0, INB)
import _dl
from _seed import rng_for
D = f'{INB}/P2'
E = pd.read_csv(f'{REV}/data/canonical_edges.tsv', sep='\t'); E = E[E.edge_type == 'miRNA_target']
NET = E.groupby('source').target.apply(set).to_dict()
TS = pd.read_csv(os.path.join(_dl.CACHE, 'common/targetscan/ts80_human_family_sites.tsv.gz'), sep='\t')
FI = pd.read_csv(io.BytesIO(zipfile.ZipFile(os.path.join(_dl.CACHE, 'common/targetscan/miR_Family_Info.txt.zip')).read('miR_Family_Info.txt')), sep='\t')
FI = FI[FI['Species ID'] == 9606]
def seeds(node):
    base = node.replace('hsa-', '')
    rx = re.compile(r'^hsa-' + re.escape(base) + r'(-[35]p)?$')
    f = FI[FI['MiRBase ID'].str.match(rx)]
    c = f[f['Family Conservation?'] >= 1]
    return sorted(set((c if len(c) else f)['Seed+m8']))
SITE = {s: set(g[(g.cons8 + g.cons7m8) > 0].gene) for s, g in TS.groupby('family')}
TOT = TS.assign(tot=TS.cons_total + TS.noncons_total).groupby(['family', 'gene'], as_index=False).tot.max()   # one value per gene (max over transcripts)
UTRL = None
utr_zip = os.path.join(_dl.CACHE, 'common/targetscan/UTR_Sequences.txt.zip')
if os.path.exists(utr_zip):
    z = zipfile.ZipFile(utr_zip); lens = {}
    with z.open(z.namelist()[0]) as f:
        for line in io.TextIOWrapper(f, encoding='utf-8', errors='replace'):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 5 and p[3] == '9606': lens[p[2]] = max(lens.get(p[2], 0), len(p[4].replace('-', '')))
    UTRL = pd.Series(lens)
names = open(f'{REV}/analyses/six_node_pattern/S1/node_names.txt').read().split('\n')
BH = pd.read_csv(f'{REV}/analyses/six_node_pattern/S1/obs/nolegacy_bhat6_composite_instances.tsv', sep='\t', header=None, names=['m1', 'm2', 't1', 't2', 'g1', 'g2'])
for c in BH.columns: BH[c] = BH[c].map(lambda i: names[i])
def dmed(a, b, key=''):
    rng = rng_for(key)
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or len(b) < 3: return np.nan, np.nan, np.nan, len(a), len(b)
    d = np.median(a) - np.median(b)
    se = float(np.std([np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(1000)], ddof=1))
    return float(d), se, float(mannwhitneyu(a, b, alternative='two-sided').pvalue), len(a), len(b)
rows = []
for mp in sorted(glob.glob(f'{D}/de/*.meta.json')):
    m = json.load(open(mp)); t = pd.read_csv(mp.replace('.meta.json', '.tsv'), sep='\t').set_index('gene')
    members = m['mirna'].split(';'); sgn = -1.0 if m['direction'] == 'gain' else 1.0
    ex = t[t.expressed.astype(str) == 'True'].logFC.dropna(); eff = sgn * ex
    net_n = pd.Series(0, index=eff.index); ts_n = pd.Series(0, index=eff.index); tot = pd.Series(0.0, index=eff.index)
    seedset = []
    for mem in members:
        net_n += eff.index.isin(NET.get(mem, set())).astype(int)
        sd = seeds(mem); seedset += sd
        hit = set().union(*[SITE.get(s, set()) for s in sd]) if sd else set()
        ts_n += eff.index.isin(hit).astype(int)
    for s in sorted(set(seedset)):
        x = TOT[TOT.family == s].set_index('gene').tot; tot = tot.add(x.reindex(eff.index).fillna(0), fill_value=0)
    six = set(BH[BH.m1.isin(members) & BH.m2.isin(members)].g1)
    base = dict(gse=m['gse'], cluster=m['tf'], cell=m['cell'], cls=m['cls'], direction=m['direction'], members=';'.join(members),
                seeds=';'.join(sorted(set(seedset))), n_expressed=len(eff))
    for lab, n in (('network edges', net_n), ('TargetScan', ts_n)):
        d, se, p, n2, n1 = dmed(eff[n >= 2], eff[n == 1], f"{m['gse']}|{m['cell']}|{lab}")
        rows.append(dict(base, test=f'>=2 vs 1 members ({lab})', diff=d, se=se, wilcoxon_p=p, n_group=n2, n_ref=n1, n_zero=int((n == 0).sum())))
    d, se, p, n6, n1 = dmed(eff[eff.index.isin(six)], eff[net_n == 1], f"{m['gse']}|{m['cell']}|six")
    rows.append(dict(base, test='six-node class vs 1-member (network edges)', diff=d, se=se, wilcoxon_p=p, n_group=n6, n_ref=n1, n_zero=np.nan))
    if UTRL is not None:
        g = pd.DataFrame({'eff': eff, 'two': (ts_n >= 2).astype(int), 'tot': tot, 'utr': np.log10(UTRL.reindex(eff.index))})
        g = g[(ts_n >= 1) & g.utr.notna() & np.isfinite(g.utr)]
        if len(g) > 20 and g.two.nunique() == 2:
            fit = sm.OLS(g.eff, sm.add_constant(g[['two', 'tot', 'utr']])).fit()
            rows.append(dict(base, test='site-adjusted: >=2 members coefficient (TargetScan)', diff=float(fit.params['two']), se=float(fit.bse['two']),
                             wilcoxon_p=float(fit.pvalues['two']), n_group=int(g.two.sum()), n_ref=int((g.two == 0).sum()), n_zero=np.nan))
R = pd.DataFrame(rows); R.to_csv(f'{D}/p2_datasets.tsv', sep='\t', index=False)
pin = []
for grp, cl in (('human primary', 'human primary'), ('human secondary', 'human secondary'), ('mouse', 'mouse')):
    for test in R.test.unique():
        for r in R[(R.cls == cl) & (R.test == test)].itertuples(): pin.append(dict(group=f'{grp}|{test}', yi=r.diff, sei=r.se))
PT = pd.DataFrame()
if pin:
    pd.DataFrame(pin).to_csv(f'{D}/p2_pool_in.tsv', sep='\t', index=False)
    subprocess.run(['Rscript', f'{INB}/_rma.R', f'{D}/p2_pool_in.tsv', f'{D}/p2_pooled.tsv'], check=True)
    PT = pd.read_csv(f'{D}/p2_pooled.tsv', sep='\t')
L = []; P = L.append
for grp in ('human primary', 'human secondary', 'mouse'):
    x = R[R.cls == grp].drop_duplicates(['gse', 'cell'])
    P(f'{grp} datasets: {len(x)}' + (': ' + '; '.join(f'{r.gse} {r.cell} [{r.cluster}; {r.direction}; members {r.members}]' for r in x.itertuples()) if len(x) else ''))
def get(grp, test):
    x = PT[PT.group == f'{grp}|{test}'] if len(PT) else PT
    return x.iloc[0] if len(x) and x.iloc[0].k > 0 else None
for grp in ('human primary', 'human secondary', 'mouse'):
    for test in R.test.unique():
        x = get(grp, test)
        if x is None: P(f'{grp} | {test}: no dataset'); continue
        P(f'{grp} | {test}: k = {int(x.k)}; pooled {x.est:+.4f} (95 % CI {x.ci_lb:+.4f} to {x.ci_ub:+.4f}), p {x.p:.3g}; tau2 {x.tau2:.3g}')
a, b = get('human primary', '>=2 vs 1 members (network edges)'), get('human primary', '>=2 vs 1 members (TargetScan)')
if a is None or b is None: dec = 'NOT FOUND (no human primary-cluster dataset)'
elif (a.ci_ub < 0) or (b.ci_ub < 0): dec = 'OPPOSITE'
elif a.ci_lb > 0 and b.ci_lb > 0: dec = 'SUPPORTED'
else: dec = 'NOT SUPPORTED'
P(f'decision (human primary clusters, >=2 vs 1 members under network edges and TargetScan): {dec}')
for grp in ('human secondary', 'mouse'):
    a2, b2 = get(grp, '>=2 vs 1 members (network edges)'), get(grp, '>=2 vs 1 members (TargetScan)')
    d2 = ('NOT FOUND' if a2 is None or b2 is None else 'OPPOSITE' if (a2.ci_ub < 0 or b2.ci_ub < 0) else 'SUPPORTED' if (a2.ci_lb > 0 and b2.ci_lb > 0) else 'NOT SUPPORTED')
    P(f'same rule applied to the {grp} datasets (secondary, not the decision): {d2}')
s = get('human primary', 'site-adjusted: >=2 members coefficient (TargetScan)')
P('site-adjusted effect persists: ' + ('NOT ESTIMATED' if s is None else ('yes' if s.ci_lb > 0 else 'no') + f' (pooled coefficient {s.est:+.4f}, 95 % CI {s.ci_lb:+.4f} to {s.ci_ub:+.4f})'))
for r in R.itertuples():
    P(f'dataset {r.gse} {r.cell} [{r.cls}; {r.cluster}; {r.direction}] {r.test}: diff {r.diff:+.4f} (SE {r.se:.4f}), p {r.wilcoxon_p:.3g}; n {r.n_group} vs {r.n_ref}')
open(f'{D}/p2_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
