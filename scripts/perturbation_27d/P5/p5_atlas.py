#!/usr/bin/env python3
"""P5: mature miR-29a/b/c-3p and miR-101-3p in purified cell types, fibroblasts against mammary epithelial cells, in two
atlases (SETTINGS.md, P5).  FANTOM5: CPM as provided; whole-cell samples only (nuclear/cytoplasmic fractions excluded).
microRNAome (Bioconductor 1.30.0): counts -> CPM = count / total miRNA counts of the sample x 1e6.
Groups: fibroblasts (description or CellType names a fibroblast); mammary fibroblasts; mammary epithelial cells.
Group value = mean over samples; ratio = fibroblast mean / mammary-epithelial mean.
Decision per miRNA: FIBROBLAST-ENRICHED if ratio >= 2 in every atlas with both groups; EPITHELIAL-ENRICHED if <= 0.5 in
every such atlas; otherwise NO CONSISTENT DIFFERENCE.
Writes p5_groups.tsv (one row per atlas x miRNA x group), p5_other_cell_types.tsv (one row per atlas x miRNA x cell
type, mean CPM) and p5_summary.txt (numbered lines).
Correction (SETTINGS.md change 16, made after the P5 result had been delivered): the version first delivered missed
two kinds of sample.  FANTOM5 marks some fractions "(cytosolic)" or "(nuclear)" rather than "(... fraction)"; these are
now excluded as well, which removes two dermal-fibroblast fractions from the fibroblasts.  The microRNAome CellType
"iPSC_fibroblast" (Class "Stem": induced pluripotent stem cells) names a fibroblast but is not one; its 32 samples are
now excluded from the fibroblasts.  --as-delivered reproduces the first version and writes *_as_delivered files."""
import os, sys, re
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _dl
OUT = os.path.dirname(os.path.abspath(__file__))
MIRS = ['hsa-miR-29a-3p', 'hsa-miR-29b-3p', 'hsa-miR-29c-3p', 'hsa-miR-101-3p']
F = os.path.join(_dl.CACHE, 'P5/fantom5'); M = os.path.join(_dl.CACHE, 'P5/microRNAome')
L = []; P = L.append; rows = []; other = []
AD = '--as-delivered' in sys.argv; SUF = '_as_delivered' if AD else ''

# ---------------- FANTOM5
S = pd.read_csv(f'{F}/human.srna.samples.tsv', sep='\t')
X = pd.read_csv(f'{F}/human.srna.cpm.txt', sep='\t', index_col=0)
S = S[S.name.isin(X.columns)]
frac = S.description.str.contains(r'\((nuclear|cytoplasmic) fraction\)' if AD else r'\((?:nuclear|cytoplasmic|cytosolic)(?: fraction)?\)', case=False)
W = S[~frac].copy()
W['cell'] = W.description.str.replace(r',\s*(donor|pool|rep|biol_rep|tech_rep).*$', '', regex=True).str.strip()
groups = {'fibroblasts (all)': W.description.str.match(r'Fibroblast', case=False),
          'mammary fibroblasts': W.description.str.startswith('Fibroblast - Mammary'),
          'mammary epithelial cells': W.description.str.startswith('Mammary Epithelial Cell')}
P(f'FANTOM5: {len(S)} samples in the CPM table; {int(frac.sum())} subcellular-fraction samples excluded; '
  + '; '.join(f'{g}: {int(m.sum())}' for g, m in groups.items()))
missing = [m for m in MIRS if m not in X.index]
if missing: P(f'FANTOM5: NOT FOUND in the CPM table: {missing}')
for mir in [m for m in MIRS if m in X.index]:
    for g, m in groups.items():
        v = X.loc[mir, W[m].name]
        rows.append(dict(atlas='FANTOM5', mirna=mir, group=g, n=len(v), mean_cpm=float(v.mean()), median_cpm=float(v.median())))
    for c, idx in W.groupby('cell').groups.items():
        other.append(dict(atlas='FANTOM5', mirna=mir, cell_type=c, n=len(idx), mean_cpm=float(X.loc[mir, W.loc[idx, 'name']].mean())))

# ---------------- microRNAome
D = pd.read_csv(f'{M}/microRNAome_4mirs.tsv', sep='\t')
cpm = D[MIRS].div(D.total_miRNA_counts, axis=0) * 1e6
ct = D.CellType.astype(str)
ipsc = D.Class.astype(str).eq('Stem') & ct.str.contains('Fibroblast', case=False)
groups = {'fibroblasts (all)': ct.str.contains('Fibroblast', case=False) & (True if AD else ~ipsc),
          'mammary fibroblasts': ct.eq('Fibroblast_breast'),
          'mammary epithelial cells': ct.eq('Breast_epithelial_cell')}
P(f'microRNAome 1.30.0: {len(D)} samples; ' + '; '.join(f'{g}: {int(m.sum())}' for g, m in groups.items())
  + f'; samples with zero total miRNA counts: {int((D.total_miRNA_counts == 0).sum())}'
  + ('' if AD else f'; excluded from the fibroblasts: {int(ipsc.sum())} samples of Class "Stem" '
     f'(CellType {", ".join(sorted(set(ct[ipsc])))})'))
for mir in MIRS:
    for g, m in groups.items():
        v = cpm.loc[m, mir]
        rows.append(dict(atlas='microRNAome', mirna=mir, group=g, n=int(m.sum()), mean_cpm=float(v.mean()), median_cpm=float(v.median())))
    for c, idx in D.groupby('CellType').groups.items():
        other.append(dict(atlas='microRNAome', mirna=mir, cell_type=c, n=len(idx), mean_cpm=float(cpm.loc[idx, mir].mean())))

G = pd.DataFrame(rows); G.to_csv(f'{OUT}/p5_groups{SUF}.tsv', sep='\t', index=False)
pd.DataFrame(other).to_csv(f'{OUT}/p5_other_cell_types{SUF}.tsv', sep='\t', index=False)
for mir in MIRS:
    ratios = []
    for atlas in ('FANTOM5', 'microRNAome'):
        g = G[(G.atlas == atlas) & (G.mirna == mir)].set_index('group')
        if len(g) == 0 or g.loc['mammary epithelial cells', 'n'] == 0 or g.loc['fibroblasts (all)', 'n'] == 0:
            P(f'{mir} {atlas}: not both groups'); continue
        fa, fm, ep = g.loc['fibroblasts (all)', 'mean_cpm'], g.loc['mammary fibroblasts', 'mean_cpm'], g.loc['mammary epithelial cells', 'mean_cpm']
        r = fa / ep if ep > 0 else np.inf; rm = fm / ep if ep > 0 else np.inf; ratios.append(r)
        P(f'{mir} {atlas}: fibroblasts (all, n = {int(g.loc["fibroblasts (all)", "n"])}) mean CPM {fa:.4g}; mammary fibroblasts (n = '
          f'{int(g.loc["mammary fibroblasts", "n"])}) {fm:.4g}; mammary epithelial cells (n = {int(g.loc["mammary epithelial cells", "n"])}) {ep:.4g}; '
          f'ratio fibroblast/mammary epithelial {r:.3g} (mammary fibroblast/mammary epithelial {rm:.3g})')
    dec = ('FIBROBLAST-ENRICHED' if ratios and all(x >= 2 for x in ratios) else 'EPITHELIAL-ENRICHED' if ratios and all(x <= 0.5 for x in ratios)
           else 'NO CONSISTENT DIFFERENCE')
    P(f'{mir} decision ({len(ratios)} atlases with both groups): {dec}')
open(f'{OUT}/p5_summary{SUF}.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
