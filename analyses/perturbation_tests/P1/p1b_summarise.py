#!/usr/bin/env python3
"""P1b summary (SETTINGS.md P1b and change 6).
Per contrast (P1b/de/*.tsv): COL1A1, COL3A1 and the COL1A1/COL3A1 pseudo-gene log2FC and SE; the sign-aligned repression
effect (-log2FC for gain, +log2FC for loss); baseline expression of COL1A1 and COL3A1 in the controls (RNA-seq TPM,
arrays percentile); the target shift = median effect of expressed genes with >= 1 conserved 8mer/7mer-m8 site of the
family minus that of expressed genes with no site of any type for the family (TargetScan 8.0 Summary Counts, human),
SE from 1,000 bootstrap resamples of genes within each class (seed 20250908).
miR-29 family: seed AGCACCA (miR-29-3p); miR-101: seed ACAGUAC (miR-101-3p.1, the canonical miR-101-3p).
Decision (collagen pseudo-gene effect, primary): REML meta-regression on fibroblast lineage (1) vs carcinoma (0):
SUPPORTED if the fibroblast coefficient > 0 with p < 0.05; REFUTED if < 0 with p < 0.05; INCONCLUSIVE otherwise;
NOT FOUND if either class has no dataset.  The target shift is reported by the same rule (secondary).
Positive control: EZH2 log2FC pooled by REML over the miR-101 gain datasets (expected negative).
Writes p1b_contrasts.tsv, p1b_meta_regression.tsv and p1b_summary.txt (numbered lines)."""
import os, sys, glob, json, subprocess
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _dl
from _seed import rng_for
D = f'{INB}/P1b'
TS = pd.read_csv(os.path.join(_dl.CACHE, 'common/targetscan/ts80_human_family_sites.tsv.gz'), sep='\t')
SEED = {'miR-29': 'AGCACCA', 'miR-101': 'ACAGUAC'}
def site_sets(seed):
    t = TS[TS.family == seed]
    site = set(t[(t.cons8 + t.cons7m8) > 0].gene)
    anysite = set(t[(t.cons_total + t.noncons_total) > 0].gene)
    return site, anysite
SS = {k: site_sets(v) for k, v in SEED.items()}
def boot_shift(a, b, n=1000, key=''):
    rng = rng_for(key)
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 3 or len(b) < 3: return np.nan
    est = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return float(np.std(est, ddof=1))
rows = []
for mp in sorted(glob.glob(f'{D}/de/*.meta.json')):
    m = json.load(open(mp)); t = pd.read_csv(mp.replace('.meta.json', '.tsv'), sep='\t').set_index('gene')
    fam = 'miR-101' if 'miR-101' in m['tf'] else 'miR-29'
    sgn = -1.0 if m.get('direction', 'gain') == 'gain' else 1.0
    g = lambda s, c: float(t.loc[s, c]) if s in t.index and pd.notna(t.loc[s, c]) else np.nan
    r = dict(gse=m['gse'], mirna_set=fam, cell=m['cell'], cls=m['cls'], tag=m.get('tag', 'primary'), direction=m.get('direction', ''), route=m['route'], platform=m['platform'],
             n_pert=m['n_pert'], n_ctrl=m['n_ctrl'])
    for s, k in (('COL1A1', 'COL1A1'), ('COL3A1', 'COL3A1'), ('COL1A1_COL3A1', 'pair'), ('EZH2', 'EZH2')):
        r[f'{k}_log2FC'], r[f'{k}_SE'] = g(s, 'logFC'), g(s, 'SE')
        r[f'{k}_repression'] = sgn * r[f'{k}_log2FC']
    for s in ('COL1A1', 'COL3A1'):
        r[f'{s}_ctrl'] = (f"{g(s, 'ctrl_TPM'):.3g} TPM" if 'ctrl_TPM' in t.columns and np.isfinite(g(s, 'ctrl_TPM')) else f"percentile {g(s, 'ctrl_pct'):.1f}")
    ex = t[t.expressed.astype(str) == 'True'] if 'expressed' in t.columns else t
    site, anysite = SS[fam]
    eff = sgn * ex.logFC.dropna()
    a = eff[eff.index.isin(site)]; b = eff[~eff.index.isin(anysite)]
    r['n_site_genes'], r['n_nosite_genes'] = len(a), len(b)
    r['target_shift'] = float(a.median() - b.median()) if len(a) and len(b) else np.nan
    r['target_shift_SE'] = boot_shift(a.values, b.values, key=f"{m['gse']}|{m['cell']}|{m.get('tag', 'primary')}")
    r['note'] = m.get('note', '')[:200]
    rows.append(r)
C = pd.DataFrame(rows); C.to_csv(f'{D}/p1b_contrasts.tsv', sep='\t', index=False)
L = []; P = L.append
ALL29 = C[C.mirna_set == 'miR-29']
M29 = ALL29[ALL29.tag == 'primary']
# sensitivity set: the force-keep version where it exists (RNA-seq datasets), otherwise the primary contrast
SENS = pd.concat([g[g.tag == 'sensitivity'] if (g.tag == 'sensitivity').any() else g[g.tag == 'primary'] for _, g in ALL29.groupby(['gse', 'cell'])])
for cls in ('fibroblast lineage', 'carcinoma', 'other'):
    x = M29[M29.cls == cls]
    P(f'miR-29 datasets, {cls}: {len(x)}' + (': ' + '; '.join(f'{r.gse} {r.cell} ({r.direction})' for r in x.itertuples()) if len(x) else ''))
mr_in = []
for DS, tagname in ((M29, ''), (SENS, ' [sensitivity: collagens force-kept in RNA-seq]')):
    for col, se, lab in (('pair_repression', 'pair_SE', 'collagen pseudo-gene (primary)'), ('target_shift', 'target_shift_SE', 'target shift (secondary)'),
                         ('COL1A1_repression', 'COL1A1_SE', 'COL1A1 alone'), ('COL3A1_repression', 'COL3A1_SE', 'COL3A1 alone')):
        for r in DS[DS.cls.isin(['fibroblast lineage', 'carcinoma'])].itertuples():
            mr_in.append(dict(group=lab + tagname, yi=getattr(r, col), sei=getattr(r, se), mod=1 if r.cls == 'fibroblast lineage' else 0))
pd.DataFrame(mr_in).to_csv(f'{D}/p1b_mr_in.tsv', sep='\t', index=False)
subprocess.run(['Rscript', f'{INB}/_rma.R', f'{D}/p1b_mr_in.tsv', f'{D}/p1b_meta_regression.tsv'], check=True)
MR = pd.read_csv(f'{D}/p1b_meta_regression.tsv', sep='\t')
MI = pd.DataFrame(mr_in)
for lab in MR.group.unique():
    x = MR[MR.group == lab].set_index('term')
    mi = MI[(MI.group == lab) & np.isfinite(MI.yi) & np.isfinite(MI.sei)]
    nf, nc = int((mi['mod'] == 1).sum()), int((mi['mod'] == 0).sum())
    if 'mod' not in x.index or not np.isfinite(x.loc['mod', 'est']):
        P(f'{lab}: meta-regression not estimable (fibroblast k = {nf}, carcinoma k = {nc}) -> NOT FOUND'); continue
    b, p, i0 = x.loc['mod', 'est'], x.loc['mod', 'p'], x.loc['intercept', 'est']
    dec = 'NOT FOUND' if (nf == 0 or nc == 0) else ('SUPPORTED' if (b > 0 and p < 0.05) else 'REFUTED' if (b < 0 and p < 0.05) else 'INCONCLUSIVE')
    P(f'{lab}: carcinoma mean {i0:+.3f}; fibroblast minus carcinoma {b:+.3f} (95 % CI {x.loc["mod", "ci_lb"]:+.3f} to {x.loc["mod", "ci_ub"]:+.3f}), '
      f'z-test p {p:.3g}; k = {int(x.loc["mod", "k"])} with a value (fibroblast {nf}, carcinoma {nc}); tau2 {x.loc["mod", "tau2"]:.3g} -> {dec}'
      + ('' if 'sensitivity' not in lab else ' (sensitivity; not the decision)'))
for cls in ('fibroblast lineage', 'carcinoma', 'other'):
    x = M29[M29.cls == cls]
    if len(x): P(f'baseline in controls, {cls}: COL1A1 ' + '; '.join(f'{r.gse} {r.cell}: {r.COL1A1_ctrl}' for r in x.itertuples())
                 + ' | COL3A1 ' + '; '.join(f'{r.gse} {r.cell}: {r.COL3A1_ctrl}' for r in x.itertuples()))
M101 = C[(C.mirna_set == 'miR-101') & (C.tag == 'primary')]
if len(M101):
    pd.DataFrame([dict(group='EZH2', yi=r.EZH2_log2FC, sei=r.EZH2_SE) for r in M101.itertuples()]).to_csv(f'{D}/p1b_pc_in.tsv', sep='\t', index=False)
    subprocess.run(['Rscript', f'{INB}/_rma.R', f'{D}/p1b_pc_in.tsv', f'{D}/p1b_positive_control.tsv'], check=True)
    pc = pd.read_csv(f'{D}/p1b_positive_control.tsv', sep='\t').iloc[0]
    P(f'positive control, miR-101 gain -> EZH2: k = {int(pc.k)}; pooled log2FC {pc.est:+.3f} (95 % CI {pc.ci_lb:+.3f} to {pc.ci_ub:+.3f}), p {pc.p:.3g}; '
      + ('negative, as expected' if pc.ci_ub < 0 else 'not significantly negative'))
for r in C.itertuples():
    P(f'contrast {r.gse} {r.mirna_set} {r.cell} [{r.tag}; {r.cls}; {r.direction}; {r.route}; {r.n_pert} vs {r.n_ctrl}]: COL1A1 {r.COL1A1_log2FC:+.3f} (SE {r.COL1A1_SE:.3f}); '
      f'COL3A1 {r.COL3A1_log2FC:+.3f} (SE {r.COL3A1_SE:.3f}); pair repression {r.pair_repression:+.3f} (SE {r.pair_SE:.3f}); EZH2 {r.EZH2_log2FC:+.3f}; '
      f'target shift {r.target_shift:+.3f} (SE {r.target_shift_SE:.3f}; {r.n_site_genes} site vs {r.n_nosite_genes} no-site genes); '
      f'baseline COL1A1 {r.COL1A1_ctrl}, COL3A1 {r.COL3A1_ctrl}')
open(f'{D}/p1b_summary.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
