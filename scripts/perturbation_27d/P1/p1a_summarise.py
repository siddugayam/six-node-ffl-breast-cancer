#!/usr/bin/env python3
"""P1a summary: per contrast (P1a/de/*.tsv): TF, COL1A1 and COL3A1 log2FC with SE; knockdown verification
(TF log2FC <= -0.74); COL1A1 expression in the controls (RNA-seq TPM >= 10; arrays percentile >= 50); positive control
(TRRUST v2 Activation targets of the TF vs expressed genes that are not TRRUST targets of the TF in any mode: difference
of median log2FC, two-sided Wilcoxon rank-sum).  Eligible primary contrasts are pooled per TF and cell class by REML
(_rma.R).  Decision rule of the analysis plan on the fibroblast-lineage pool: SUPPORTED if the pooled fold change and its 95 % CI
lie within 0.80-1.25; REFUTED if the CI lies wholly below 0.80 or wholly above 1.25; INCONCLUSIVE otherwise; NOT FOUND if
no eligible fibroblast-lineage dataset.  ETS1 carries the decision; NFKB1, RELA and SP1 are 'descriptive'.
Writes p1a_contrasts.tsv, p1a_pooled.tsv and p1a_summary.txt (numbered lines).
Optional: --first-pass leaves out the contrasts added by the second reading of the hits (P1a/p1a_triage_reread.tsv)
and writes the same files with the suffix _first_pass (SETTINGS change 12: both versions are reported)."""
import os, sys, glob, json, subprocess
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REV = os.path.dirname(INB)
D = f'{INB}/P1a'
FIRST = '--first-pass' in sys.argv; SUF = '_first_pass' if FIRST else ''
REREAD = set(pd.read_csv(f'{D}/p1a_triage_reread.tsv', sep='\t').gse) if FIRST else set()
TR = pd.read_csv(f'{REV}/data/db/trrust_human.tsv', sep='\t', header=None, names=['tf', 'target', 'mode', 'pmid'])
rows = []
for mp in sorted(glob.glob(f'{D}/de/*.meta.json')):
    m = json.load(open(mp))
    if m['gse'] in REREAD: continue
    t = pd.read_csv(mp.replace('.meta.json', '.tsv'), sep='\t').set_index('gene')
    tf = m['tf']
    def g(sym, col):
        return float(t.loc[sym, col]) if sym in t.index and pd.notna(t.loc[sym, col]) else np.nan
    r = dict(gse=m['gse'], tf=tf, cell=m['cell'], cls=m['cls'], tag=m['tag'], route=m['route'], platform=m['platform'],
             n_pert=m['n_pert'], n_ctrl=m['n_ctrl'])
    for sym in (tf, 'COL1A1', 'COL3A1'):
        k = 'TF' if sym == tf else sym
        r[f'{k}_log2FC'], r[f'{k}_SE'], r[f'{k}_P'] = g(sym, 'logFC'), g(sym, 'SE'), g(sym, 'P')
    r['kd_verified'] = 'yes' if r['TF_log2FC'] <= -0.74 else ('no' if np.isfinite(r['TF_log2FC']) else 'TF not measured')
    if m['route'] == 'ncbi':
        v = g('COL1A1', 'ctrl_TPM'); r['COL1A1_ctrl'] = f'{v:.3g} TPM'; r['col1a1_expressed'] = 'yes' if v >= 10 else 'no'
    else:
        v = g('COL1A1', 'ctrl_pct'); r['COL1A1_ctrl'] = f'percentile {v:.1f}'; r['col1a1_expressed'] = 'yes' if v >= 50 else 'no'
    act = set(TR[(TR.tf == tf) & (TR['mode'] == 'Activation')].target); anyt = set(TR[TR.tf == tf].target)
    ex = t[t.expressed == True] if t.expressed.dtype == bool else t[t.expressed.astype(str) == 'True']
    a = ex[ex.index.isin(act)].logFC.dropna(); b = ex[~ex.index.isin(anyt)].logFC.dropna()
    r['posctrl_n_targets'] = len(a); r['posctrl_diff_median'] = float(a.median() - b.median()) if len(a) else np.nan
    r['posctrl_wilcoxon_p'] = float(mannwhitneyu(a, b, alternative='two-sided').pvalue) if len(a) >= 3 else np.nan
    r['eligible'] = 'yes' if (r['kd_verified'] == 'yes' and r['col1a1_expressed'] == 'yes') else 'no'
    r['note'] = m.get('note', '')
    rows.append(r)
C = pd.DataFrame(rows); C.to_csv(f'{D}/p1a_contrasts{SUF}.tsv', sep='\t', index=False)
L = []; P = L.append
pool_in = []
for (tf, cls), g in C[(C.tag == 'primary') & (C.eligible == 'yes')].groupby(['tf', 'cls']):
    for gene in ('COL1A1', 'COL3A1'):
        for _, x in g.iterrows(): pool_in.append(dict(group=f'{tf}|{cls}|{gene}', yi=x[f'{gene}_log2FC'], sei=x[f'{gene}_SE']))
PT = pd.DataFrame()
if pool_in:
    pd.DataFrame(pool_in).to_csv(f'{D}/p1a_pool_in{SUF}.tsv', sep='\t', index=False)
    subprocess.run(['Rscript', f'{INB}/_rma.R', f'{D}/p1a_pool_in{SUF}.tsv', f'{D}/p1a_pooled{SUF}.tsv'], check=True)
    PT = pd.read_csv(f'{D}/p1a_pooled{SUF}.tsv', sep='\t')
def decide(est, lb, ub):
    if not np.isfinite(est): return 'NOT FOUND'
    f, fl, fu = 2 ** est, 2 ** lb, 2 ** ub
    if fl >= 0.80 and fu <= 1.25: return 'SUPPORTED'
    if fu < 0.80 or fl > 1.25: return 'REFUTED'
    return 'INCONCLUSIVE'
for tf in ('ETS1', 'NFKB1', 'RELA', 'SP1'):
    for gene in ('COL1A1', 'COL3A1'):
        for cls in ('fibroblast lineage', 'other'):
            x = PT[PT.group == f'{tf}|{cls}|{gene}'] if len(PT) else PT
            if len(x) == 0 or x.k.iloc[0] == 0:
                el = C[(C.tag == 'primary') & (C.eligible == 'yes') & (C.tf == tf) & (C.cls == cls)]
                why = (f'{len(el)} eligible dataset(s), none with a {gene} value (removed by filterByExpr or not measured: ' + ', '.join(el.gse) + ')') if len(el) else 'no eligible dataset'
                P(f'{tf} -> {gene}, {cls}: {why}' + (' -> NOT FOUND' if cls == 'fibroblast lineage' else '')); continue
            x = x.iloc[0]; dec = decide(x.est, x.ci_lb, x.ci_ub)
            lab = ('decision' if tf == 'ETS1' else 'descriptive') if cls == 'fibroblast lineage' else 'no decision (other cell types)'
            P(f'{tf} -> {gene}, {cls}: k = {int(x.k)}; pooled log2FC {x.est:+.3f} (95 % CI {x.ci_lb:+.3f} to {x.ci_ub:+.3f}); fold change {2**x.est:.3f} '
              f'({2**x.ci_lb:.3f}-{2**x.ci_ub:.3f}); tau2 {x.tau2:.3g}; rule: {dec} [{lab}]')
for _, r in C.iterrows():
    P(f"contrast {r.gse} {r.tf} {r.cell} [{r.tag}; {r.cls}; {r.route}; {r.n_pert} vs {r.n_ctrl}]: TF log2FC {r.TF_log2FC:+.2f} (knockdown verified: {r.kd_verified}); "
      f"COL1A1 {r.COL1A1_log2FC:+.3f} (SE {r.COL1A1_SE:.3f}); COL3A1 {r.COL3A1_log2FC:+.3f} (SE {r.COL3A1_SE:.3f}); COL1A1 in controls {r.COL1A1_ctrl} "
      f"(expressed: {r.col1a1_expressed}); positive control: {int(r.posctrl_n_targets)} activation targets, median shift {r.posctrl_diff_median:+.3f}, "
      f"Wilcoxon p {r.posctrl_wilcoxon_p:.2g}; eligible: {r.eligible}")
open(f'{D}/p1a_summary{SUF}.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
