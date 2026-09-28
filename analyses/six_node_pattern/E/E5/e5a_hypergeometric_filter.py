#!/usr/bin/env python3
"""
E5(a): which TF list entered the hypergeometric miRNA-TF pair filter, and are any of the 22 TF-typed
nodes outside the Lambert census among the retained TFs?

Input: the author-supplied zip 'hypergeometrictestoutput.zip' (2026-09-27; Hypergeometric_Test_2.xlsx,
Ntf.txt, Nmir.txt, tf-gene.txt, mirna-gene.txt, mirna-tf.txt, value.txt, common-i.txt), unpacked outside
the project; its md5s are printed.  Project file used for comparison: data/db/trrust_human.tsv (TRRUST v2).
Workbook layout: left block = tested pairs (miRNA, TF, Nmir, Ntf, i = common genes, P-value); right block =
the same p-values ranked with the Benjamini-Hochberg critical value, adjusted p and 'Significant using an
FDR of 0.05?' (Yes = retained).
Output: e5a_hypergeometric_filter.txt, e5a_tested_pairs.csv (one row per tested miRNA-TF pair, with the
workbook p-value, the recomputed lower- and upper-tail hypergeometric probabilities and the retained flag).
usage: e5a_hypergeometric_filter.py <unpacked zip dir>
"""
import sys, os, hashlib, csv, collections
import numpy as np, pandas as pd
from scipy.stats import hypergeom
D = sys.argv[1]; H = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
out = []
P = lambda *a: (out.append(' '.join(str(x) for x in a)), print(*a))
for f in sorted(os.listdir(D)):
    if f.endswith(('.txt', '.xlsx')): P(f'md5 {hashlib.md5(open(os.path.join(D, f), "rb").read()).hexdigest()}  {f}')
TF_list = {l.split('\t')[0] for l in open(f'{D}/Ntf.txt') if l.strip()}
tg = {tuple(l.rstrip('\n').split('\t')[:2]) for l in open(f'{D}/tf-gene.txt') if l.strip()}
tr = [l.rstrip('\n').split('\t') for l in open(f'{REV}/data/db/trrust_human.tsv') if l.strip()]
trTF = {r[0] for r in tr}; trp = {(r[0], r[1]) for r in tr}
P(f'TF list entering the filter (Ntf.txt = regulators of tf-gene.txt): {len(TF_list)} TFs; identical to the {len(trTF)} '
  f'regulators of TRRUST v2 human (data/db/trrust_human.tsv): {TF_list == trTF}; TF-target pairs {len(tg)} vs TRRUST distinct '
  f'pairs {len(trp)}; pairs differing only by Excel date conversion of target symbols: {sorted(tg - trp)} / {sorted(trp - tg)}')
d = pd.read_excel(f'{D}/Hypergeometric_Test_2.xlsx', header=None).iloc[2:]
Lb = d[[1, 2, 3, 4, 5, 6]].dropna(how='all').copy(); Lb.columns = ['miRNA', 'TF', 'Nmir', 'Ntf', 'common', 'p']
for c in ['Nmir', 'Ntf', 'common', 'p']: Lb[c] = pd.to_numeric(Lb[c])
Lb['TF'] = Lb.TF.astype(str); Lb['miRNA'] = Lb.miRNA.astype(str)
Rb = d[[8, 9, 10, 11, 12, 13, 14]].dropna(how='all').copy(); Rb.columns = ['rank', 'p', 'miRNA', 'TF', 'crit', 'padj', 'sig']
Rb['TF'] = Rb.TF.astype(str); Rb['miRNA'] = Rb.miRNA.astype(str)
yes = Rb[Rb.sig.astype(str).str.strip().str.lower() == 'yes']
P(f'tested pairs {len(Lb)} ({Lb.TF.nunique()} TF entries); retained (Yes) {len(yes)} pairs with {yes.TF.nunique()} distinct TF entries, '
  f'of which non-symbol entries (Excel date serials): {sorted(t for t in yes.TF.unique() if t.isdigit())}')
T22 = 'APEX1 BMI1 BRCA1 CREBBP CTNNB1 EP300 EZH2 HDAC1 HDAC2 HDAC3 HDAC4 HDAC9 ILF3 MEN1 MKL1 MTA1 NCOR1 NF1 RB1 SIRT1 SUZ12 VHL'.split()
P(f'of the 22: in the filter TF list {sum(t in TF_list for t in T22)}; in a tested pair {sum(t in set(Lb.TF) for t in T22)}; '
  f'among the retained TFs {len(set(T22) & set(yes.TF))}: {sorted(set(T22) & set(yes.TF))}; not retained: {sorted(set(T22) - set(yes.TF))}')
# retained TF entries against the filter TF list (the pair list, mirna-tf.txt, is not restricted to it)
rin = sorted(t for t in yes.TF.unique() if t in TF_list); rout = sorted(t for t in yes.TF.unique() if t not in TF_list)
P(f'retained TF entries in the 795-TF filter list (Ntf.txt): {len(rin)}; not in it: {len(rout)} ({sum(not t.isdigit() for t in rout)} symbols '
  f'and {sum(t.isdigit() for t in rout)} Excel date serial), in {int((~yes.TF.isin(TF_list)).sum())} of the {len(yes)} retained pairs: {rout}; '
  f'tested pairs with a TF entry outside the list: {int((~Lb.TF.isin(TF_list)).sum())} of {len(Lb)} ({Lb.TF[~Lb.TF.isin(TF_list)].nunique()} entries)')
# the p-value formula
N = 256
lo = hypergeom.cdf(Lb.common, N, Lb.Nmir, Lb.Ntf); up = hypergeom.sf(Lb.common - 1, N, Lb.Nmir, Lb.Ntf)
ok = Lb.p.notna()
P(f'workbook P-value = lower-tail cumulative hypergeometric P(X <= i) with a universe of N = {N} genes: max |difference| over '
  f'the {int(ok.sum())} rows with a P-value = {np.nanmax(np.abs(lo[ok] - Lb.p[ok])):.2e} ({int((~ok).sum())} rows have no P-value)')
key = set(yes.miRNA + '|' + yes.TF); Lb['retained'] = (Lb.miRNA + '|' + Lb.TF).isin(key)
Y = Lb[Lb.retained]; ex = Y.Nmir * Y.Ntf / N
P(f'retained pairs matched to the input block: {len(Y)}; sharing 0 targets: {int((Y.common == 0).sum())}; sharing >= 1: '
  f'{int((Y.common >= 1).sum())}; sharing more targets than expected (Nmir*Ntf/N): {int((Y.common > ex).sum())}; '
  f'upper-tail (enrichment) p < 0.05: {int((hypergeom.sf(Y.common - 1, N, Y.Nmir, Y.Ntf) < 0.05).sum())}')
P(f'all tested pairs with upper-tail (enrichment) p < 0.05 at N = {N}: {int((up < 0.05).sum())} of {len(Lb)}')
Lb.assign(p_lower_tail_recomputed=lo, p_upper_tail_enrichment=up).to_csv(f'{H}/e5a_tested_pairs.csv', index=False)
open(f'{H}/e5a_hypergeometric_filter.txt', 'w').write('\n'.join(out) + '\n')
