#!/usr/bin/env python3
"""R1 (analyses/perturbation_sensitivity): rebuild of results/v6/atac_encode_compartment_specificity_tests.csv, for which no saved script
exists (it was written on 2026-09-10 at 10:25, one minute after scripts/09_regulatory_evidence/encode_accessibility/atac_06_encode_compartment.py last changed, by a
step that was not kept).  The rule is read off the stored file and checked against it:
  input   results/v6/atac_encode_region_accessible.csv (atac_06's 0/1 matrix: query region x ENCODE sample)
  samples epithelial = MCF-7, MCF 10A, breast epithelium (5); fibroblast = IMR-90 and fibroblasts of mammary gland,
          dermis and lung (5); the compartments of atac_06 (EPI / FIB lists)
  windows TSS+/-100bp = the <gene>_TSS regions; TSS+/-1kb = the <gene>_prom regions (atac_06 lines 24-25)
  test    two-sided Fisher exact test on [[fib accessible, fib not], [epi accessible, epi not]] per gene and window
  FDR     Benjamini-Hochberg across all 34 tests (both windows together)
  order   the genes in the matrix order, TSS+/-100bp rows first
Writes atac_encode_compartment_specificity_tests_rebuilt.csv here and reports the comparison with the stored file."""
import os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
HERE = os.path.dirname(os.path.abspath(__file__)); REV = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
A = pd.read_csv(f'{REV}/results/v6/atac_encode_region_accessible.csv', index_col=0)
EPI = ('MCF-7', 'MCF 10A', 'breast epithelium'); FIB = ('fibroblast of mammary gland', 'fibroblast of dermis', 'fibroblast of lung', 'IMR-90')
epi = [c for c in A.columns if c.split(' (')[0] in EPI]; fib = [c for c in A.columns if c.split(' (')[0] in FIB]
genes = list(dict.fromkeys(r.rsplit('_', 1)[0] for r in A.index))
rows = []
for window, suf in (('TSS+/-100bp', 'TSS'), ('TSS+/-1kb', 'prom')):
    part = []
    for g in genes:
        x = A.loc[f'{g}_{suf}']; fa, ea = int(x[fib].sum()), int(x[epi].sum())
        p = fisher_exact([[fa, len(fib) - fa], [ea, len(epi) - ea]], alternative='two-sided')[1]
        part.append(dict(gene=g, window=window, fib_accessible=fa, fib_n=len(fib), epi_accessible=ea, epi_n=len(epi), p_fisher=p))
    rows += part
P = np.array([r['p_fisher'] for r in rows]); m = len(P); o = np.argsort(P, kind='stable')
q = np.minimum.accumulate((P[o] * m / np.arange(1, m + 1))[::-1])[::-1]; fdr = np.empty(m); fdr[o] = np.minimum(q, 1)
for r, f in zip(rows, fdr): r['fdr_BH'] = f
D = pd.DataFrame(rows); D.to_csv(f'{HERE}/atac_encode_compartment_specificity_tests_rebuilt.csv', index=False)
S = pd.read_csv(f'{REV}/results/v6/atac_encode_compartment_specificity_tests.csv', float_precision='round_trip')
same_bytes = open(f'{HERE}/atac_encode_compartment_specificity_tests_rebuilt.csv', 'rb').read() == open(f'{REV}/results/v6/atac_encode_compartment_specificity_tests.csv', 'rb').read()
same_keys = (S[['gene', 'window', 'fib_accessible', 'fib_n', 'epi_accessible', 'epi_n']].values == D[['gene', 'window', 'fib_accessible', 'fib_n', 'epi_accessible', 'epi_n']].values).all()
dp = float(np.max(np.abs(S.p_fisher.values - D.p_fisher.values))); dq = float(np.max(np.abs(S.fdr_BH.values - D.fdr_BH.values)))
L = [f'samples: fibroblast {len(fib)}, epithelial {len(epi)}; rows {len(D)} (stored {len(S)})',
     f'genes, windows and counts identical to the stored file: {"yes" if same_keys else "NO"}',
     f'largest absolute difference: p_fisher {dp:.3g}, fdr_BH {dq:.3g}; byte-identical to the stored file: {"yes" if same_bytes else "no"}',
     'COL1A1 and COL3A1 at TSS+/-100bp: ' + '; '.join(f'{r.gene} {r.fib_accessible}/{r.fib_n} vs {r.epi_accessible}/{r.epi_n}, p {r.p_fisher:.4g}, FDR {r.fdr_BH:.4g}'
                                                      for r in D[(D.window == 'TSS+/-100bp') & D.gene.isin(['COL1A1', 'COL3A1'])].itertuples())]
open(f'{HERE}/r1_summary.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
