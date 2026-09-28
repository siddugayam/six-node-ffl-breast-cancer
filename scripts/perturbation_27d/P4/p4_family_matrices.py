#!/usr/bin/env python3
"""Requested P4 checks (a) and (b), made traceable: (a) the matrices of each family row of
results/v3/seqreg_ext_occlusion_calibrated.csv (families of scripts/v3/44_seqreg_ext_synthesis.py lines 14-16) at each
promoter; (b) the COL1A1 calls ranked 9-18 by the Enformer effect (d_fib_cage, most negative = rank 1) with their position
and span relative to the exon-1 donor (+222/+223).  Writes p4_family_matrices.txt (numbered lines)."""
import os
import pandas as pd
OUT = os.path.dirname(os.path.abspath(__file__)); REV = os.path.dirname(os.path.dirname(OUT))
a = pd.read_csv(f'{REV}/results/v3/seqreg_ext_occlusion_allsites.csv')
FAM = {"NF-kB": ["NFKB1", "NFKB2", "RELA", "RELB", "REL"], "SP": ["SP1", "SP2", "SP3"], "ETS": ["ETS1", "ETS2", "ELK1", "ELK4", "GABPA"]}
L = ['source: results/v3/seqreg_ext_occlusion_allsites.csv (Enformer, 1,626 calls); families: scripts/v3/44_seqreg_ext_synthesis.py lines 14-16']
for reg in ('COL1A1', 'COL3A1'):
    g = a[a.region == reg]
    for f, tfs in FAM.items():
        s = g[g.tf_name.isin(tfs)]
        L.append(f'(a) {reg} {f}: {len(s)} calls; matrices: ' + ('; '.join(f'{m} x{n}' for m, n in s.groupby("motif").size().items()) or 'none')
                 + f'; TFs present: {sorted(s.tf_name.unique())}; family members with no call: {sorted(set(tfs) - set(s.tf_name))}')
g = a[a.region == 'COL1A1'].copy(); g['rank'] = g.d_fib_cage.rank(method='min').astype(int)
for r in g.sort_values('d_fib_cage').iloc[8:18].itertuples():
    L.append(f'(b) COL1A1 rank {r.rank}: {r.motif} ({r.tf_name}) at TSS{r.rel_to_TSS:+d}, span {r.rel_to_TSS:+d}..{r.rel_to_TSS + r.site_len - 1:+d}, '
             f'{r.pct_change_fib_cage:+.2f} %; covers the donor (+222/+223): {"yes" if r.rel_to_TSS <= 223 and r.rel_to_TSS + r.site_len - 1 >= 222 else "no"}')
open(f'{OUT}/p4_family_matrices.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(f'{i + 1}: {l}' for i, l in enumerate(L)))
