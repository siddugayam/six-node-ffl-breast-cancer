#!/usr/bin/env python3
"""E7: direct (ADE) and mediated (ACME) effects with 95 % percentile-bootstrap CIs (2,000 resamples) for
the four TF -> COL1A1 arms (ETS1, NFKB1, RELA, SP1) and the three miR-29 -> COL1A1 arms, one row per arm and
stromal estimate, from results/v3/deconv_mediation_extended.csv (scripts/08_stroma_and_cell_types/deconvolution_and_mediation/46_deconv2_mediation.R).
The stored run has 45 stromal estimates per arm; the 42 estimates of the manuscript leave out three that repeat
another estimate's results to the last digit.  Rows are flagged so that the 42-estimate subset can be
recounted: 'exact_duplicate_of' marks the three estimates whose results equal another
estimate's to the last digit: MCPcounter_Fibroblasts_recomputed (= MCPcounter_Fibroblasts), NNLS_Wu_CAFs
(= NNLS_Wu641_full_CAFs) and QPROG_Wu_CAFs (= QPROG_Wu641_full_CAFs), the later name alphabetically being
flagged; 'non_RNA' the six non-transcriptomic estimates.
inconsistent = ACME and ADE of opposite sign (the stored column); ADE_CI_excludes_zero = ade_lo > 0 or
ade_hi < 0.  Read-only."""
import os, pandas as pd
REV = '/path/to/revision'; H = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(f'{REV}/results/v3/deconv_mediation_extended.csv')
meta = pd.read_csv(f'{REV}/results/v3/deconv_mediator_metadata.csv').set_index('mediator')
arms = ['ETS1', 'NFKB1', 'RELA', 'SP1', 'hsa-miR-29a', 'hsa-miR-29b', 'hsa-miR-29c']
d = d[(d.y == 'COL1A1') & d.x.isin(arms)].copy()
NONRNA = {'ABSOLUTE_nontumour', 'PATH_stroma_frozen', 'PATH_stroma_allslides', 'METH_EpiDISH_Fib', 'METH_EpiDISH_Epi', 'METH_EpiDISH_IC'}
# exact duplicates: identical ACME, ADE and CIs on every arm
key = d.groupby('mediator').apply(lambda g: tuple(g.sort_values('x')[['ACME', 'ADE', 'acme_lo', 'acme_hi', 'ade_lo', 'ade_hi']].round(15).values.ravel()))
first = {}; dup = {}
for m in sorted(key.index):
    k = key[m]
    if k in first: dup[m] = first[k]
    else: first[k] = m
print('exact duplicates (estimate -> duplicated estimate):', dup)
d['estimate_family'] = d.mediator.map(meta['family']); d['non_RNA'] = d.mediator.isin(NONRNA)
d['exact_duplicate_of'] = d.mediator.map(dup).fillna('')
d['ADE_CI_excludes_zero'] = (d.ade_lo > 0) | (d.ade_hi < 0)
d['ACME_CI_excludes_zero'] = (d.acme_lo > 0) | (d.acme_hi < 0)
cols = ['x', 'y', 'mediator', 'estimate_family', 'non_RNA', 'exact_duplicate_of', 'n', 'ACME', 'acme_lo', 'acme_hi', 'ADE', 'ade_lo',
        'ade_hi', 'prop_mediated', 'prop_lo', 'prop_hi', 'sobel_p', 'inconsistent', 'ADE_CI_excludes_zero', 'ACME_CI_excludes_zero']
d = d[cols].rename(columns={'x': 'arm_regulator', 'y': 'target', 'mediator': 'stromal_estimate'}).sort_values(['arm_regulator', 'stromal_estimate'])
d.to_csv(f'{H}/e7_mediation_intervals_45_per_arm.csv', index=False)
rows = []
for a, g in d.groupby('arm_regulator'):
    for lab, gg in (('all 45', g), ('42 = 45 minus the 3 exact duplicates', g[g.exact_duplicate_of == '']), ('39 transcriptomic', g[~g.non_RNA])):
        rows.append(dict(arm=f'{a} -> COL1A1', estimate_set=lab, n_estimates=len(gg), n_inconsistent=int(gg.inconsistent.sum()),
                         n_inconsistent_with_ADE_CI_excluding_zero=int((gg.inconsistent & gg.ADE_CI_excludes_zero).sum()),
                         n_ADE_CI_excluding_zero=int(gg.ADE_CI_excludes_zero.sum())))
S = pd.DataFrame(rows); S.to_csv(f'{H}/e7_inconsistent_counts.csv', index=False)
pd.set_option('display.width', 200); print(S.to_string(index=False))
