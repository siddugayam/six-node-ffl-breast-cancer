#!/usr/bin/env python3
"""Stability of the composite evidence ranking of network nodes (Results 3.3; Note S3).

The composite score of a node is the mean of its rank-normalised scores in the evidence domains in which it has
data (topology T, edge evidence E, expression D, replication R, multi-omic M, clinical C, functional F); a node is
ranked only if it has data in at least three domains. The thirty prioritised nodes are the ten highest-ranked
transcription factors, miRNAs and genes. Three perturbations test whether that list depends on the equal weights
or on any single domain:

  A  Leave-one-domain-out: each of the seven domains is removed in turn (the three-domain rule is kept).
  B  Random weights: 10,000 weight vectors drawn from a flat Dirichlet distribution over the seven domains; the
     composite of a node is the weighted mean over the domains in which it has data.
  C  Median of the available domains in place of the mean.

Per perturbation the script reports the overlap of the thirty-node list with the baseline list, the rank
correlation with the baseline ranking (Kendall tau over the ranked nodes) and, for every baseline node, the fraction
of perturbations in which it stays among the ten highest-ranked nodes of its class.

Inputs : results/prioritisation_and_stroma/node_prioritisation_full.csv and top10_per_class.csv. The rank columns of
         node_prioritisation_full.csv count the four nodes without a composite score first; ranks are therefore
         recomputed here from the composite score (Table S1 gives the ranks of the paper).
Outputs: results/prioritisation_and_stroma/composite_rank_robustness_summary.csv (one row per perturbation scheme
         and domain), composite_rank_robustness_nodes.csv (one row per baseline node), and
         composite_rank_robustness_baseline_check.csv (the baseline recomputed from the domain scores).
Seed   : 20 (numpy default_rng).  Run: python3 52_composite_rank_robustness.py
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import kendalltau

ROOT = "/path/to/six-node-ffl-breast-cancer"
OUT = ROOT + "/results/prioritisation_and_stroma"
DOMAINS = ['T_topology', 'E_evidence', 'D_expression', 'R_replication', 'M_multiomic', 'C_clinical', 'F_functional']
MIN_DOMAINS = 3
TOP = 10
N_DRAWS = 10000
CLASSES = ['TF', 'miRNA', 'Gene']

df = pd.read_csv(os.path.join(OUT, 'node_prioritisation_full.csv'))
X = df[DOMAINS].to_numpy(float)
avail = ~np.isnan(X)
X0 = np.where(avail, X, 0.0)
types = df['type'].to_numpy()
names = df['name'].to_numpy()


def composite(weights, keep=None):
    """Weighted mean over the available (and kept) domains; NaN when fewer than MIN_DOMAINS remain."""
    a = avail if keep is None else avail & keep[None, :]
    w = np.asarray(weights, float)[None, :] * a
    den = w.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        s = (X0 * w).sum(axis=1) / den
    s[a.sum(axis=1) < MIN_DOMAINS] = np.nan
    return s


def top_set(score):
    """Boolean mask of the TOP highest-scoring nodes within each class."""
    mask = np.zeros(len(score), bool)
    for c in CLASSES:
        idx = np.where((types == c) & ~np.isnan(score))[0]
        order = idx[np.argsort(-score[idx], kind='stable')]
        mask[order[:TOP]] = True
    return mask


def compare(score, base_score, base_top):
    m = ~np.isnan(score) & ~np.isnan(base_score)
    tau = kendalltau(score[m], base_score[m])[0]
    t = top_set(score)
    return dict(n_ranked=int((~np.isnan(score)).sum()), kendall_tau=tau, top30_shared=int((t & base_top).sum()),
                top30_jaccard=(t & base_top).sum() / (t | base_top).sum()), t


base = composite(np.ones(len(DOMAINS)))
base_top = top_set(base)
chk = pd.DataFrame(dict(name=names, type=types, composite_published=df['priority'], composite_recomputed=base))
chk['abs_difference'] = (chk.composite_published - chk.composite_recomputed).abs()
chk.to_csv(os.path.join(OUT, 'composite_rank_robustness_baseline_check.csv'), index=False)
assert np.nanmax(chk.abs_difference) < 1e-9, 'baseline composite differs from the published priority'
published_top = pd.read_csv(os.path.join(OUT, 'top10_per_class.csv'))
assert set(names[base_top]) == set(published_top.name), 'baseline top-ten lists differ from the published lists'

rows, stay = [], {}
t_count = np.zeros(len(df))

# A: leave one domain out
for j, d in enumerate(DOMAINS):
    keep = np.ones(len(DOMAINS), bool)
    keep[j] = False
    r, t = compare(composite(np.ones(len(DOMAINS)), keep), base, base_top)
    rows.append(dict(scheme='leave_one_domain_out', setting=d, **r))
    stay['without_' + d] = t

# C: median of the available domains
med = np.nanmedian(np.where(avail, X, np.nan), axis=1)
med[avail.sum(axis=1) < MIN_DOMAINS] = np.nan
r, t = compare(med, base, base_top)
rows.append(dict(scheme='median_of_domains', setting='median', **r))
stay['median'] = t

# B: random weights
rng = np.random.default_rng(20)
taus, shared, jac, ranks = [], [], [], np.full((N_DRAWS, len(df)), np.nan)
for i in range(N_DRAWS):
    w = rng.dirichlet(np.ones(len(DOMAINS)))
    s = composite(w)
    r, t = compare(s, base, base_top)
    taus.append(r['kendall_tau']); shared.append(r['top30_shared']); jac.append(r['top30_jaccard'])
    t_count += t
    for c in CLASSES:
        idx = np.where((types == c) & ~np.isnan(s))[0]
        order = idx[np.argsort(-s[idx], kind='stable')]
        ranks[i, order] = np.arange(1, len(order) + 1)
taus, shared, jac = map(np.array, (taus, shared, jac))
rows.append(dict(scheme='random_weights', setting=f'{N_DRAWS} Dirichlet(1,...,1) draws, median',
                 n_ranked=int((~np.isnan(base)).sum()), kendall_tau=np.median(taus),
                 top30_shared=np.median(shared), top30_jaccard=np.median(jac)))
rows.append(dict(scheme='random_weights', setting='2.5th percentile', n_ranked=int((~np.isnan(base)).sum()),
                 kendall_tau=np.percentile(taus, 2.5), top30_shared=np.percentile(shared, 2.5),
                 top30_jaccard=np.percentile(jac, 2.5)))
rows.append(dict(scheme='random_weights', setting='97.5th percentile', n_ranked=int((~np.isnan(base)).sum()),
                 kendall_tau=np.percentile(taus, 97.5), top30_shared=np.percentile(shared, 97.5),
                 top30_jaccard=np.percentile(jac, 97.5)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'composite_rank_robustness_summary.csv'), index=False)

# per-node table: the thirty baseline nodes plus every node that reached the top ten of its class in some draw
frac_top = t_count / N_DRAWS
nodes = pd.DataFrame(dict(name=names, type=types, composite=base,
                          baseline_rank_in_class=pd.Series(base).groupby(types).rank(ascending=False, method='first').to_numpy(),
                          in_baseline_top30=base_top, fraction_of_random_weights_in_top10_of_class=frac_top))
for k, t in stay.items():
    nodes['in_top10_' + k] = t
nodes['loo_domains_in_top10'] = nodes[[c for c in nodes if c.startswith('in_top10_without_')]].sum(axis=1)
nodes['median_rank_random_weights'] = np.nanmedian(ranks, axis=0)
nodes['rank_2.5pct_random_weights'] = np.nanpercentile(ranks, 2.5, axis=0)
nodes['rank_97.5pct_random_weights'] = np.nanpercentile(ranks, 97.5, axis=0)
keep_rows = nodes[nodes.in_baseline_top30 | (nodes.fraction_of_random_weights_in_top10_of_class > 0)]
keep_rows.sort_values(['type', 'baseline_rank_in_class']).to_csv(
    os.path.join(OUT, 'composite_rank_robustness_nodes.csv'), index=False)

print(pd.DataFrame(rows).round(3).to_string(index=False))
b = nodes[nodes.in_baseline_top30].sort_values(['type', 'baseline_rank_in_class'])
print(b[['name', 'type', 'baseline_rank_in_class', 'loo_domains_in_top10', 'in_top10_median',
         'fraction_of_random_weights_in_top10_of_class']].round(3).to_string(index=False))
