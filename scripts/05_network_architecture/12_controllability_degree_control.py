#!/usr/bin/env python3
"""(B cont.) Is the driver-node depletion among FFL hubs anything more than a degree effect?

Driver frequency correlates with degree at Spearman -0.54, and FFL hubs are high-degree by
construction, so the raw Fisher test is confounded.  Two controls:
  1. logistic regression  driver ~ log(degree+1) + node_type + is_FFL_hub
  2. degree-matched resampling: for each FFL hub draw a non-hub of the closest degree
     (10,000 resamples) and compare driver rates

Also records the structural constraints that make several of the topology results trivial,
so they are not over-read.
"""
import os, sys, csv, json, collections
import numpy as np
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES

rng = np.random.default_rng(11)
nt = load_nodes(); E = load_edges()
ctrl = list(csv.DictReader(open(f'{RES}/systems_controllability_nodes.csv')))
NODES = [r['node'] for r in ctrl]
deg = np.array([float(r['degree_total']) for r in ctrl])
outd = np.array([float(r['degree_out']) for r in ctrl])
drv = np.array([r['is_driver_one_config'] == 'True' for r in ctrl])
hub = np.array([r['is_FFL_hub'] == 'True' for r in ctrl])
ffl = np.array([float(r['ffl_participation']) for r in ctrl])
typ = [r['type'] for r in ctrl]

out = {}
# ---- 1. logistic regression
import statsmodels.api as sm
X = np.column_stack([np.log1p(deg),
                     [1.0 if t == 'TF' else 0.0 for t in typ],
                     [1.0 if t == 'Gene' else 0.0 for t in typ],
                     hub.astype(float)])
X = sm.add_constant(X)
m = sm.Logit(drv.astype(float), X).fit(disp=0)
names = ['const', 'log(degree+1)', 'is_TF', 'is_Gene', 'is_FFL_hub']
print(m.summary2().tables[1].to_string())
out['logistic_driver'] = {n: dict(coef=float(m.params[i]), se=float(m.bse[i]),
                                  z=float(m.tvalues[i]), p=float(m.pvalues[i]))
                          for i, n in enumerate(names)}
print(f"\nis_FFL_hub coefficient after adjusting for degree and node type: "
      f"{m.params[4]:+.3f} (p = {m.pvalues[4]:.3g})")

# also with FFL participation as a continuous predictor instead of the hub flag
X2 = sm.add_constant(np.column_stack([np.log1p(deg),
                                      [1.0 if t == 'TF' else 0.0 for t in typ],
                                      [1.0 if t == 'Gene' else 0.0 for t in typ],
                                      np.log1p(ffl)]))
m2 = sm.Logit(drv.astype(float), X2).fit(disp=0)
out['logistic_driver_continuous_ffl'] = dict(
    coef=float(m2.params[4]), se=float(m2.bse[4]), z=float(m2.tvalues[4]),
    p=float(m2.pvalues[4]))
print(f"log(FFL participation+1) coefficient after the same adjustment: "
      f"{m2.params[4]:+.3f} (p = {m2.pvalues[4]:.3g})")

# ---- 2. degree-matched resampling
hub_i = np.where(hub)[0]; non_i = np.where(~hub)[0]
obs = drv[hub_i].mean()
NB = 10000
nullr = np.empty(NB)
for b in range(NB):
    pick = []
    for i in hub_i:
        d = np.abs(deg[non_i] - deg[i])
        cand = non_i[d <= max(2.0, d.min())]
        pick.append(rng.choice(cand))
    nullr[b] = drv[pick].mean()
p = (np.sum(nullr <= obs) + 1) / (NB + 1)
print(f'\ndegree-matched null: FFL hubs {obs:.3f} driver rate vs matched non-hubs '
      f'{nullr.mean():.3f} +/- {nullr.std(ddof=1):.3f}, one-sided p = {p:.4g} ({NB} resamples)')
out['degree_matched'] = dict(observed=float(obs), null_mean=float(nullr.mean()),
                             null_sd=float(nullr.std(ddof=1)), p_one_sided=float(p),
                             n_resamples=NB, n_hubs=int(hub.sum()))

# ---- 3. structural constraints that make some results trivial
od = collections.Counter(); idg = collections.Counter()
for e in E:
    od[e['source']] += 1; idg[e['target']] += 1
sinks = [n for n in nt if od[n] == 0]
sources = [n for n in nt if idg[n] == 0]
cons = dict(
    n_out_degree_zero=len(sinks),
    out_degree_zero_by_type=dict(collections.Counter(nt[n] for n in sinks)),
    n_in_degree_zero=len(sources),
    in_degree_zero_by_type=dict(collections.Counter(nt[n] for n in sources)),
    genes_with_any_out_edge=[n for n in nt if nt[n] == 'Gene' and od[n] > 0],
    edge_class_counts={f"{k[0]}->{k[1]}": v for k, v in
                       collections.Counter((e['source_type'], e['target_type'])
                                           for e in E).items()},
    note=('206 of the 207 genes have out-degree 0 because the assembled network contains '
          'essentially no gene->anything edges (one STRING edge, COL1A1->COL3A1). The '
          'bow-tie OUT layer being made of genes, and genes never appearing in the strongly '
          'connected core, are therefore consequences of how the network was built, not '
          'discoveries about breast cancer regulation. Likewise the 5 "critical" driver '
          'nodes are exactly the 5 miRNAs with in-degree 0, which must be drivers in every '
          'control configuration.'))
out['structural_constraints'] = cons
print('\nSTRUCTURAL CONSTRAINTS:', json.dumps(cons, indent=2))

json.dump(out, open(f'{RES}/systems_controllability_degree_control.json', 'w'), indent=2)
rows = [dict(test='logistic: is_FFL_hub | degree + node type',
             estimate=round(out['logistic_driver']['is_FFL_hub']['coef'], 4),
             se=round(out['logistic_driver']['is_FFL_hub']['se'], 4),
             statistic=round(out['logistic_driver']['is_FFL_hub']['z'], 3),
             p=out['logistic_driver']['is_FFL_hub']['p'],
             detail='log-odds of being a driver node'),
        dict(test='logistic: log(FFL participation+1) | degree + node type',
             estimate=round(out['logistic_driver_continuous_ffl']['coef'], 4),
             se=round(out['logistic_driver_continuous_ffl']['se'], 4),
             statistic=round(out['logistic_driver_continuous_ffl']['z'], 3),
             p=out['logistic_driver_continuous_ffl']['p'],
             detail='log-odds of being a driver node'),
        dict(test='degree-matched resampling of FFL hubs',
             estimate=round(obs - nullr.mean(), 4), se=round(nullr.std(ddof=1), 4),
             statistic='', p=float(p),
             detail=f'driver rate {obs:.3f} vs degree-matched {nullr.mean():.3f}')]
with open(f'{RES}/systems_controllability_degree_control.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\nwrote systems_controllability_degree_control.{csv,json}')
