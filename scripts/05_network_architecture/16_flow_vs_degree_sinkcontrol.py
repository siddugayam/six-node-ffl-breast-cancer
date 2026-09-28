#!/usr/bin/env python3
"""(C cont.) The flow-vs-degree comparison, controlled for the structural sinks.

206 of the 587 nodes (all genes bar COL1A1) have out-degree 0 because the assembled network
contains essentially no gene->anything edges. Such a node can never forward a signal, so its
directed random-walk throughput is exactly 0 and it is tied with the other 210 zero-throughput
nodes. Any statement that "information flow ranks nodes differently from degree" has to be
re-checked on the 381 nodes that can actually forward signal, otherwise the effect is just the
sink block.

Two further caveats made explicit here:
  * current-flow betweenness and information centrality are computed on the UNDIRECTED
    projection (Newman 2005 Soc Networks 27:39-54; Stephenson & Zelen 1989), so they are
    direction-blind and will credit a pure sink with routing capacity it does not have.
  * the "bottleneck" list (high flow, unremarkable degree) is therefore recomputed restricted
    to nodes with out-degree > 0.

Ma'ayan et al. 2005 Science 309:1078-1083; Berger & Iyengar 2009 Bioinformatics 25:2466-2472.
"""
import os, sys, csv, json
import numpy as np
from scipy.stats import spearmanr, kendalltau

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import RES

rows = list(csv.DictReader(open(f'{RES}/systems_information_flow.csv')))
node = [r['node'] for r in rows]
typ = [r['type'] for r in rows]
V = {k: np.array([float(r[k]) for r in rows]) for k in
     ('degree_total', 'degree_out', 'current_flow_betweenness', 'information_centrality',
      'rw_throughput_betweenness', 'rw_flux_in', 'shortest_path_betweenness_directed',
      'ffl_participation')}
od = V['degree_out']
can_forward = od > 0
print(f'{len(rows)} nodes; {int((~can_forward).sum())} with out-degree 0 (structural sinks); '
      f'{int(can_forward.sum())} can forward signal')
print(f'nodes with directed RW throughput exactly 0: {int((V["rw_throughput_betweenness"]==0).sum())}')

meas = ['current_flow_betweenness', 'information_centrality', 'rw_throughput_betweenness',
        'rw_flux_in', 'shortest_path_betweenness_directed', 'ffl_participation']
out = []
for m in meas:
    for lab, msk in [('all nodes', np.ones(len(rows), bool)),
                     ('out-degree > 0 only', can_forward)]:
        d = V['degree_total'][msk]; v = V[m][msk]
        s = spearmanr(d, v); k = kendalltau(d, v)
        t1 = set(np.argsort(-d)[:20]); t2 = set(np.argsort(-v)[:20])
        out.append(dict(measure=m, node_set=lab, n=int(msk.sum()),
                        spearman_vs_degree=round(float(s.statistic), 4),
                        spearman_p=float(s.pvalue),
                        kendall_vs_degree=round(float(k.statistic), 4),
                        top20_overlap=len(t1 & t2),
                        top20_jaccard=round(len(t1 & t2) / len(t1 | t2), 3)))
        print(f'{m:<36} {lab:<20} n={int(msk.sum()):<4} rho={s.statistic:+.3f} '
              f'top20 overlap={len(t1 & t2)}/20')
with open(f'{RES}/systems_flow_vs_degree_sinkcontrol.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)

# ---- bottlenecks among nodes that can actually forward signal
idxf = np.where(can_forward)[0]
rk_deg = {}; rk_cfb = {}; rk_rw = {}
for name, key, store in [('deg', 'degree_total', rk_deg),
                         ('cfb', 'current_flow_betweenness', rk_cfb),
                         ('rw', 'rw_throughput_betweenness', rk_rw)]:
    order = idxf[np.argsort(-V[key][idxf])]
    for r, i in enumerate(order):
        store[i] = r + 1
gain = sorted(idxf, key=lambda i: rk_deg[i] - rk_rw[i], reverse=True)
brows = []
print('\nnodes that carry far more directed signal than their degree predicts '
      '(out-degree > 0 only):')
for i in gain[:20]:
    brows.append(dict(node=node[i], type=typ[i], degree_total=int(V['degree_total'][i]),
                      degree_out=int(od[i]), rank_degree=rk_deg[i],
                      rank_rw_throughput=rk_rw[i],
                      rank_current_flow_betweenness=rk_cfb[i],
                      rank_gain_degree_minus_rw=rk_deg[i] - rk_rw[i],
                      rw_throughput=float(V['rw_throughput_betweenness'][i])))
    print(f"  {node[i]:<14} {typ[i]:<6} deg={int(V['degree_total'][i]):<4} "
          f"rank_deg={rk_deg[i]:<4} rank_RW={rk_rw[i]:<4} gain={rk_deg[i]-rk_rw[i]:+d}")
with open(f'{RES}/systems_flow_bottlenecks_forwarding.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(brows[0].keys())); w.writeheader(); w.writerows(brows)

# ---- degree hubs that cannot route anything
deg_order = np.argsort(-V['degree_total'])[:30]
dead = [(node[i], typ[i], int(V['degree_total'][i])) for i in deg_order if od[i] == 0]
print(f'\ntop-30 degree hubs with out-degree 0 (cannot forward any signal): '
      f'{len(dead)} -> {[d[0] for d in dead]}')
json.dump(dict(n_nodes=len(rows), n_sinks=int((~can_forward).sum()),
               n_zero_rw_throughput=int((V['rw_throughput_betweenness'] == 0).sum()),
               top30_degree_hubs_that_are_sinks=[d[0] for d in dead]),
          open(f'{RES}/systems_flow_vs_degree_sinkcontrol.json', 'w'), indent=2)
print('\nwrote systems_flow_vs_degree_sinkcontrol.{csv,json} and '
      'systems_flow_bottlenecks_forwarding.csv')
