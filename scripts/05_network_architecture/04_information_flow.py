#!/usr/bin/env python3
"""(C) INFORMATION FLOW through the regulatory network.

Iyengar's programme treats a signalling network as a device that routes information; the
nodes that matter are those through which flow is forced, which need not be the degree hubs
(Ma'ayan et al. 2005 Science 309:1078; Berger & Iyengar 2009 Bioinformatics 25:2466-2472).

Measures computed
  1. Current-flow (random-walk) betweenness   - Newman 2005 Soc Networks 27:39;
     Brandes & Fleischer 2005.  Undirected projection (connected: single WCC).
  2. Information centrality  - Stephenson & Zelen 1989; networkx
     current_flow_closeness_centrality.
  3. Directed random-walk signal flux: expected number of visits to v by a damped random
     walk started at s, V = (I - (1-lambda) P)^-1 with lambda = 0.15.  Row-stochastic P;
     the 206 out-degree-0 nodes are absorbing.  Convergence: the Neumann series converges
     because the spectral radius of (1-lambda)P is <= 1-lambda = 0.85 < 1; the closed-form
     inverse is used, so the result is exact.
  4. Shortest-path betweenness (directed) for reference - what the manuscript's hub table used.

Then: do these rank nodes differently from degree?
"""
import os, sys, csv, json, collections
import numpy as np
import networkx as nx
from scipy.stats import spearmanr, kendalltau

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES, FIG, ffl_cores

nt = load_nodes()
E = load_edges()
P = [(e['source'], e['target']) for e in E]
G = nx.DiGraph(); G.add_nodes_from(nt); G.add_edges_from(P)
U = G.to_undirected()
assert nx.is_connected(U), 'undirected projection not connected'
NODES = sorted(G.nodes()); idx = {n: i for i, n in enumerate(NODES)}; N = len(NODES)
print(f'N={N} directed M={G.number_of_edges()} undirected M={U.number_of_edges()}')

# 1-2 undirected current-flow measures
print('current-flow betweenness ...', flush=True)
cfb = nx.current_flow_betweenness_centrality(U, normalized=True, solver='full')
print('information centrality ...', flush=True)
infc = nx.current_flow_closeness_centrality(U, solver='full')
print('shortest-path betweenness (directed) ...', flush=True)
spb = nx.betweenness_centrality(G, normalized=True)
spb_u = nx.betweenness_centrality(U, normalized=True)

# 3 directed damped random-walk flux
LAM = 0.15
A = np.zeros((N, N))
for u, v in P:
    A[idx[u], idx[v]] = 1.0
rs = A.sum(1)
Pm = np.zeros_like(A)
nz = rs > 0
Pm[nz] = A[nz] / rs[nz, None]
V = np.linalg.inv(np.eye(N) - (1 - LAM) * Pm)      # V[s,v] = expected visits to v from s
spec = np.max(np.abs(np.linalg.eigvals((1 - LAM) * Pm)))
print(f'spectral radius of (1-lambda)P = {spec:.4f} (<1 required for convergence)')
np.fill_diagonal_ = None
Vt = V.copy()
np.fill_diagonal(Vt, 0.0)                          # exclude the trivial s==v term
flux_in = Vt.sum(0)                                # total signal arriving at v from all sources
flux_out = Vt.sum(1)                               # total signal v can deliver
# throughput betweenness: expected visits to v on walks that do not start or end at v,
# weighted by v's probability of forwarding rather than absorbing
forward = np.where(rs > 0, 1.0 - LAM, 0.0)
rw_between = flux_in * forward

deg = dict(G.degree()); ind = dict(G.in_degree()); outd = dict(G.out_degree())
cores = ffl_cores(set(P), nt)
part = collections.Counter()
for R, M, T, c in cores:
    part[R] += 1; part[M] += 1; part[T] += 1
ctrl = {}
with open(f'{RES}/systems_controllability_nodes.csv') as fh:
    for r in csv.DictReader(fh):
        ctrl[r['node']] = r

rows = []
for n in NODES:
    i = idx[n]
    rows.append(dict(node=n, type=nt[n], degree_total=deg[n], degree_in=ind[n],
                     degree_out=outd[n], ffl_participation=part[n],
                     current_flow_betweenness=cfb[n], information_centrality=infc[n],
                     rw_flux_in=float(flux_in[i]), rw_flux_out=float(flux_out[i]),
                     rw_throughput_betweenness=float(rw_between[i]),
                     shortest_path_betweenness_directed=spb[n],
                     shortest_path_betweenness_undirected=spb_u[n],
                     is_driver=ctrl.get(n, {}).get('is_driver_one_config', ''),
                     deletion_class=ctrl.get(n, {}).get('deletion_class', '')))
flds = list(rows[0].keys())
with open(f'{RES}/systems_information_flow.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=flds); w.writeheader(); w.writerows(rows)
print('wrote systems_information_flow.csv')

# ---------------- comparison of rankings
meas = ['degree_total', 'current_flow_betweenness', 'information_centrality',
        'rw_throughput_betweenness', 'rw_flux_in', 'shortest_path_betweenness_directed',
        'ffl_participation']
M_ = {m: np.array([r[m] for r in rows], float) for m in meas}
cmp_rows = []
for m in meas:
    if m == 'degree_total':
        continue
    s = spearmanr(M_['degree_total'], M_[m]); k = kendalltau(M_['degree_total'], M_[m])
    top_d = set(np.argsort(-M_['degree_total'])[:20]); top_m = set(np.argsort(-M_[m])[:20])
    j = len(top_d & top_m) / len(top_d | top_m)
    cmp_rows.append(dict(measure=m, spearman_vs_degree=float(s.statistic),
                         spearman_p=float(s.pvalue), kendall_vs_degree=float(k.statistic),
                         top20_overlap=len(top_d & top_m), top20_jaccard=round(j, 3)))
    print(f'{m:<38} rho_vs_degree={s.statistic:+.3f} top20 overlap={len(top_d&top_m)}/20 J={j:.2f}')
with open(f'{RES}/systems_flow_vs_degree_comparison.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(cmp_rows[0].keys())); w.writeheader(); w.writerows(cmp_rows)

# top-20 tables side by side
top = {}
for m in meas:
    order = np.argsort(-M_[m])[:20]
    top[m] = [NODES[i] for i in order]
with open(f'{RES}/systems_top20_by_measure.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['rank'] + meas)
    for r in range(20):
        w.writerow([r + 1] + [top[m][r] for m in meas])
print('\nTop 10 by current-flow betweenness:', top['current_flow_betweenness'][:10])
print('Top 10 by degree            :', top['degree_total'][:10])
print('Top 10 by information centr.:', top['information_centrality'][:10])
print('Top 10 by RW throughput     :', top['rw_throughput_betweenness'][:10])

# "bottlenecks": high flow, unremarkable degree
rk = {m: {NODES[i]: int(r) for r, i in enumerate(np.argsort(-M_[m]))} for m in meas}
bott = sorted(NODES, key=lambda n: rk['degree_total'][n] - rk['current_flow_betweenness'][n],
              reverse=True)[:20]
with open(f'{RES}/systems_flow_bottlenecks.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['node', 'type', 'degree_total', 'rank_degree', 'rank_current_flow_betweenness',
                'rank_gain', 'current_flow_betweenness'])
    for n in bott:
        w.writerow([n, nt[n], deg[n], rk['degree_total'][n] + 1,
                    rk['current_flow_betweenness'][n] + 1,
                    rk['degree_total'][n] - rk['current_flow_betweenness'][n], cfb[n]])
print('\nFlow bottlenecks (rank gain over degree):',
      [(n, deg[n], rk['degree_total'][n] + 1, rk['current_flow_betweenness'][n] + 1) for n in bott[:10]])

# where do the focal nodes sit?
FOCUS = ['hsa-miR-29a', 'hsa-miR-29b', 'hsa-miR-29c', 'COL1A1', 'COL3A1', 'SP1', 'RELA',
         'NFKB1', 'MYC', 'STAT3', 'HIF1A', 'TP53', 'VEGFA']
frows = []
for n in FOCUS:
    if n not in idx:
        print('MISSING', n); continue
    frows.append([n, nt[n], deg[n], rk['degree_total'][n] + 1,
                  rk['current_flow_betweenness'][n] + 1, rk['information_centrality'][n] + 1,
                  rk['rw_throughput_betweenness'][n] + 1, rk['ffl_participation'][n] + 1])
    print(f'{n:<14} deg={deg[n]:<4} rank_deg={rk["degree_total"][n]+1:<4} '
          f'rank_CFB={rk["current_flow_betweenness"][n]+1:<4} '
          f'rank_IC={rk["information_centrality"][n]+1:<4} '
          f'rank_RW={rk["rw_throughput_betweenness"][n]+1:<4}')
with open(f'{RES}/systems_focus_node_ranks.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['node', 'type', 'degree', 'rank_degree', 'rank_current_flow_betweenness',
                'rank_information_centrality', 'rank_rw_throughput', 'rank_ffl_participation'])
    w.writerows(frows)

json.dump(dict(spectral_radius=float(spec), lambda_damping=LAM,
               n_nodes=N, n_edges=G.number_of_edges()),
          open(f'{RES}/systems_information_flow_meta.json', 'w'), indent=2)
