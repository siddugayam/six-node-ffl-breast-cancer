#!/usr/bin/env python3
"""(A) ARCHITECTURE of the canonical miRNA-TF-gene network.

  * degree distributions (exported for Clauset-Shalizi-Newman fitting in R/poweRlaw)
  * bow-tie decomposition (Broder et al. 2000 style) on the strongly connected core
  * strongly connected components: number and size distribution
  * flow hierarchy (Luo & Magee 2011, as implemented in networkx)
  * global reaching centrality (Mones, Vicsek & Vicsek 2012)
  * trophic levels / trophic incoherence (MacKay, Johnson & Sansom 2020)
  * Ravasz-Barabasi hierarchical-modularity test  C(k) ~ k^-1
  * motif-level clustering: FFL closure of regulatory triples, vs degree-preserving null

Refs: Ma'ayan et al. 2005 Science 309:1078 (bow-tie / regulatory motifs in signalling
networks); Berger & Iyengar 2009 Bioinformatics 25:2466-2472; Azeloglu & Iyengar 2015 CSH
Perspect Biol 7:a005934.
"""
import os, sys, json, csv, collections, random
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import (load_nodes, load_edges, RES, FIG, LOG, ffl_cores, ffl_unique_count)

random.seed(1); np.random.seed(1)
OUT = {}

nt = load_nodes()
E = load_edges()
P = [(e['source'], e['target']) for e in E]
G = nx.DiGraph()
G.add_nodes_from(nt.keys())
G.add_edges_from(P)
print(f"DiGraph: N={G.number_of_nodes()} M={G.number_of_edges()}")
OUT['N'] = G.number_of_nodes(); OUT['M'] = G.number_of_edges()

U = G.to_undirected()          # simple undirected projection
OUT['M_undirected'] = U.number_of_edges()

# ------------------------------------------------------------------ 1. degrees
deg = dict(G.degree()); ind = dict(G.in_degree()); outd = dict(G.out_degree())
with open(f'{RES}/systems_degree_table.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['node', 'type', 'degree_total', 'degree_in', 'degree_out',
                                    'degree_undirected'])
    for n in sorted(G.nodes()):
        w.writerow([n, nt.get(n, 'Gene'), deg[n], ind[n], outd[n], U.degree(n)])
print('wrote systems_degree_table.csv')
for lab, d in [('total', deg), ('in', ind), ('out', outd)]:
    v = np.array(list(d.values()))
    OUT[f'deg_{lab}_mean'] = float(v.mean()); OUT[f'deg_{lab}_max'] = int(v.max())
    OUT[f'deg_{lab}_n_zero'] = int((v == 0).sum())

# ------------------------------------------------------------------ 2. SCC / bow-tie
sccs = sorted(nx.strongly_connected_components(G), key=len, reverse=True)
OUT['n_SCC'] = len(sccs)
OUT['SCC_sizes_top10'] = [len(s) for s in sccs[:10]]
OUT['n_SCC_size1'] = sum(1 for s in sccs if len(s) == 1)
LSCC = set(sccs[0]); OUT['LSCC_size'] = len(LSCC)

lscc_rep = next(iter(LSCC))
anc = nx.ancestors(G, lscc_rep); desc = nx.descendants(G, lscc_rep)
IN = anc - LSCC
OUT_ = desc - LSCC
rest = set(G.nodes()) - LSCC - IN - OUT_
# tubes: reachable from IN and reaching OUT without passing through LSCC
tubes = set(); in_tend = set(); out_tend = set(); disc = set()
Gm = G.copy(); Gm.remove_nodes_from(LSCC)
from_IN = set()
for s in IN:
    if s in Gm:
        from_IN |= nx.descendants(Gm, s)
to_OUT = set()
Grev = Gm.reverse(copy=False)
for t in OUT_:
    if t in Grev:
        to_OUT |= nx.descendants(Grev, t)
for n in rest:
    if n in from_IN and n in to_OUT:
        tubes.add(n)
    elif n in from_IN:
        out_tend.add(n)      # tendril hanging off IN
    elif n in to_OUT:
        in_tend.add(n)       # tendril feeding OUT
    else:
        disc.add(n)
OUT['bowtie'] = dict(IN=len(IN), CORE=len(LSCC), OUT=len(OUT_), TUBES=len(tubes),
                     TENDRILS_from_IN=len(out_tend), TENDRILS_to_OUT=len(in_tend),
                     OTHER=len(disc))
print('bow-tie:', OUT['bowtie'])
with open(f'{RES}/systems_bowtie_membership.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['node', 'type', 'bowtie_component', 'degree_total'])
    for n in sorted(G.nodes()):
        c = ('CORE' if n in LSCC else 'IN' if n in IN else 'OUT' if n in OUT_ else
             'TUBE' if n in tubes else 'TENDRIL_from_IN' if n in out_tend else
             'TENDRIL_to_OUT' if n in in_tend else 'OTHER')
        w.writerow([n, nt.get(n, 'Gene'), c, deg[n]])
# node-type composition of each bow-tie sector
comp = collections.defaultdict(collections.Counter)
for n in G.nodes():
    c = ('CORE' if n in LSCC else 'IN' if n in IN else 'OUT' if n in OUT_ else 'PERIPHERY')
    comp[c][nt.get(n, 'Gene')] += 1
OUT['bowtie_types'] = {k: dict(v) for k, v in comp.items()}
print('bow-tie composition:', OUT['bowtie_types'])

# fraction of edges inside the core
e_in_core = sum(1 for u, v in G.edges() if u in LSCC and v in LSCC)
OUT['edges_within_core'] = e_in_core
OUT['edges_within_core_frac'] = e_in_core / G.number_of_edges()

# ------------------------------------------------------------------ 3. hierarchy
OUT['flow_hierarchy'] = float(nx.flow_hierarchy(G))
OUT['global_reaching_centrality'] = float(nx.global_reaching_centrality(G))
# reciprocity
OUT['reciprocity'] = float(nx.reciprocity(G))

# ---- trophic levels (MacKay, Johnson & Sansom 2020 Roy Soc Open Sci)
# solve  (Lambda) h = v  on the largest weakly connected component
wcc = max(nx.weakly_connected_components(G), key=len)
OUT['n_WCC'] = nx.number_weakly_connected_components(G)
OUT['LWCC_size'] = len(wcc)
Gw = G.subgraph(wcc)
nodes_w = sorted(Gw.nodes()); idx = {n: i for i, n in enumerate(nodes_w)}
nW = len(nodes_w)
A = np.zeros((nW, nW))
for u, v in Gw.edges():
    A[idx[u], idx[v]] = 1.0
kin = A.sum(0); kout = A.sum(1)
u_vec = kin + kout
v_vec = kin - kout
Lam = np.diag(u_vec) - (A + A.T)
# Lambda is singular (constant vector); solve least-squares then centre
h, *_ = np.linalg.lstsq(Lam, v_vec, rcond=None)
h = h - h.min()
num = 0.0
for u_, v_ in Gw.edges():
    num += (h[idx[v_]] - h[idx[u_]] - 1) ** 2
F0 = np.sqrt(num / Gw.number_of_edges())      # trophic incoherence parameter
OUT['trophic_incoherence_F0'] = float(F0)
OUT['trophic_level_range'] = [float(h.min()), float(h.max())]
with open(f'{RES}/systems_trophic_levels.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['node', 'type', 'trophic_level', 'degree_in', 'degree_out'])
    for n in nodes_w:
        w.writerow([n, nt.get(n, 'Gene'), round(float(h[idx[n]]), 6), ind[n], outd[n]])
# mean trophic level by node class
tl_by = collections.defaultdict(list)
for n in nodes_w:
    tl_by[nt.get(n, 'Gene')].append(float(h[idx[n]]))
OUT['trophic_level_mean_by_type'] = {k: float(np.mean(v)) for k, v in tl_by.items()}
print('trophic incoherence F0 =', round(F0, 4), OUT['trophic_level_mean_by_type'])

# ------------------------------------------------------------------ 4. clustering
OUT['transitivity_undirected'] = float(nx.transitivity(U))
OUT['avg_clustering_undirected'] = float(nx.average_clustering(U))
OUT['avg_clustering_directed'] = float(nx.average_clustering(G))
# Ravasz-Barabasi C(k) ~ k^-1 test on undirected projection
cl = nx.clustering(U)
by_k = collections.defaultdict(list)
for n, c in cl.items():
    k = U.degree(n)
    if k >= 2:
        by_k[k].append(c)
ks = sorted(by_k)
Ck = np.array([np.mean(by_k[k]) for k in ks]); kk = np.array(ks, float)
m = Ck > 0
if m.sum() >= 5:
    b, a = np.polyfit(np.log10(kk[m]), np.log10(Ck[m]), 1)
    r = np.corrcoef(np.log10(kk[m]), np.log10(Ck[m]))[0, 1]
    OUT['Ck_scaling_exponent'] = float(b)
    OUT['Ck_scaling_r'] = float(r)
    OUT['Ck_n_bins'] = int(m.sum())
with open(f'{RES}/systems_Ck_scaling.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['k', 'mean_clustering', 'n_nodes'])
    for k in ks:
        w.writerow([k, round(float(np.mean(by_k[k])), 6), len(by_k[k])])

# ------------------------------------------------------------------ 5. motif clustering
Pset = set(P)
cores = ffl_cores(Pset, nt)
n_ffl = ffl_unique_count(cores)
OUT['ffl_cores_unique'] = n_ffl
OUT['ffl_cores_raw'] = len(cores)
# regulatory triples eligible to close into an FFL:
#   R->M with (miRNA,TF) or (TF,miRNA) typing, and R->T for some non-miRNA T
out_adj = collections.defaultdict(set)
for a, b in Pset:
    out_adj[a].add(b)
open_paths = 0; closed_paths = 0
for R, tg in out_adj.items():
    tR = nt.get(R, 'Gene')
    for M in tg:
        tM = nt.get(M, 'Gene')
        if not ((tR == 'miRNA' and tM == 'TF') or (tR == 'TF' and tM == 'miRNA')):
            continue
        mt = out_adj.get(M, set())
        for T in mt:
            if T in (R, M) or nt.get(T, 'Gene') == 'miRNA':
                continue
            if T in tg:
                closed_paths += 1
            else:
                open_paths += 1
OUT['ffl_closure_observed'] = closed_paths / (closed_paths + open_paths)
OUT['ffl_closure_counts'] = [closed_paths, open_paths]

# degree-preserving, edge-type-preserving, reciprocity-preserving null for closure & FFL count
def rewire_null(pairs, edge_types, n_swaps_per_edge=20, rng=None):
    """Double-edge swap within edge_type class; keeps in/out degree per class.
    Mutual TF<->miRNA pairs are the entire source of reciprocity here, and swapping
    within class does not create/destroy them systematically, so reciprocity is
    monitored and reported rather than hard-constrained."""
    rng = rng or random
    byt = collections.defaultdict(list)
    for (a, b) in pairs:
        byt[edge_types[(a, b)]].append([a, b])
    cur = set(pairs)
    for t, lst in byt.items():
        if len(lst) < 4:
            continue
        target = n_swaps_per_edge * len(lst)
        done = 0; tries = 0
        while done < target and tries < 40 * target:
            tries += 1
            i, j = rng.randrange(len(lst)), rng.randrange(len(lst))
            if i == j:
                continue
            a, b = lst[i]; c, d = lst[j]
            if a == d or c == b or a == c or b == d:
                continue
            if (a, d) in cur or (c, b) in cur:
                continue
            cur.discard((a, b)); cur.discard((c, d))
            cur.add((a, d)); cur.add((c, b))
            lst[i] = [a, d]; lst[j] = [c, b]
            done += 1
    return cur

etypes = {(e['source'], e['target']): e['edge_type'] for e in E}
NNULL = int(os.environ.get('NNULL', '200'))
null_ffl = []; null_close = []; null_core = []; null_F0 = []; null_gr = []
rng = random.Random(7)
for i in range(NNULL):
    Q = rewire_null(P, etypes, 20, rng)
    null_ffl.append(ffl_unique_count(ffl_cores(Q, nt)))
    Gq = nx.DiGraph(); Gq.add_nodes_from(nt); Gq.add_edges_from(Q)
    s = max(nx.strongly_connected_components(Gq), key=len)
    null_core.append(len(s))
    null_gr.append(nx.global_reaching_centrality(Gq))
    if (i + 1) % 50 == 0:
        print('null', i + 1)
null_ffl = np.array(null_ffl); null_core = np.array(null_core); null_gr = np.array(null_gr)
OUT['null_n'] = NNULL
OUT['null_ffl_mean'] = float(null_ffl.mean()); OUT['null_ffl_sd'] = float(null_ffl.std(ddof=1))
OUT['ffl_zscore'] = float((n_ffl - null_ffl.mean()) / null_ffl.std(ddof=1))
OUT['ffl_p_emp'] = float((np.sum(null_ffl >= n_ffl) + 1) / (NNULL + 1))
OUT['null_core_mean'] = float(null_core.mean()); OUT['null_core_sd'] = float(null_core.std(ddof=1))
OUT['core_p_emp'] = float((np.sum(null_core >= len(LSCC)) + 1) / (NNULL + 1))
OUT['null_grc_mean'] = float(null_gr.mean()); OUT['null_grc_sd'] = float(null_gr.std(ddof=1))
OUT['grc_p_emp'] = float((np.sum(null_gr >= OUT['global_reaching_centrality']) + 1) / (NNULL + 1))
print('FFL z =', round(OUT['ffl_zscore'], 2), 'core obs/null =', len(LSCC), null_core.mean())

with open(f'{RES}/systems_architecture_summary.json', 'w') as fh:
    json.dump(OUT, fh, indent=2, sort_keys=True)
rows = []
def flat(prefix, d):
    for k, v in d.items():
        if isinstance(v, dict):
            flat(f'{prefix}{k}.', v)
        else:
            rows.append([f'{prefix}{k}', json.dumps(v) if isinstance(v, list) else v])
flat('', OUT)
with open(f'{RES}/systems_architecture_summary.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['metric', 'value']); w.writerows(sorted(rows))
np.save(f'{RES}/systems_null_ffl_counts.npy', null_ffl)
print(json.dumps(OUT, indent=2, sort_keys=True))
