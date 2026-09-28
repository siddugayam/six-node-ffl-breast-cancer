#!/usr/bin/env python3
"""(A cont.) Is the "hierarchical regulatory structure" real, or an artefact of how the
network was assembled?

In the deposited network 206 of the 207 genes have out-degree 0 (the sole exception is the
single STRING edge COL1A1->COL3A1). Genes therefore CANNOT be in a strongly connected
component and MUST fall in the bow-tie OUT layer, whatever the biology. Any statement of the
form "TFs and miRNAs sit above genes in a regulatory hierarchy" is thus true by construction.

The non-trivial question is whether there is hierarchy AMONG THE REGULATORS, where edges can
run in both directions. This script repeats the architecture analysis on the induced
subnetwork of TF and miRNA nodes only, and reports, for both the full and the regulator-only
network:

  * SCC count and size distribution, bow-tie sectors
  * flow hierarchy (fraction of edges not on any cycle; Luo & Magee 2011)
  * global reaching centrality (Mones, Vicsek & Vicsek 2012 PLoS ONE 7:e33799)
  * trophic incoherence F0 (MacKay, Johnson & Sansom 2020 R Soc Open Sci 7:201138)
  * the depth of the condensation DAG = the number of genuine hierarchical layers
  * the same quantities on 200 degree-preserving curveball randomisations of the regulator
    subnetwork (Strona et al. 2014 Nat Commun 5:4114), so that "hierarchy" is judged against
    what the degree sequence alone produces.
"""
import os, sys, csv, json, collections
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES

rng = np.random.default_rng(20260909)
NRAND = int(sys.argv[1]) if len(sys.argv) > 1 else 200

nt = load_nodes()
E = load_edges()
P = [(e['source'], e['target']) for e in E]


def trophic_F0_and_levels(G):
    wcc = max(nx.weakly_connected_components(G), key=len)
    Gw = G.subgraph(wcc)
    ns = sorted(Gw.nodes()); ix = {n: i for i, n in enumerate(ns)}
    n = len(ns)
    A = np.zeros((n, n))
    for u, v in Gw.edges():
        A[ix[u], ix[v]] = 1.0
    kin = A.sum(0); kout = A.sum(1)
    Lam = np.diag(kin + kout) - (A + A.T)
    h, *_ = np.linalg.lstsq(Lam, kin - kout, rcond=None)
    h = h - h.min()
    num = sum((h[ix[v]] - h[ix[u]] - 1) ** 2 for u, v in Gw.edges())
    return float(np.sqrt(num / Gw.number_of_edges())), {n: float(h[ix[n]]) for n in ns}


def dag_depth(G):
    """Longest path (in nodes) through the condensation DAG = number of hierarchical layers."""
    C = nx.condensation(G)
    return int(nx.dag_longest_path_length(C)) + 1, C.number_of_nodes()


def describe(G, label):
    sccs = sorted(nx.strongly_connected_components(G), key=len, reverse=True)
    L = sccs[0]
    rep = next(iter(L))
    IN = nx.ancestors(G, rep) - L
    OU = nx.descendants(G, rep) - L
    other = set(G.nodes()) - L - IN - OU
    F0, _ = trophic_F0_and_levels(G)
    depth, ncomp = dag_depth(G)
    d = dict(network=label, N=G.number_of_nodes(), M=G.number_of_edges(),
             n_SCC=len(sccs), largest_SCC=len(L),
             frac_nodes_in_largest_SCC=len(L) / G.number_of_nodes(),
             n_SCC_singleton=sum(1 for s in sccs if len(s) == 1),
             bowtie_IN=len(IN), bowtie_CORE=len(L), bowtie_OUT=len(OU),
             bowtie_other=len(other),
             flow_hierarchy=float(nx.flow_hierarchy(G)),
             global_reaching_centrality=float(nx.global_reaching_centrality(G)),
             trophic_incoherence_F0=F0,
             condensation_depth_layers=depth,
             reciprocity=float(nx.reciprocity(G)),
             transitivity=float(nx.transitivity(G.to_undirected())))
    return d


G_full = nx.DiGraph(); G_full.add_nodes_from(nt); G_full.add_edges_from(P)
REG = {n for n, t in nt.items() if t in ('TF', 'miRNA')}
P_reg = [(u, v) for u, v in P if u in REG and v in REG]
G_reg = nx.DiGraph(); G_reg.add_nodes_from(sorted(REG)); G_reg.add_edges_from(P_reg)

print(f'full network      : {G_full.number_of_nodes()} nodes, {G_full.number_of_edges()} edges')
print(f'regulator-only    : {G_reg.number_of_nodes()} nodes '
      f'({sum(1 for n in REG if nt[n]=="TF")} TF + {sum(1 for n in REG if nt[n]=="miRNA")} miRNA), '
      f'{G_reg.number_of_edges()} edges')

od = collections.Counter(); idg = collections.Counter()
for u, v in P:
    od[u] += 1; idg[v] += 1
genes = [n for n in nt if nt[n] == 'Gene']
print(f'genes with out-degree 0: {sum(1 for g in genes if od[g] == 0)}/{len(genes)}')

rows = [describe(G_full, 'full (587 nodes)'), describe(G_reg, 'regulators only (TF+miRNA)')]
for r in rows:
    print('\n' + r['network'])
    for k, v in r.items():
        if k != 'network':
            print(f'   {k:<32} {v}')

# ---------------------------------------------------------------- null for the regulators
byc = collections.defaultdict(list)
for e in E:
    if e['source'] in REG and e['target'] in REG:
        byc[e['edge_type']].append((e['source'], e['target']))


def curveball(edges):
    adj = collections.defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
    keys = list(adj)
    if len(keys) < 2:
        return list(edges)
    for _ in range(5 * len(keys)):
        a, b = rng.choice(len(keys), 2, replace=False)
        A_, B_ = keys[a], keys[b]
        sa, sb = adj[A_], adj[B_]
        com = sa & sb
        ea = list(sa - com); eb = list(sb - com)
        pool = ea + eb
        if not pool:
            continue
        rng.shuffle(pool)
        adj[A_] = com | set(pool[:len(ea)]); adj[B_] = com | set(pool[len(ea):])
    return [(u, v) for u in adj for v in adj[u]]


keys = ['largest_SCC', 'flow_hierarchy', 'global_reaching_centrality',
        'trophic_incoherence_F0', 'condensation_depth_layers', 'transitivity']
acc = collections.defaultdict(list)
for i in range(NRAND):
    ed = []
    for k, es in byc.items():
        ed += curveball(es)
    Gq = nx.DiGraph(); Gq.add_nodes_from(sorted(REG)); Gq.add_edges_from(ed)
    d = describe(Gq, 'null')
    for k in keys:
        acc[k].append(d[k])
    if (i + 1) % 50 == 0:
        print(f'  null {i+1}/{NRAND}', flush=True)

obs = rows[1]
nrows = []
print('\nregulator subnetwork vs degree-preserving null '
      f'({NRAND} curveball randomisations within edge class):')
for k in keys:
    v = np.array(acc[k], float)
    o = float(obs[k])
    sd = v.std(ddof=1)
    z = (o - v.mean()) / sd if sd > 0 else float('nan')
    p_hi = (np.sum(v >= o) + 1) / (NRAND + 1)
    p_lo = (np.sum(v <= o) + 1) / (NRAND + 1)
    p = min(1.0, 2 * min(p_hi, p_lo))
    nrows.append(dict(metric=k, observed=o, null_mean=float(v.mean()), null_sd=float(sd),
                      z=float(z), p_two_sided=float(p), n_randomisations=NRAND))
    print(f'  {k:<32} obs={o:<10.4f} null={v.mean():.4f} +/- {sd:.4f}  z={z:+.2f}  p={p:.4g}')

with open(f'{RES}/systems_regulator_subnetwork.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(f'{RES}/systems_regulator_subnetwork_nulls.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(nrows[0].keys())); w.writeheader(); w.writerows(nrows)
json.dump(dict(observed=rows, nulls=nrows,
               n_genes=len(genes),
               n_genes_out_degree_zero=sum(1 for g in genes if od[g] == 0)),
          open(f'{RES}/systems_regulator_subnetwork.json', 'w'), indent=2)
print('\nwrote systems_regulator_subnetwork.{csv,json} and _nulls.csv')
