#!/usr/bin/env python3
"""
Topology and MCC hubs of the merged 3-to-6-node composite FFL network, as in Bhat et al. (2024)
section 2.7 (Cytoscape NetworkAnalyzer + cytoHubba MCC).  Input: sif/<graph>_composite_3to6_merged.sif
from build_bhat_networks.py.  The network is analysed as undirected, as NetworkAnalyzer does with
"treat the network as undirected".

Definitions (NetworkAnalyzer, undirected):
  degree_edges     edges touching the node, counted as in the SIF: a reciprocal M<->T or T<->T pair is 2 edges
  neighbours       distinct adjacent nodes; every metric below uses the simple graph (reciprocal pairs collapsed)
  avg_shortest_path_length   mean distance to every node reachable from n
  closeness        1 / avg_shortest_path_length
  betweenness      sum over pairs s,t (s,t != n) of sigma_st(n)/sigma_st, divided by (N-1)(N-2)/2 with N the
                   size of n's connected component
  clustering       2 e_n / (k_n (k_n - 1)); 0 when k_n < 2
  neighbourhood_connectivity  mean number of neighbours of n's neighbours
  topological_coefficient     mean over m sharing >= 1 neighbour with n of J(n,m), divided by k_n; J(n,m) = shared
                   neighbours + 1 if n and m are adjacent; NA when k_n < 2
  MCC (cytoHubba)  sum over maximal cliques C containing n of (|C| - 1)!; equals the degree when no two neighbours
                   of n are adjacent
"""
import os, sys, math, collections
import networkx as nx
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))


def load(tag, stem):
    ty = dict(l.rstrip('\n').split('\t') for l in open(os.path.join(HERE, 'sif', f'{tag}_node_attributes.tsv')).readlines()[1:])
    E = [l.rstrip('\n').split('\t') for l in open(os.path.join(HERE, 'sif', f'{tag}_{stem}.sif')) if l.strip()]
    return ty, E


def topological_coefficient(G, n):
    Nn = set(G[n]); k = len(Nn)
    if k < 2: return None
    J = collections.Counter()
    for u in Nn:
        for m in G[u]:
            if m != n: J[m] += 1
    if not J: return None
    vals = [J[m] + (1 if m in Nn else 0) for m in J]
    return sum(vals) / len(vals) / k


def analyse(tag, stem='composite_3to6_merged'):
    ty, E = load(tag, stem)
    G = nx.Graph()
    deg_edges = collections.Counter()
    for s, lab, t in E:
        deg_edges[s] += 1; deg_edges[t] += 1
        G.add_edge(s, t)
    bet = {}
    for comp in nx.connected_components(G):
        bet.update(nx.betweenness_centrality(G.subgraph(comp), normalized=True))
    clo = nx.closeness_centrality(G, wf_improved=False)
    clu = nx.clustering(G)
    aspl = {}
    for n, d in nx.all_pairs_shortest_path_length(G):
        r = [v for k, v in d.items() if k != n]
        aspl[n] = sum(r) / len(r) if r else 0.0
    mcc = collections.Counter()
    ncl = collections.Counter()
    for C in nx.find_cliques(G):
        f = math.factorial(len(C) - 1)
        for v in C: mcc[v] += f; ncl[v] += 1
    rows = []
    for n in G.nodes:
        nb = list(G[n])
        tc = topological_coefficient(G, n)
        rows.append(dict(node=n, type=ty[n], degree_edges=deg_edges[n], neighbours=len(nb),
                         avg_shortest_path_length=aspl[n], closeness=clo[n], betweenness=bet[n],
                         clustering=clu[n], neighbourhood_connectivity=sum(G.degree(u) for u in nb) / len(nb),
                         topological_coefficient='NA' if tc is None else tc, MCC=mcc[n], maximal_cliques=ncl[n]))
    rows.sort(key=lambda r: (-r['MCC'], -r['neighbours'], r['node']))
    for i, r in enumerate(rows, 1): r['MCC_rank'] = i
    hdr = list(rows[0].keys())
    fmt = lambda v: f'{v:.6g}' if isinstance(v, float) else str(v)
    with open(os.path.join(HERE, f'topology_{tag}_{stem}.csv'), 'w') as f:
        f.write(','.join(hdr) + '\n')
        for r in rows: f.write(','.join(fmt(r[h]) for h in hdr) + '\n')
    comps = sorted((len(c) for c in nx.connected_components(G)), reverse=True)
    big = G.subgraph(max(nx.connected_components(G), key=len))
    summ = dict(graph=tag, network=stem, sif_edges=len(E), nodes=G.number_of_nodes(), undirected_edges=G.number_of_edges(),
                nodes_by_type=dict(collections.Counter(ty[n] for n in G.nodes)),
                sif_edges_by_label=dict(collections.Counter(lab for _, lab, _ in E)),
                connected_components=len(comps), largest_component=comps[0],
                density=nx.density(G), avg_neighbours=2 * G.number_of_edges() / G.number_of_nodes(),
                clustering_coefficient_mean=sum(clu.values()) / len(clu),
                characteristic_path_length_largest_component=nx.average_shortest_path_length(big),
                diameter_largest_component=nx.diameter(big),
                maximal_cliques=sum(1 for _ in nx.find_cliques(G)),
                max_clique_size=max(len(c) for c in nx.find_cliques(G)))
    return rows, summ


if __name__ == '__main__':
    import json
    out = []
    for tag in ('nolegacy', 'nolegacy_nostring'):
        rows, summ = analyse(tag)
        out.append(summ)
        print(f"\n== {tag}: {summ['nodes']} nodes {summ['nodes_by_type']}, {summ['sif_edges']} SIF edges "
              f"({summ['undirected_edges']} undirected), components {summ['connected_components']}, "
              f"CPL {summ['characteristic_path_length_largest_component']:.3f}, diameter {summ['diameter_largest_component']}, "
              f"mean clustering {summ['clustering_coefficient_mean']:.3f}, max clique {summ['max_clique_size']}")
        print('   top 10 MCC:')
        for r in rows[:10]:
            print(f"   {r['MCC_rank']:>3} {r['node']:18s} {r['type']:6s} MCC={r['MCC']:<8} neighbours={r['neighbours']:<4} "
                  f"betweenness={r['betweenness']:.4f} closeness={r['closeness']:.4f}")
        for t in ('TF', 'miRNA', 'Gene'):
            top = sorted([r for r in rows if r['type'] == t], key=lambda r: (-r['neighbours'], r['node']))[:5]
            print(f'   top 5 {t} by neighbours: ' + ', '.join(f"{r['node']} ({r['neighbours']})" for r in top))
    json.dump(out, open(os.path.join(HERE, 'topology_network_summary.json'), 'w'), indent=1)
