#!/usr/bin/env python3
"""F3: node compendium section a (Note S5 source) recomputed on the analysed network.
Method copied from scripts/v5/40_node_compendium_assemble.py section 1: networkx DiGraph on
data/canonical_edges.tsv (all 587 nodes added first), in/out-degree, directed betweenness (raw and
normalised) and betweenness of the undirected simple graph.  Here the 30 miRNA_miRNA rows are dropped
(analysed network, 6,829 edges).  The stored values are recomputed with the same code on the full
deposit and checked against results/v5/node_compendium_table.csv first.  Read-only."""
import os, csv
import pandas as pd, networkx as nx
REV = '/path/to/revision'
H = os.path.dirname(os.path.abspath(__file__))
E = pd.read_csv(f'{REV}/data/canonical_edges.tsv', sep='\t')
NODES = pd.read_csv(f'{REV}/data/canonical_nodes.tsv', sep='\t')['name']
thirty = [r['name'] for r in csv.DictReader(open(f'{REV}/results/v5/top10_per_class.csv'))]
def topo(E):
    G = nx.DiGraph(); G.add_nodes_from(NODES)
    for s, t in zip(E.source, E.target): G.add_edge(s, t)
    U = nx.Graph(); U.add_nodes_from(G.nodes()); U.add_edges_from(G.edges())
    return dict(in_degree=dict(G.in_degree()), out_degree=dict(G.out_degree()),
                betweenness_directed=nx.betweenness_centrality(G, normalized=False),
                betweenness_directed_norm=nx.betweenness_centrality(G, normalized=True),
                betweenness_undirected=nx.betweenness_centrality(U, normalized=False)), G.number_of_edges()
full, n0 = topo(E); ana, n1 = topo(E[E.edge_type != 'miRNA_miRNA'])
C = pd.read_csv(f'{REV}/results/v5/node_compendium_table.csv').set_index('name')
cols = ['in_degree', 'out_degree', 'betweenness_directed', 'betweenness_directed_norm', 'betweenness_undirected']
ok = all(abs(full[c][n] - float(C.loc[n, c])) < 1e-6 * max(1.0, abs(full[c][n])) for n in thirty for c in cols if c in C.columns)
print(f'graphs: deposit {n0} edges, analysed network {n1} edges; stored compendium values reproduced on the deposit: {ok}')
rows = []
for n in thirty:
    r = dict(node=n, type=C.loc[n, 'type'] if 'type' in C.columns else '')
    changed = []
    for c in cols:
        a, b = full[c][n], ana[c][n]; r[f'{c}_deposit'] = a; r[f'{c}_analysed'] = b
        if abs(a - b) > 1e-9 * max(1.0, abs(a)): changed.append(c)
    r['changed'] = ';'.join(changed); rows.append(r)
pd.DataFrame(rows).to_csv(f'{H}/f3_compendium_topology_thirty.csv', index=False)
ch = [r for r in rows if r['changed']]
print(f'{len(ch)} of 30 nodes change in at least one value:')
for r in ch:
    print('  ', r['node'], '|', ', '.join(f"{c}: {r[c + '_deposit']:.6g} -> {r[c + '_analysed']:.6g}" for c in r['changed'].split(';')))
