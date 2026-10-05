#!/usr/bin/env python3
"""Section T of SETTINGS.md: topology of the twelve Bhat FFL networks.

T1  Composition: instances, distinct node sets, nodes by type, and edges by type and by source layer (from the
    deposited files).
T2  Degree and hubs. Each network is analysed as undirected, with the NetworkAnalyzer definitions of the Bhat topology
    script:
    - degree counted as edges (a reciprocal pair counts 2), and neighbours;
    - betweenness within each connected component, normalised;
    - closeness and clustering;
    - cytoHubba MCC: the sum over the maximal cliques C containing a node of (|C| - 1)!.
    The same measures are computed on the analysed network (587 nodes, 6,829 edges).
T3  Jaccard index of the node sets and of the edge sets for every pair of networks. An undirected link is one unordered
    pair; a directed arc is an ordered pair.
T4  Overlap with the typed cores (Table S4 in the paper): the share of each network's nodes that belong to any typed
    core, and of its edges that some typed core uses. A core uses R->M, M->T and R->T, plus M->R for a composite.
T5  The share of distinct instance node sets that are D1-D4 modules of the same size on the census graph. The graph is
    the census script's build_graph(True), with its TF<->miRNA contraction; the test is its is_ffl on the induced
    subgraph.

usage: python3 topology.py <analysis root> <network deposit folder> <original SIF folder> <census script> <output folder>
"""
import os, sys, csv, json, math, itertools, importlib.util, collections
import numpy as np
import networkx as nx
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bhat_common as bc

ANALYSED = 'analysed_network'
METRICS = {'MCC': lambda r: (-r['MCC'], -r['neighbours'], r['node']),
           'neighbours': lambda r: (-r['neighbours'], r['node']),
           'betweenness': lambda r: (-r['betweenness'], r['node'])}


def metrics(E, ty):
    """E: (source, target) lines of one network, each line one edge. Node rows and a graph summary."""
    G = nx.Graph(); deg = collections.Counter()
    for s, t in E:
        deg[s] += 1; deg[t] += 1; G.add_edge(s, t)
    bet = {}
    for comp in nx.connected_components(G):
        bet.update(nx.betweenness_centrality(G.subgraph(comp), normalized=True))
    clo = nx.closeness_centrality(G, wf_improved=False)
    clu = nx.clustering(G)
    mcc = collections.Counter(); ncl = 0; maxc = 0
    for C in nx.find_cliques(G):
        f = math.factorial(len(C) - 1); ncl += 1; maxc = max(maxc, len(C))
        for v in C: mcc[v] += f
    rows = [dict(node=v, type=ty[v], degree_edges=deg[v], neighbours=G.degree(v), betweenness=bet[v],
                 closeness=clo[v], clustering=clu[v], MCC=mcc[v]) for v in G.nodes]
    comps = sorted((len(c) for c in nx.connected_components(G)), reverse=True)
    summ = dict(nodes=G.number_of_nodes(), edges_as_lines=len(E), undirected_edges=G.number_of_edges(),
                components=len(comps), largest_component=comps[0], density=nx.density(G),
                mean_clustering=sum(clu.values()) / len(clu), maximal_cliques=ncl, max_clique_size=maxc)
    for m in ('degree_edges', 'neighbours'):
        x = np.array([r[m] for r in rows], dtype=float)
        q = np.percentile(x, [0, 25, 50, 75, 100])            # numpy default (linear interpolation)
        summ.update({f'{m}_min': q[0], f'{m}_q1': q[1], f'{m}_median': q[2], f'{m}_q3': q[3], f'{m}_max': q[4],
                     f'{m}_mean': x.mean(), f'{m}_share_eq_1': float((x == 1).mean())})
    return rows, summ


def top10(rows):
    return {m: sorted(rows, key=k)[:10] for m, k in METRICS.items()}


def write_csv(path, rows, hdr=None):
    hdr = hdr or list(rows[0].keys())
    fmt = lambda v: f'{v:.6g}' if isinstance(v, float) else str(v)
    with open(path, 'w', newline='') as f:
        f.write(','.join(hdr) + '\n')
        for r in rows: f.write(','.join(fmt(r.get(h, '')) for h in hdr) + '\n')


def main(root, netdir, sifdir, census_script, out):
    os.makedirs(out, exist_ok=True)
    N = bc.networks(root, netdir, sifdir)
    nodes, analysed, inst, edges, vtab = N['nodes'], N['analysed'], N['inst'], N['edges'], N['vtab']
    log = []; say = lambda *a: (print(*a), log.append(' '.join(map(str, a))))
    say('networks rebuilt byte-identical to the deposit (md5 of all', len(N['manifest']), 'listed files)')

    # T1 composition
    types = ['TF_miRNA', 'miRNA_target', 'TF_target', 'gene_gene', 'miRNA_miRNA']
    layers = ['analysed network', 'census layer: TRRUST TF-target', 'census layer: STRING gene-gene',
              'census layer: miRNA pair within 10 kb']
    T1 = []
    for net in bc.NETWORKS:
        vt = collections.Counter(r['type'] for r in vtab[net]); et = collections.Counter(r['edge_type'] for r in edges[net])
        ly = collections.Counter(r['layer'] for r in edges[net])
        assert set(et) <= set(types) and set(ly) <= set(layers)
        T1.append(dict(network=net, instances=len(inst[net]), node_sets=len({frozenset(i) for i in inst[net]}),
                       nodes=len(vtab[net]), TF=vt['TF'], miRNA=vt['miRNA'], Gene=vt['Gene'], edges=len(edges[net]),
                       **{f'type_{t}': et[t] for t in types},
                       **{'layer_' + l.replace('census layer: ', '').replace(' ', '_'): ly[l] for l in layers}))
    write_csv(os.path.join(out, 'T1_composition.csv'), T1)
    say('\nT1 composition')
    for r in T1:
        say(f"  {r['network']:20s} instances {r['instances']:6d} sets {r['node_sets']:6d} nodes {r['nodes']:4d} "
            f"(TF {r['TF']}, miRNA {r['miRNA']}, Gene {r['Gene']}) edges {r['edges']:5d}")

    # T2 degree and hubs
    graphs = {ANALYSED: (list(analysed.keys()), nodes)}
    for net in bc.NETWORKS: graphs[net] = ([(r['source'], r['target']) for r in edges[net]], nodes)
    node_rows, summ_rows, top = [], [], {}
    for g, (E, ty) in graphs.items():
        rows, s = metrics(E, ty)
        for r in rows: node_rows.append(dict(network=g, **r))
        summ_rows.append(dict(network=g, **s)); top[g] = top10(rows)
    write_csv(os.path.join(out, 'T2_node_metrics.csv'),
              sorted(node_rows, key=lambda r: (list(graphs).index(r['network']), -r['MCC'], -r['neighbours'], r['node'])))
    write_csv(os.path.join(out, 'T2_degree_summary.csv'), summ_rows)
    T2top = []
    for g in graphs:
        for m in METRICS:
            ref = {r['node'] for r in top[ANALYSED][m]}
            for i, r in enumerate(top[g][m], 1):
                T2top.append(dict(network=g, metric=m, rank=i, node=r['node'], type=r['type'], value=r[m],
                                  neighbours=r['neighbours'], in_analysed_network_top10=int(r['node'] in ref)))
    write_csv(os.path.join(out, 'T2_top10.csv'), T2top)
    ov = []
    for g in bc.NETWORKS:
        ov.append(dict(network=g, **{f'{m}_top10_shared_with_analysed_network':
                                     len({r['node'] for r in top[g][m]} & {r['node'] for r in top[ANALYSED][m]}) for m in METRICS}))
    write_csv(os.path.join(out, 'T2_top10_overlap_with_analysed_network.csv'), ov)
    say('\nT2 degree (edges) and neighbours; top 10 by MCC')
    for s in summ_rows:
        g = s['network']
        say(f"  {g:20s} nodes {s['nodes']:4d} undirected edges {s['undirected_edges']:5d} components {s['components']:3d} "
            f"degree min/q1/median/q3/max {s['degree_edges_min']:.0f}/{s['degree_edges_q1']:g}/{s['degree_edges_median']:g}/"
            f"{s['degree_edges_q3']:g}/{s['degree_edges_max']:.0f} mean {s['degree_edges_mean']:.2f} share deg 1 "
            f"{s['degree_edges_share_eq_1']:.3f}; neighbours median {s['neighbours_median']:g} share 1 {s['neighbours_share_eq_1']:.3f}")
        say('      MCC top 10: ' + ', '.join(r['node'] for r in top[g]['MCC']))
        say('      neighbours top 10: ' + ', '.join(r['node'] for r in top[g]['neighbours']))
        say('      betweenness top 10: ' + ', '.join(r['node'] for r in top[g]['betweenness']))
    for r in ov: say(f"  {r['network']:20s} top-10 shared with the analysed network: " + ', '.join(f'{m} {r[m + "_top10_shared_with_analysed_network"]}' for m in METRICS))

    # T3 Jaccard between networks
    V = {net: {r['node'] for r in vtab[net]} for net in bc.NETWORKS}
    E = {net: {bc.edge_key(r['source'], r['target'], r['edge_type']) for r in edges[net]} for net in bc.NETWORKS}
    T3 = []
    for a, b in itertools.combinations(bc.NETWORKS, 2):
        T3.append(dict(network_1=a, network_2=b, nodes_shared=len(V[a] & V[b]), nodes_union=len(V[a] | V[b]),
                       jaccard_nodes=len(V[a] & V[b]) / len(V[a] | V[b]), edges_shared=len(E[a] & E[b]),
                       edges_union=len(E[a] | E[b]), jaccard_edges=len(E[a] & E[b]) / len(E[a] | E[b])))
    write_csv(os.path.join(out, 'T3_jaccard.csv'), T3)
    say('\nT3 Jaccard (nodes / edges)')
    for r in T3: say(f"  {r['network_1']:20s} {r['network_2']:20s} {r['jaccard_nodes']:.3f} / {r['jaccard_edges']:.3f}")

    # T4 overlap with the typed cores
    s2 = list(csv.DictReader(open(os.path.join(root, 'results/v5/tables/TableS2_ffl_cores.csv'), newline='')))
    assert len(s2) == 1649
    cn, ce = set(), set()
    for r in s2:
        R, M, T = r['regulator'], r['intermediate'], r['target']
        cn |= {R, M, T}; ce |= {(R, M), (M, T), (R, T)}
        if r['class'] == 'Composite-FFL': ce.add((M, R))
    assert ce <= set(analysed), 'a typed-core edge is not in the analysed network'
    T4 = []
    for net in bc.NETWORKS:
        en = [bc.edge_key(r['source'], r['target'], r['edge_type']) for r in edges[net]]
        T4.append(dict(network=net, nodes=len(V[net]), nodes_in_a_typed_core=len(V[net] & cn),
                       share_nodes=len(V[net] & cn) / len(V[net]), edges=len(en),
                       edges_used_by_a_typed_core=sum(e in ce for e in en), share_edges=sum(e in ce for e in en) / len(en),
                       nodes_not_in_any_typed_core=len(V[net] - cn)))
    write_csv(os.path.join(out, 'T4_typed_core_overlap.csv'), T4)
    say(f'\nT4 typed cores: {len(s2)} cores, {len(cn)} nodes, {len(ce)} edges')
    for r in T4:
        say(f"  {r['network']:20s} nodes in a core {r['nodes_in_a_typed_core']}/{r['nodes']} = {r['share_nodes']:.3f}; "
            f"edges used by a core {r['edges_used_by_a_typed_core']}/{r['edges']} = {r['share_edges']:.3f}")

    # T5 D1-D4 modules on the contracted census graph
    spec = importlib.util.spec_from_file_location('census', census_script)
    cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
    cnodes, cedges, nrecip = cs.build_graph(True)
    removed = set(N['census']) - set(cedges)
    assert set(cedges) <= set(N['census']) and len(removed) == nrecip and len(N['census']) == 9226
    assert all(cnodes[a] == 'miRNA' and cnodes[b] == 'TF' and (b, a) in cedges for a, b in removed)
    say(f'\nT5 census graph: {len(N["census"])} arcs, {nrecip} miRNA->TF arcs removed by the contraction, {len(cedges)} left')
    outs, ins = collections.defaultdict(set), collections.defaultdict(set)
    for a, b in cedges: outs[a].add(b); ins[b].add(a)
    T5 = []; why = collections.OrderedDict()
    for net in bc.NETWORKS:
        sets = {frozenset(i) for i in inst[net]}; npass = 0; why[net] = collections.Counter()
        for S in sets:
            S = set(S)
            if cs.is_ffl(S, {v: outs[v] & S for v in S}, {v: ins[v] & S for v in S}): npass += 1; continue
            # why a set fails: the node types and arc types of every 2-cycle in the induced subgraph, if any
            cyc = sorted({' / '.join(sorted(f'{cnodes[a]}->{cnodes[b]} ({cedges[(a, b)]})' for a, b in ((u, v), (v, u))))
                          for u in S for v in outs[u] & S if u in outs[v]})
            why[net]['2-cycle: ' + '; '.join(cyc) if cyc else 'no 2-cycle'] += 1
        T5.append(dict(network=net, node_sets=len(sets), D1_D4_modules=npass, share=npass / len(sets)))
    write_csv(os.path.join(out, 'T5_d1d4.csv'), T5)
    write_csv(os.path.join(out, 'T5_d1d4_failures.csv'),
              [dict(network=net, reason=k, node_sets=v) for net in why for k, v in why[net].most_common()])
    for r in T5:
        say(f"  {r['network']:20s} node sets passing D1-D4: {r['D1_D4_modules']} of {r['node_sets']} = {r['share']:.3f}")
        for k, v in why[r['network']].most_common(): say(f'      fails: {v:5d}  {k}')
    json.dump(dict(T1=T1, T2_summary=summ_rows, T2_overlap=ov, T4=T4, T5=T5), open(os.path.join(out, 'topology_summary.json'), 'w'), indent=1)
    open(os.path.join(out, 'topology.log'), 'w').write('\n'.join(log) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:6])
