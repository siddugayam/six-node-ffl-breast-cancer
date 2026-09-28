#!/usr/bin/env python3
"""
08_motif_node_sets.py
=====================
Build the node sets used for functional enrichment of each FFL motif class.

GRAPH
-----
G_full  = canonical_edges.tsv  +  layer_TF_target.tsv (directed)
                               +  layer_gene_gene.tsv (as deposited, directed field ignored -> directed as given)
                               +  layer_miRNA_miRNA.tsv (10 kb co-transcription, added in BOTH directions)
          Reciprocal TF<->miRNA pairs are RETAINED (needed to type Composite-FFLs).
G_dag   = G_full with reciprocal TF<->miRNA pairs contracted to the transcriptional arc
          (TF->miRNA kept, miRNA->TF dropped).  This is the graph used by 03_ffl_census.py
          and is required because the generalised FFL definition demands acyclicity.

MOTIF CLASSES
-------------
3-node (classical, Zhang et al. 2015 nomenclature; enumerated EXHAUSTIVELY on G_full):
    core = (R, M, T) with R->M, R->T, M->T, T not a miRNA
    3-TF   : R = TF,    M = miRNA   (TF is the master regulator)
    3-miR  : R = miRNA, M = TF      (miRNA is the master regulator)
    3-Comp : the TF<->miRNA pair is RECIPROCAL, i.e. the core is both a TF-FFL and a
             miRNA-FFL on the same three nodes.
    NOTE ON EXCLUSIVITY.  These classes are defined INCLUSIVELY: a core with a reciprocal
    TF<->miRNA pair is counted in 3-TF, in 3-miR AND in 3-Comp.  On the merged canonical
    network 1223 of the TF-miRNA pairs are reciprocal, so an exclusive definition (the one
    used in 02x_ffl_core_crosscheck.py) empties the 3-TF class down to 9 cores / 14 nodes
    and makes a comparative enrichment meaningless.  Exclusive counts are reported
    alongside for transparency.

n-node generalised FFL (n = 4,5,6; Kashtan et al. 2004 "topological generalisation"),
enumerated on G_dag with RAND-ESU (Wernicke 2006).  A vertex-induced connected subgraph S
on n vertices is an n-node FFL iff
    (D1) S acyclic; (D2) exactly one source and one sink;
    (D3) every vertex lies on a directed source->sink path;
    (D4) >= 2 internally vertex-disjoint source->sink paths.
Exhaustive enumeration is impossible at n=5,6 (6.3e8 / 3.2e10 connected induced subgraphs),
so RAND-ESU sampling is used and the NODE UNION is reported with a rarefaction curve so the
reader can see whether the node set has saturated.  Sampling is unbiased over subgraphs.

OUTPUT
------
results/motif_node_sets.tsv        motif_set, node, type
results/motif_node_set_summary.tsv motif_set, n_ffl, n_nodes, n_TF, n_Gene, n_miRNA, n_protein_coding
results/motif_rarefaction.tsv      motif_set, n_ffl_seen, n_nodes_union, n_protein_coding_union
"""
import csv, collections, json, os, random, sys, time

REV = '/path/to/revision'
LOG = open(f'{REV}/logs/08_motif_node_sets.log', 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + '\n')
    LOG.flush()


# ---------------------------------------------------------------- graph build
def build():
    nodes = {r['name']: r['type'] for r in
             csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}
    edges = {}
    prov = collections.Counter()
    for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t'):
        edges[(r['source'], r['target'])] = r['edge_type']
        prov['canonical'] += 1

    def add(path, sc, tc, et, undirected=False, filt=None):
        n = 0
        if not os.path.exists(path):
            log(f'  WARNING missing {path}')
            return 0
        for r in csv.DictReader(open(path), delimiter='\t'):
            if filt and not filt(r):
                continue
            a, b = r.get(sc), r.get(tc)
            if a not in nodes or b not in nodes or a == b:
                continue
            for e in ([(a, b), (b, a)] if undirected else [(a, b)]):
                if e not in edges:
                    edges[e] = et
                    n += 1
        prov[et + '_added'] += n
        return n

    add(f'{REV}/data/layer_TF_target.tsv', 'source', 'target', 'TF_target')
    add(f'{REV}/data/layer_gene_gene.tsv', 'source', 'target', 'gene_gene')
    add(f'{REV}/data/layer_miRNA_miRNA.tsv', 'miRNA_1', 'miRNA_2', 'miRNA_miRNA',
        undirected=True,
        filt=lambda r: str(r.get('threshold_10kb', 'TRUE')).upper() in ('TRUE', '1'))
    log('graph provenance:', dict(prov))
    return nodes, edges


def contract(nodes, edges):
    e2 = dict(edges)
    n = 0
    for (a, b) in list(e2):
        if nodes.get(a) == 'TF' and nodes.get(b) == 'miRNA' and (b, a) in e2:
            if e2.pop((b, a), None) is not None:
                n += 1
    return e2, n


# ------------------------------------------------------- 3-node classical FFL
def three_node(nodes, edges):
    Eset = set(edges)
    out = collections.defaultdict(set)
    for a, b in Eset:
        out[a].add(b)
    tR = lambda x: nodes.get(x, 'Gene')
    res = collections.defaultdict(list)
    excl = collections.Counter()
    for R in list(out):
        for M in out[R]:
            for T in (out[R] & out[M]):
                if T == R or T == M or tR(T) == 'miRNA':
                    continue
                recip = (M, R) in Eset
                if tR(R) == 'TF' and tR(M) == 'miRNA':
                    res['3-TF'].append((R, M, T))
                    if recip:
                        res['3-Comp'].append((R, M, T))
                        excl['Composite-FFL(exclusive)'] += 1
                    else:
                        excl['TF-FFL(exclusive)'] += 1
                elif tR(R) == 'miRNA' and tR(M) == 'TF':
                    res['3-miR'].append((R, M, T))
                    if not recip:
                        excl['miRNA-FFL(exclusive)'] += 1
                else:
                    continue
    log('  3-node EXCLUSIVE core counts (crosscheck convention):', dict(excl))
    return res


# ---------------------------------------------------------------- ESU + FFL test
def is_ffl(S, out_s, in_s):
    src = [v for v in S if not in_s[v]]
    if len(src) != 1:
        return None
    snk = [v for v in S if not out_s[v]]
    if len(snk) != 1:
        return None
    s, t = src[0], snk[0]
    colour = {}

    def visit(u):
        stack = [(u, iter(out_s[u]))]
        colour[u] = 1
        while stack:
            node, it = stack[-1]
            adv = False
            for v in it:
                c = colour.get(v)
                if c == 1:
                    return False
                if c is None:
                    colour[v] = 1
                    stack.append((v, iter(out_s[v])))
                    adv = True
                    break
            if not adv:
                colour[node] = 2
                stack.pop()
        return True

    for v in S:
        if colour.get(v) is None and not visit(v):
            return None
    reach = {s}
    st = [s]
    while st:
        u = st.pop()
        for v in out_s[u]:
            if v not in reach:
                reach.add(v)
                st.append(v)
    back = {t}
    st = [t]
    while st:
        u = st.pop()
        for v in in_s[u]:
            if v not in back:
                back.add(v)
                st.append(v)
    if not S <= (reach & back):
        return None

    def find_path(banned):
        prev = {s: None}
        st = [s]
        while st:
            u = st.pop()
            if u == t:
                break
            for v in out_s[u]:
                if v in prev or (v in banned and v != t):
                    continue
                prev[v] = u
                st.append(v)
        if t not in prev:
            return None
        p, u = [], t
        while u is not None:
            p.append(u)
            u = prev[u]
        return p[::-1]

    p1 = find_path(set())
    if p1 is None:
        return None
    if find_path(set(p1[1:-1])) is None:
        return None
    return (s, t)


def esu_ffl(nodes, edges, k, probs, seed, callback, max_seconds=None):
    out = collections.defaultdict(set)
    inn = collections.defaultdict(set)
    adj = collections.defaultdict(set)
    for a, b in edges:
        out[a].add(b)
        inn[b].add(a)
        adj[a].add(b)
        adj[b].add(a)
    nl = sorted([v for v in nodes if adj[v]])
    idx = {v: i for i, v in enumerate(nl)}
    rng = random.Random(seed)
    stats = collections.Counter()
    t0 = time.time()
    stop = [False]

    def extend(sub, ext, vidx, excl):
        if stop[0]:
            return
        if len(sub) == k:
            S = set(sub)
            os_ = {v: out[v] & S for v in S}
            is_ = {v: inn[v] & S for v in S}
            stats['subgraphs'] += 1
            if max_seconds and stats['subgraphs'] % 200000 == 0 and time.time() - t0 > max_seconds:
                stop[0] = True
            r = is_ffl(S, os_, is_)
            if r:
                stats['ffl'] += 1
                callback(sub, r)
            return
        d = len(sub)
        p = probs[d]
        ext = set(ext)
        while ext:
            w = ext.pop()
            if p < 1.0 and rng.random() >= p:
                continue
            newext = ext | {u for u in adj[w] if idx[u] > vidx and u not in excl}
            extend(sub + [w], newext, vidx, excl | {w} | adj[w])
            if stop[0]:
                return

    for v in nl:
        if probs[0] < 1.0 and rng.random() >= probs[0]:
            continue
        i = idx[v]
        extend([v], {u for u in adj[v] if idx[u] > i}, i, {v} | adj[v])
        if stop[0]:
            break
    scale = 1.0
    for p in probs:
        scale /= p
    stats['elapsed'] = round(time.time() - t0, 1)
    stats['stopped_early'] = int(stop[0])
    return stats, scale


# ---------------------------------------------------------------- main
def main():
    nodes, edges = build()
    log(f'G_full : {len(nodes)} nodes, {len(edges)} directed edges')
    edges_dag, nrec = contract(nodes, edges)
    log(f'G_dag  : {len(edges_dag)} directed edges after contracting {nrec} reciprocal TF<->miRNA pairs')

    node_sets = {}
    ffl_counts = {}
    rare_rows = []

    # ---- 3-node, exhaustive
    tn = three_node(nodes, edges)
    for cls in ('3-miR', '3-TF', '3-Comp'):
        cores = tn.get(cls, [])
        ffl_counts[cls] = len(cores)
        ns = set()
        seen = 0
        for (R, M, T) in cores:
            before = len(ns)
            ns.update((R, M, T))
            seen += 1
            if seen % 25 == 0 or before != len(ns):
                pass
        node_sets[cls] = ns
        log(f'{cls}: {len(cores)} FFL cores, {len(ns)} distinct nodes')
        # rarefaction
        acc = set()
        for i, (R, M, T) in enumerate(cores, 1):
            acc.update((R, M, T))
            if i % max(1, len(cores) // 40) == 0 or i == len(cores):
                rare_rows.append((cls, i, len(acc),
                                  sum(1 for x in acc if nodes[x] != 'miRNA')))

    # ---- 4/5/6-node generalised FFL, RAND-ESU
    plans = {
        4: dict(probs=[1.0, 1.0, 1.0, 0.30], seed=11, max_seconds=900),
        5: dict(probs=[1.0, 1.0, 0.60, 0.20, 0.06], seed=12, max_seconds=1500),
        6: dict(probs=[1.0, 1.0, 0.50, 0.12, 0.05, 0.05], seed=13, max_seconds=2400),
    }
    for k in (4, 5, 6):
        pl = plans[k]
        cls = f'{k}-node'
        acc = set()
        comp = collections.Counter()
        nffl = [0]

        def cb(sub, r, _acc=acc, _comp=comp, _n=nffl):
            _n[0] += 1
            _acc.update(sub)
            tc = collections.Counter(nodes[v] for v in sub)
            _comp[(tc['TF'], tc['miRNA'], tc['Gene'])] += 1
            if _n[0] % 250 == 0:
                rare_rows.append((cls, _n[0], len(_acc),
                                  sum(1 for x in _acc if nodes[x] != 'miRNA')))

        st, scale = esu_ffl(nodes, edges_dag, k, pl['probs'], pl['seed'], cb,
                            max_seconds=pl['max_seconds'])
        node_sets[cls] = acc
        ffl_counts[cls] = st['ffl']
        rare_rows.append((cls, st['ffl'], len(acc),
                          sum(1 for x in acc if nodes[x] != 'miRNA')))
        log(f'{cls}: probs={pl["probs"]} sampled_subgraphs={st["subgraphs"]} '
            f'sampled_FFL={st["ffl"]} scale={scale:.1f} '
            f'est_total_FFL={st["ffl"] * scale:,.0f} nodes_union={len(acc)} '
            f'elapsed={st["elapsed"]}s stopped_early={st["stopped_early"]}')
        json.dump({'k': k, 'probs': pl['probs'], 'seed': pl['seed'],
                   'sampled_subgraphs': st['subgraphs'], 'sampled_ffl': st['ffl'],
                   'scale': scale, 'est_ffl': st['ffl'] * scale,
                   'elapsed_s': st['elapsed'], 'stopped_early': bool(st['stopped_early']),
                   'nodes_union': len(acc),
                   'composition': {f'TF{a}_miR{b}_G{c}': v for (a, b, c), v in comp.items()}},
                  open(f'{REV}/results/motif_enum_k{k}.json', 'w'), indent=1)

    # ---- author exemplars
    for tag, f in (('exemplar_4node', '4-TF.csv'), ('exemplar_5node', '5-TF.csv'),
                   ('exemplar_6node', '6-TF.csv')):
        rows = list(csv.DictReader(open(
            f'/path/to/revision/analyses/original_submission_code/node_attributes/{f}')))
        node_sets[tag] = {r['name'] for r in rows}
        ffl_counts[tag] = 1
        log(f'{tag}: {len(rows)} nodes ('
            f'{sum(1 for r in rows if r["Type"] != "miRNA")} protein-coding)')

    # ---- write
    with open(f'{REV}/results/motif_node_sets.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['motif_set', 'node', 'type'])
        for cls in ('3-miR', '3-TF', '3-Comp', '4-node', '5-node', '6-node',
                    'exemplar_4node', 'exemplar_5node', 'exemplar_6node'):
            for n in sorted(node_sets[cls]):
                w.writerow([cls, n, nodes.get(n, 'Gene')])

    with open(f'{REV}/results/motif_node_set_summary.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['motif_set', 'n_ffl', 'n_nodes', 'n_TF', 'n_Gene', 'n_miRNA',
                    'n_protein_coding'])
        for cls in ('3-miR', '3-TF', '3-Comp', '4-node', '5-node', '6-node',
                    'exemplar_4node', 'exemplar_5node', 'exemplar_6node'):
            ns = node_sets[cls]
            c = collections.Counter(nodes.get(x, 'Gene') for x in ns)
            w.writerow([cls, ffl_counts[cls], len(ns), c['TF'], c['Gene'], c['miRNA'],
                        c['TF'] + c['Gene']])
            log(f'SUMMARY {cls}: n_ffl={ffl_counts[cls]} nodes={len(ns)} '
                f'TF={c["TF"]} Gene={c["Gene"]} miRNA={c["miRNA"]} '
                f'protein_coding={c["TF"] + c["Gene"]}')

    with open(f'{REV}/results/motif_rarefaction.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['motif_set', 'n_ffl_seen', 'n_nodes_union', 'n_protein_coding_union'])
        w.writerows(rare_rows)

    # ---- edge list of the graphs actually used (for the topology step in R)
    with open(f'{REV}/results/graph_used_edges.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['source', 'target', 'edge_type'])
        for (a, b), et in sorted(edges.items()):
            w.writerow([a, b, et])
    log('wrote results/motif_node_sets.tsv, motif_node_set_summary.tsv, '
        'motif_rarefaction.tsv, graph_used_edges.tsv')


if __name__ == '__main__':
    main()
