#!/usr/bin/env python3
"""
60_census_sensitivity.py
========================
Sensitivity of the n-node feed-forward-loop census to three construction choices
that can reasonably be questioned:

  (1) whether STRING functional associations are admitted as directed regulatory
      arcs.  Methods 2.2 states they are "used only for module connectivity and
      never to define a regulatory arm of an FFL"; the published census script
      (scripts/03_ffl_census/03_ffl_census.py) nevertheless adds every row of
      data/layer_gene_gene.tsv as a DIRECTED arc, and all 834 arcs that survive
      de-duplication are STRING associations.  This script quantifies the cost.
  (2) the miRNA polycistronic clustering window (3 / 10 / 50 kb, or none).
  (3) the evidence tier admitted for miRNA->target edges (all / experimentally
      supported / low-throughput-validated only).  This is the "validated-only
"""
import csv, collections, json, os, random, sys, time

REV = '/path/to/revision'
OUT = os.path.join(REV, 'results', 'v5')
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- graph build
def load_nodes():
    return {r['name']: r['type'] for r in
            csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}


def load_tiers():
    t = {}
    for r in csv.DictReader(open(f'{REV}/data/edge_evidence_tier.tsv'), delimiter='\t'):
        t[(r['source'], r['target'])] = r['tier']
    return t


def build_graph(string_mode='as_directed', mirna_thresh='threshold_10kb',
                mirna_tier='all', drop_author_mirmir=False):
    """string_mode: 'as_directed' (as published) | 'exclude' | 'trrust_only'
       mirna_thresh: 'threshold_3kb' | 'threshold_10kb' | 'threshold_50kb' | None
       mirna_tier  : 'all' | 'experimental' (strong+weak) | 'strong'"""
    nodes = load_nodes()
    tiers = load_tiers()
    keep = {'all': None,
            'experimental': {'strong', 'weak'},
            'strong': {'strong'}}[mirna_tier]

    edges, prov = {}, {}
    n_drop_tier = 0
    for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t'):
        e = (r['source'], r['target'])
        if r['edge_type'] == 'miRNA_target' and keep is not None:
            if tiers.get(e, 'predicted_only') not in keep:
                n_drop_tier += 1
                continue
        if r['edge_type'] == 'miRNA_miRNA' and drop_author_mirmir:
            continue
        edges[e] = r['edge_type']
        prov[e] = 'canonical'

    def add(path, sc, tc, et, undirected=False, filt=None, tag=''):
        n = 0
        if not os.path.exists(path):
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
                    prov[e] = tag
                    n += 1
        return n

    n_tf = add(f'{REV}/data/layer_TF_target.tsv', 'source', 'target', 'TF_target',
               tag='TRRUST_TF_target')

    n_gg = 0
    if string_mode == 'as_directed':
        n_gg = add(f'{REV}/data/layer_gene_gene.tsv', 'source', 'target', 'gene_gene',
                   tag='gene_gene_layer')
    elif string_mode == 'trrust_only':
        n_gg = add(f'{REV}/data/layer_gene_gene.tsv', 'source', 'target', 'gene_gene',
                   filt=lambda r: r['evidence'] == 'TRRUST', tag='gene_gene_TRRUST')
    # 'exclude' -> add nothing

    n_mm = 0
    if mirna_thresh is not None:
        n_mm = add(f'{REV}/data/layer_miRNA_miRNA.tsv', 'miRNA_1', 'miRNA_2',
                   'miRNA_miRNA', undirected=True,
                   filt=lambda r: str(r.get(mirna_thresh, 'TRUE')).upper() in ('TRUE', '1'),
                   tag='miRBase_cluster')

    nrecip = 0
    for (a, b) in list(edges):
        if nodes.get(a) == 'TF' and nodes.get(b) == 'miRNA' and (b, a) in edges:
            if edges.pop((b, a), None) is not None:
                nrecip += 1
                prov.pop((b, a), None)

    meta = dict(n_edges=len(edges), n_recip_contracted=nrecip,
                n_added_TF_target=n_tf, n_added_gene_gene=n_gg,
                n_added_miRNA_miRNA=n_mm, n_miRNA_target_dropped_by_tier=n_drop_tier,
                by_type=dict(collections.Counter(edges.values())))
    return nodes, edges, meta


# ------------------------------------------------------------------ D1-D4 test
def is_ffl(S, out_s, in_s):
    src = [v for v in S if not in_s[v]]
    snk = [v for v in S if not out_s[v]]
    if len(src) != 1 or len(snk) != 1:
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
                reach.add(v); st.append(v)
    back = {t}
    st = [t]
    while st:
        u = st.pop()
        for v in in_s[u]:
            if v not in back:
                back.add(v); st.append(v)
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
            p.append(u); u = prev[u]
        return p[::-1]

    p1 = find_path(set())
    if p1 is None:
        return None
    if find_path(set(p1[1:-1])) is None:
        return None
    return (s, t)


def esu(nodes_list, adj_u, k, probs, rng, callback):
    idx = {v: i for i, v in enumerate(nodes_list)}
    for v in nodes_list:
        if probs[0] < 1.0 and rng.random() >= probs[0]:
            continue
        ext = {u for u in adj_u[v] if idx[u] > idx[v]}
        _extend([v], ext, v, adj_u, idx, k, probs, rng, callback)


def _extend(sub, ext, v, adj_u, idx, k, probs, rng, callback):
    if len(sub) == k:
        callback(sub); return
    d = len(sub)
    ext = set(ext)
    while ext:
        w = ext.pop()
        if probs[d] < 1.0 and rng.random() >= probs[d]:
            continue
        excl = set(sub) | {u for s in sub for u in adj_u[s]}
        newext = ext | {u for u in adj_u[w] if idx[u] > idx[v] and u not in excl}
        _extend(sub + [w], newext, v, adj_u, idx, k, probs, rng, callback)


def census(nodes, edges, k, probs, seed=1):
    out = collections.defaultdict(set)
    inn = collections.defaultdict(set)
    adj_u = collections.defaultdict(set)
    for a, b in edges:
        out[a].add(b); inn[b].add(a); adj_u[a].add(b); adj_u[b].add(a)
    nl = [v for v in sorted(nodes) if adj_u[v]]
    rng = random.Random(seed)
    stats = collections.Counter()
    comp = collections.Counter()
    classdiv = collections.Counter()

    def _cls(a, b):
        t = edges[(a, b)]
        if t == 'TF_target':
            return 'TF->TF' if nodes[b] == 'TF' else 'TF->Gene'
        if t == 'miRNA_target':
            return 'miRNA->TF' if nodes[b] == 'TF' else 'miRNA->Gene'
        return {'TF_miRNA': 'TF->miRNA', 'gene_gene': 'Gene->Gene',
                'miRNA_miRNA': 'miRNA-miRNA'}.get(t, t)

    def cb(sub):
        S = set(sub)
        os_ = {v: {u for u in out[v] if u in S} for v in S}
        is_ = {v: {u for u in inn[v] if u in S} for v in S}
        stats['subgraphs'] += 1
        if is_ffl(S, os_, is_):
            stats['ffl'] += 1
            tc = collections.Counter(nodes[v] for v in S)
            comp[(tc['TF'], tc['miRNA'], tc['Gene'])] += 1
            ec = frozenset(_cls(a, b) for a in S for b in os_[a])
            classdiv[len(ec)] += 1

    esu(nl, adj_u, k, probs, rng, cb)
    scale = 1.0
    for p in probs:
        scale /= p
    return stats, comp, classdiv, scale, len(nl)


# --------------------------------------------------------------------- driver
VARIANTS = collections.OrderedDict([
    ('V0_published',
     dict(string_mode='as_directed', mirna_thresh='threshold_10kb', mirna_tier='all',
          label='as published (STRING admitted as directed arcs, 10 kb clusters, all tiers)')),
    ('V1_noSTRING',
     dict(string_mode='exclude', mirna_thresh='threshold_10kb', mirna_tier='all',
          label='Methods-compliant: STRING associations excluded from the FFL graph')),
    ('V2_noSTRING_3kb',
     dict(string_mode='exclude', mirna_thresh='threshold_3kb', mirna_tier='all',
          label='STRING excluded, 3 kb polycistronic window')),
    ('V3_noSTRING_50kb',
     dict(string_mode='exclude', mirna_thresh='threshold_50kb', mirna_tier='all',
          label='STRING excluded, 50 kb polycistronic window')),
    ('V4_noSTRING_nomirmir',
     dict(string_mode='exclude', mirna_thresh=None, mirna_tier='all',
          drop_author_mirmir=True,
          label='STRING excluded and the miRNA-miRNA layer removed entirely')),
    ('V5_noSTRING_experimental',
     dict(string_mode='exclude', mirna_thresh='threshold_10kb', mirna_tier='experimental',
          label='STRING excluded; miRNA-target edges restricted to experimentally supported (strong+weak)')),
    ('V6_noSTRING_strongonly',
     dict(string_mode='exclude', mirna_thresh='threshold_10kb', mirna_tier='strong',
          label='STRING excluded; miRNA-target edges restricted to low-throughput validated (strong tier)')),
    ('V7_published_3kb',
     dict(string_mode='as_directed', mirna_thresh='threshold_3kb', mirna_tier='all',
          label='as published but 3 kb polycistronic window')),
    ('V8_published_50kb',
     dict(string_mode='as_directed', mirna_thresh='threshold_50kb', mirna_tier='all',
          label='as published but 50 kb polycistronic window')),
    ('V9_TRRUSTgenegene',
     dict(string_mode='trrust_only', mirna_thresh='threshold_10kb', mirna_tier='all',
          label='gene-gene layer restricted to directed TRRUST records only')),
])

PROBS = {3: [1.0, 1.0, 1.0],
         4: [1.0, 1.0, 1.0, 0.05],
         5: [1.0, 1.0, 0.6, 0.2, 0.1]}

if __name__ == '__main__':
    ns = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [3]
    seeds = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1]
    only = sys.argv[3].split(',') if len(sys.argv) > 3 else list(VARIANTS)

    grows, crows, clrows = [], [], []
    for vname in only:
        par = dict(VARIANTS[vname])
        label = par.pop('label')
        nodes, edges, meta = build_graph(**par)
        grows.append(dict(variant=vname, label=label, **{k: v for k, v in meta.items()
                                                         if k != 'by_type'},
                          **{f'edges_{k}': v for k, v in meta['by_type'].items()}))
        print(f'\n=== {vname}: {meta["n_edges"]} arcs  {meta["by_type"]}')
        for n in ns:
            for sd in ([1] if n == 3 else seeds):
                t0 = time.time()
                st, comp, cdiv, scale, nconn = census(nodes, edges, n, PROBS[n], sd)
                el = time.time() - t0
                crows.append(dict(variant=vname, label=label, n=n, seed=sd,
                                  probs=';'.join(str(p) for p in PROBS[n]),
                                  connected_nodes=nconn, arcs=meta['n_edges'],
                                  sampled_subgraphs=st['subgraphs'], sampled_ffl=st['ffl'],
                                  scale=scale,
                                  est_subgraphs=round(st['subgraphs'] * scale),
                                  est_ffl=round(st['ffl'] * scale),
                                  exhaustive=(n == 3), seconds=round(el, 1)))
                for nc, v in cdiv.items():
                    clrows.append(dict(variant=vname, n=n, seed=sd, n_edge_classes=nc,
                                       sampled=v, est_modules=round(v * scale)))
                print(f'   n={n} seed={sd}: sampled {st["subgraphs"]:,} sub, '
                      f'{st["ffl"]:,} FFL -> est {st["ffl"]*scale:,.0f}   ({el:.0f}s)  '
                      f'max classes {max(cdiv) if cdiv else 0}')

    def wr(path, rows):
        if not rows:
            return
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print('wrote', path)

    tag = 'n' + ''.join(str(x) for x in ns)
    wr(f'{OUT}/sensitivity_graphs.csv', grows)
    wr(f'{OUT}/sensitivity_census_{tag}.csv', crows)
    wr(f'{OUT}/sensitivity_census_classes_{tag}.csv', clrows)
