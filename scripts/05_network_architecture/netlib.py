#!/usr/bin/env python3
"""Shared network utilities for the v3 systems-pharmacology / network-control analysis.

Primary object: the canonical manuscript network
  587 nodes, 6,859 directed signed edges  (data/canonical_edges.tsv)

Sign conventions
  miRNA_target  : -1  (repression; mechanistically certain)
  TF_target     : +1 in the published file for ALL edges; the curated sign comes from
                  TRRUST mode in data/layer_TF_target.tsv (Activation +1, Repression -1,
                  Unknown/Ambiguous -> None)
  TF_miRNA      : TransmiR mode in data/layer_TF_miRNA.tsv
  miRNA_miRNA   : +1 co-transcription (genomic cluster; direction arbitrary)
  gene_gene     : +1 (single STRING edge COL1A1-COL3A1)

FFL core definition (identical to scripts/03_ffl_census/02x_ffl_core_crosscheck.py, which reproduces the
established count of 1,649 unique cores):
  R->M, R->T, M->T all present; type(T) != miRNA;
  (type(R)==miRNA and type(M)==TF) or (type(R)==TF and type(M)==miRNA)
  Composite-FFL if (M,R) is also an edge, else miRNA-FFL / TF-FFL.
  Unique count = composites/2 + miRNA-FFL + TF-FFL
"""
import csv
import os
import collections

REV = '/path/to/revision'
DATA = os.path.join(REV, 'data')
RES = os.path.join(REV, 'results', 'v3')
FIG = os.path.join(REV, 'figures', 'v3')
LOG = os.path.join(REV, 'logs', 'v3')


def load_nodes():
    nt = {}
    with open(os.path.join(DATA, 'canonical_nodes.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            nt[r['name']] = r['type']
    return nt


def load_edges():
    """Returns list of dicts with source,target,edge_type,sign_published,sign_curated."""
    curated = {}
    for fn, sk, tk, mk in [('layer_TF_target.tsv', 'source', 'target', 'mode'),
                           ('layer_TF_miRNA.tsv', 'source', 'target', 'mode')]:
        p = os.path.join(DATA, fn)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                s = {'Activation': 1, 'Repression': -1}.get((r.get(mk) or '').strip())
                if s is not None:
                    curated[(r[sk], r[tk])] = s
    E = []
    with open(os.path.join(DATA, 'canonical_edges.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            e = (r['source'], r['target'])
            et = r['edge_type']
            if et == 'miRNA_target':
                sc = -1
            elif et in ('miRNA_miRNA', 'gene_gene'):
                sc = 1
            else:
                sc = curated.get(e)          # None when TRRUST/TransmiR say Unknown/Ambiguous
            E.append(dict(source=r['source'], target=r['target'], edge_type=et,
                          sign_published=int(r['sign']), sign_curated=sc,
                          source_type=r['source_type'], target_type=r['target_type']))
    return E


def adjacency(E):
    out = collections.defaultdict(set)
    inn = collections.defaultdict(set)
    for e in E:
        out[e['source']].add(e['target'])
        inn[e['target']].add(e['source'])
    return out, inn


# ---------------------------------------------------------------- FFL cores
def ffl_cores(edge_pairs, node_type):
    """edge_pairs: set of (u,v). Returns list of (R,M,T,class)."""
    out = collections.defaultdict(set)
    for a, b in edge_pairs:
        out[a].add(b)
    cores = []
    for R, tg in out.items():
        tR = node_type.get(R, 'Gene')
        for M in tg:
            tM = node_type.get(M, 'Gene')
            if not ((tR == 'miRNA' and tM == 'TF') or (tR == 'TF' and tM == 'miRNA')):
                continue
            comp = 'Composite-FFL' if (M, R) in edge_pairs else (
                'miRNA-FFL' if tR == 'miRNA' else 'TF-FFL')
            for T in (tg & out.get(M, set())):
                if T == R or T == M:
                    continue
                if node_type.get(T, 'Gene') == 'miRNA':
                    continue
                cores.append((R, M, T, comp))
    return cores


def ffl_unique_count(cores):
    """Composite cores are counted once per reciprocal (R,M) pair."""
    n_comp = sum(1 for c in cores if c[3] == 'Composite-FFL')
    n_oth = len(cores) - n_comp
    return n_comp // 2 + n_oth


def ffl_count_fast(edge_pairs, node_type):
    return ffl_unique_count(ffl_cores(edge_pairs, node_type))


# ---------------------------------------------------------------- module sets
MIR29 = ['hsa-miR-29a', 'hsa-miR-29b', 'hsa-miR-29c']
COLLAGEN = ['COL1A1', 'COL3A1']


def mir29_module(edge_pairs, node_type):
    """miR-29 family + the collagen targets + every node on a direct edge between them."""
    core = set(MIR29) | set(COLLAGEN)
    return core


if __name__ == '__main__':
    nt = load_nodes()
    E = load_edges()
    P = {(e['source'], e['target']) for e in E}
    print('nodes', len(nt), 'edges', len(E))
    c = ffl_cores(P, nt)
    cc = collections.Counter(x[3] for x in c)
    print(dict(cc), 'unique =', ffl_unique_count(c))
