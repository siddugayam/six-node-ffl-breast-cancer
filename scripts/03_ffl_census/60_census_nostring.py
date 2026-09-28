#!/usr/bin/env python3
"""
Sensitivity re-run of the RAND-ESU n-node FFL census (scripts/03_ffl_census/03_ffl_census.py)
with STRING co-functional associations EXCLUDED from the gene-gene layer.

Motivation. build_graph() in 03_ffl_census.py adds every row of data/layer_gene_gene.tsv
as a directed gene->gene arc. 941 of those 1,638 rows are STRING undirected functional
associations (directed == FALSE). After the TF_target layer is added first, ALL 833
gene->gene arcs that actually enter the census graph are STRING associations. Treating a
, so the census must be reported with a STRING-free sensitivity analysis.

Usage: 60_census_nostring.py <n> <probs,comma> <seed>
Writes results/v5/ffl_census_nostring_k<n>_s<seed>.json
"""
import csv, collections, importlib.util, json, os, random, sys, time

REV = '/path/to/revision'
spec = importlib.util.spec_from_file_location('c3', f'{REV}/scripts/03_ffl_census/03_ffl_census.py')
C3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(C3)

def build_graph_nostring(mirna_thresh='threshold_10kb'):
    """Identical to C3.build_graph(extra_layers=True) except that layer_gene_gene.tsv
    rows with directed != TRUE (i.e. the STRING associations) are dropped."""
    nodes = {r['name']: r['type'] for r in
             csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}
    edges = {}
    for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t'):
        edges[(r['source'], r['target'])] = r['edge_type']

    def add(path, sc, tc, et, undirected=False, filt=None):
        if not os.path.exists(path): return 0
        n = 0
        for r in csv.DictReader(open(path), delimiter='\t'):
            if filt and not filt(r): continue
            a, b = r.get(sc), r.get(tc)
            if a not in nodes or b not in nodes or a == b: continue
            for e in ([(a, b), (b, a)] if undirected else [(a, b)]):
                if e not in edges: edges[e] = et; n += 1
        return n

    add(f'{REV}/data/layer_TF_target.tsv', 'source', 'target', 'TF_target')
    add(f'{REV}/data/layer_gene_gene.tsv', 'source', 'target', 'gene_gene',
        filt=lambda r: str(r.get('directed', 'TRUE')).upper() in ('TRUE', '1'))
    add(f'{REV}/data/layer_miRNA_miRNA.tsv', 'miRNA_1', 'miRNA_2', 'miRNA_miRNA',
        undirected=True,
        filt=lambda r: str(r.get(mirna_thresh, 'TRUE')).upper() in ('TRUE', '1'))
    nrecip = 0
    for (a, b) in list(edges):
        if nodes.get(a) == 'TF' and nodes.get(b) == 'miRNA' and (b, a) in edges:
            if edges.pop((b, a), None) is not None: nrecip += 1
    return nodes, edges, nrecip

C3.build_graph = lambda extra_layers=True, mirna_thresh='threshold_10kb': build_graph_nostring(mirna_thresh)

if __name__ == '__main__':
    k = int(sys.argv[1])
    probs = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1.0] * k
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    t0 = time.time()
    st, comp, found, nodes, classdiv, classsets = C3.census(k, probs, seed=seed)
    scale = 1.0
    for p_ in probs: scale /= p_
    res = {'k': k, 'seed': seed, 'probs': probs,
           'graph': 'canonical + TRRUST TF-target + TransmiR TF-miRNA + 10kb polycistron; STRING EXCLUDED',
           'n_arcs': 7198,
           'sampled_subgraphs': st['subgraphs'], 'sampled_ffl': st['ffl'], 'scale': scale,
           'est_subgraphs': st['subgraphs'] * scale, 'est_ffl': st['ffl'] * scale,
           'class_diversity': {str(a): b * scale for a, b in classdiv.items()},
           'class_sets': {'+'.join(a): b * scale for a, b in classsets.items()},
           'composition': {f"TF{a}_miR{b}_G{c}": v * scale for (a, b, c), v in comp.items()},
           'runtime_s': round(time.time() - t0, 1)}
    out = f'{REV}/results/v5/ffl_census_nostring_k{k}_s{seed}.json'
    json.dump(res, open(out, 'w'), indent=1)
    print(f"k={k} probs={probs} seed={seed} sampled={st['subgraphs']} ffl_sampled={st['ffl']} "
          f"est_ffl={res['est_ffl']:,.0f} runtime={res['runtime_s']}s")
    print('class diversity (n_classes -> est modules):',
          {a: round(b) for a, b in sorted(res['class_diversity'].items())})
    print('wrote', out)
