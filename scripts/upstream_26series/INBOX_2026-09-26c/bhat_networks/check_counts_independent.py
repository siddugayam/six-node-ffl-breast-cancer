#!/usr/bin/env python3
"""Independent re-count of bhat_table3_counts.csv (instances and distinct node sets).
Differs from build_bhat_networks.py in both steps:
  graph   rebuilt from the project data with the census script's own build_graph() (Option 1 script
          scripts/v7/03_ffl_census_nolegacymirna.py, DROP_LEGACY_MIRNA=1), edge types taken from its labels
          (TF_target, miRNA_target, TF_miRNA, gene_gene, miRNA_miRNA), STRING arcs removed for the sensitivity
          graph by the rule of INBOX_2026-09-26b/graphs/make_graphs.py
  loops   start from the target side (G1, then G2, then the TFs, then the miRNAs) and test every compulsory
          edge with explicit predicates"""
import os, sys, csv, collections, importlib.util
sys.dont_write_bytecode = True
REV = '/path/to/revision'
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['DROP_LEGACY_MIRNA'] = '1'
spec = importlib.util.spec_from_file_location('m', f'{REV}/scripts/v7/03_ffl_census_nolegacymirna.py')
m = importlib.util.module_from_spec(spec); sys.argv = ['x', '3']; spec.loader.exec_module(m)


def graph(nostring):
    nodes, edges, nrecip = m.build_graph(True)
    assert nrecip == 1223
    dep = {(r['source'], r['target']) for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t')}
    for (a, b) in dep:            # restore miRNA->TF arcs removed by the reciprocal-pair contraction
        if nodes.get(a) == 'miRNA' and nodes.get(b) == 'TF' and (b, a) in edges and (a, b) not in edges:
            edges[(a, b)] = 'miRNA_target'
    if nostring:
        sp = {(r['source'], r['target']) for r in csv.DictReader(open(f'{REV}/data/layer_gene_gene.tsv'), delimiter='\t')
              if r['evidence'] == 'STRING'}
        drop = [e for e, t in edges.items() if t == 'gene_gene' and e in sp and e not in dep]
        assert len(drop) == 833
        for e in drop: del edges[e]
    return nodes, edges


def count(nodes, edges):
    typ = nodes.get
    arc = lambda a, b, t: edges.get((a, b)) == t
    MT = lambda mm, tt: arc(mm, tt, 'miRNA_target')
    TM = lambda tt, mm: arc(tt, mm, 'TF_miRNA')
    TT = lambda a, b: arc(a, b, 'TF_target')
    Gs = sorted(v for v in nodes if typ(v) == 'Gene'); Ts = sorted(v for v in nodes if typ(v) == 'TF')
    Ms = sorted(v for v in nodes if typ(v) == 'miRNA')
    mir_of = {g: [x for x in Ms if arc(x, g, 'miRNA_target')] for g in Gs}
    tf_of = {g: [x for x in Ts if arc(x, g, 'TF_target')] for g in Gs}
    gg = lambda a, b: a != b and (arc(a, b, 'gene_gene') or arc(b, a, 'gene_gene'))
    mm = lambda a, b: a != b and (arc(a, b, 'miRNA_miRNA') or arc(b, a, 'miRNA_miRNA'))
    ggn = {g: [h for h in Gs if gg(g, h)] for g in Gs}
    cls = lambda m1, t1: (MT(m1, t1), TM(t1, m1))
    res = collections.defaultdict(set)

    def add(n, inst, mir, tf):
        if mir: res[(n, 'miRNA_FFL')].add(inst)
        if tf: res[(n, 'TF_FFL')].add(inst)
        if mir and tf: res[(n, 'Composite_FFL')].add(inst)

    for g1 in Gs:
        for t1 in tf_of[g1]:
            for m1 in mir_of[g1]:
                add(3, (m1, t1, g1), *cls(m1, t1))
        for g2 in ggn[g1]:
            for m1 in mir_of[g1]:
                m2s = [x for x in mir_of[g1] if mm(m1, x)]
                for t1 in tf_of[g2]:
                    a, b = cls(m1, t1)
                    add(4, (m1, t1, g1, g2), a, b)
                    for m2 in m2s: add(5, (m1, m2, t1, g1, g2), a, b)
            for t2 in tf_of[g2]:
                for t1 in Ts:
                    if t1 == t2 or not (TT(t1, t2) or TT(t2, t1)): continue
                    for m1 in mir_of[g1]:
                        a, b = cls(m1, t1)
                        if not (a or b): continue
                        for m2 in mir_of[g1]:
                            if not mm(m1, m2): continue
                            inst = (m1, m2, t1, t2, g1, g2)
                            add(6, inst, a and TT(t1, t2), b and TT(t2, t1))
                            if a: res[(6, 'miRNA_FFL relaxed')].add(inst)
                            if b: res[(6, 'TF_FFL relaxed')].add(inst)
                            if a and b: res[(6, 'Composite_FFL relaxed')].add(inst)
    for n in (3, 4, 5, 6):
        res[(n, 'Unique_FFL')] = res[(n, 'miRNA_FFL')] | res[(n, 'TF_FFL')]
    return res


if __name__ == '__main__':
    ref = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'bhat_table3_counts.csv'))):
        k = (r['graph'], int(r['nodes']), r['FFL_type'] + (' relaxed' if r['TF_TF_rule'].startswith('either') else ''))
        ref[k] = (int(r['No_of_FFL_instances']), int(r['distinct_node_sets']))
    bad = 0
    for tag, ns in (('nolegacy', False), ('nolegacy_nostring', True)):
        res = count(*graph(ns))
        for (n, k), s in sorted(res.items()):
            got = (len(s), len({frozenset(i) for i in s})); exp = ref[(tag, n, k)]
            bad += got != exp
            print(f"{tag:18s} {n}-node {k:22s} instances {got[0]:>7,} sets {got[1]:>7,}   "
                  f"{'MATCH' if got == exp else 'DIFFERS from ' + str(exp)}")
    print('ALL MATCH' if bad == 0 else f'{bad} DIFFER')
