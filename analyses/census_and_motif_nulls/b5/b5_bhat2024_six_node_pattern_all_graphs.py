#!/usr/bin/env python3
"""
B5: instances of the six-node FFL of Bhat et al. (2024, Current Bioinformatics 19:73-90,
doi:10.2174/1574893618666230731164002), counted in the breast-cancer network.

Pattern = the COMPULSORY (solid) edges of their Fig. 1d; their dotted edges are optional and ignored.
Roles: TF1, TF2 (TF nodes), M1, M2 (miRNA nodes), G1, G2 (Gene nodes), all distinct.
  M1 -> G1 and M2 -> G1           miRNA->Gene arcs (both miRNAs regulate the same gene)
  M1 - M2                          miRNA-miRNA link, either direction
  TF1 - TF2                        TF->TF arc (paper's class TF->TF: TF_target layer), either direction
  TF2 -> G2                        TF->Gene arc (paper's class TF->Gene)
  G1 - G2                          arc between the two Gene nodes, either direction (gene-gene layer)
  TF1 / M1 link, by FFL class:
     miRNA-FFL   M1 -> TF1         (arc in the uncontracted network)
     TF-FFL      TF1 -> M1
     composite   TF1 <-> M1        (both)
Counts are given per role assignment (instances) and per distinct six-node set, with the numbers of
distinct nodes and pairs used, as in their Table 3.  Run on two graphs exported by
nolegacy/export_graphs.py: the Table 3 census graph (incl. 28 legacy miRNA-miRNA arcs) and the graph
the Methods describe (miRNA-miRNA = 10 kb polycistronic pairs only).
"""
import sys, json, itertools, collections

CLS = ['TF->Gene', 'TF->TF', 'TF->miRNA', 'miRNA->Gene', 'miRNA->TF', 'Gene->Gene', 'miRNA-miRNA']


def load(tag):
    L = open(f'../graphs/census_{tag}_R4.txt').read().split('\n')
    nv, ne = map(int, L[0].split())
    names = [l.split('\t')[0] for l in L[1:1 + nv]]; types = [l.split('\t')[1] for l in L[1:1 + nv]]
    cls = {}
    for l in L[1 + nv:1 + nv + ne]:
        u, v, c = map(int, l.split()); cls[(u, v)] = CLS[c]
    O = open(f'../graphs/census_{tag}_orig.txt').read().split('\n')
    no = int(O[0].split()[1])
    for l in O[2:2 + no]:                          # arcs removed by TF<->miRNA contraction are miRNA->TF
        u, v = map(int, l.split()[:2])
        if (u, v) not in cls:
            assert types[u] == 'miRNA' and types[v] == 'TF'
            cls[(u, v)] = 'miRNA->TF'
    return names, types, cls


def count(tag):
    names, types, cls = load(tag)
    T = lambda t: [i for i, x in enumerate(types) if x == t]
    arcs = collections.defaultdict(set)
    for (u, v), c in cls.items():
        arcs[c].add((u, v))
    und = lambda c: {frozenset(e) for e in arcs[c]}
    mm = und('miRNA-miRNA'); tt = und('TF->TF')
    gg = {frozenset((u, v)) for (u, v) in arcs['Gene->Gene'] if types[u] == 'Gene' and types[v] == 'Gene'}
    mir_of_gene = collections.defaultdict(set)
    for (m, g) in arcs['miRNA->Gene']:
        if types[g] == 'Gene':
            mir_of_gene[g].add(m)
    tf2gene = collections.defaultdict(set)
    for (t, g) in arcs['TF->Gene']:
        if types[g] == 'Gene':
            tf2gene[t].add(g)
    m2tf = collections.defaultdict(set); tf2m = collections.defaultdict(set)
    for (m, t) in arcs['miRNA->TF']: m2tf[m].add(t)
    for (t, m) in arcs['TF->miRNA']: tf2m[m].add(t)
    tf_nb = collections.defaultdict(set)
    for e in tt:
        a, b = tuple(e); tf_nb[a].add(b); tf_nb[b].add(a)
    g_nb = collections.defaultdict(set)
    for e in gg:
        a, b = tuple(e); g_nb[a].add(b); g_nb[b].add(a)

    res = {}
    for label, tf1_of in (('miRNA_FFL (M1->TF1)', lambda m: m2tf[m]),
                          ('TF_FFL (TF1->M1)', lambda m: tf2m[m]),
                          ('composite_FFL (TF1<->M1)', lambda m: m2tf[m] & tf2m[m])):
        inst = 0; sets = set(); used = collections.defaultdict(set); ex = []
        for g1 in T('Gene'):
            ms = mir_of_gene[g1]
            for m1, m2 in itertools.permutations(sorted(ms), 2):
                if frozenset((m1, m2)) not in mm:
                    continue
                for tf1 in tf1_of(m1):
                    for tf2 in tf_nb[tf1]:
                        for g2 in g_nb[g1] & tf2gene[tf2]:
                            if g2 == g1:
                                continue
                            inst += 1
                            sets.add(frozenset((tf1, tf2, m1, m2, g1, g2)))
                            for k, v in (('miRNA', m1), ('miRNA', m2), ('TF', tf1), ('TF', tf2), ('Gene', g1), ('Gene', g2)):
                                used[k].add(v)
                            used['miRNA-miRNA'].add(frozenset((m1, m2))); used['TF-TF'].add(frozenset((tf1, tf2)))
                            used['Gene-Gene'].add(frozenset((g1, g2))); used['miRNA-Gene'].update({(m1, g1), (m2, g1)})
                            used['TF-miRNA link'].add((tf1, m1)); used['TF-Gene'].add((tf2, g2))
                            if len(ex) < 3:
                                ex.append(f"TF1={names[tf1]} TF2={names[tf2]} M1={names[m1]} M2={names[m2]} G1={names[g1]} G2={names[g2]}")
        res[label] = dict(instances=inst, distinct_six_node_sets=len(sets),
                          **{f'n_distinct_{k}': len(v) for k, v in used.items()}, examples=ex)
    return dict(graph=tag, n_arcs_contracted=sum(1 for _ in open(f'../graphs/census_{tag}_R4.txt')) - 1 - len(names),
                n_miRNA_miRNA_links=len(mm), n_TF_TF_links=len(tt), n_gene_gene_links_between_Gene_nodes=len(gg),
                classes=res)


if __name__ == '__main__':
    out = [count(t) for t in ('table3', 'nolegacy', 'table3_nostring', 'nolegacy_nostring')]
    json.dump(out, open('B5_bhat2024_six_node_pattern.json', 'w'), indent=1)
    for o in out:
        print(f"\n== graph {o['graph']}: {o['n_arcs_contracted']} census arcs; links: miRNA-miRNA {o['n_miRNA_miRNA_links']}, "
              f"TF-TF {o['n_TF_TF_links']}, Gene-Gene {o['n_gene_gene_links_between_Gene_nodes']}")
        for k, v in o['classes'].items():
            print(f"  {k:28s} instances {v['instances']:>8,}  distinct 6-node sets {v['distinct_six_node_sets']:>7,}  "
                  + '  '.join(f"{kk[11:]}={vv}" for kk, vv in v.items() if kk.startswith('n_distinct')))
            for e in v['examples']:
                print('      e.g.', e)
