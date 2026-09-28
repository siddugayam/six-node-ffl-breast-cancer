#!/usr/bin/env python3
"""
3-to-6-node FFL networks of Bhat et al. (2024, Curr Bioinform 19:73-90,
doi:10.2174/1574893618666230731164002), rebuilt on the breast-cancer census graph.
Reads the analyses/census_and_motif_nulls graph exports read-only; writes only into this folder.

Graphs
  nolegacy           Option 1 (adopted 2026-09-26): census graph without the 28 legacy miRNA-miRNA arcs
  nolegacy_nostring  sensitivity: the same without the 833 STRING gene-gene arcs

Patterns = the COMPULSORY (solid) edges of their Fig. 1; dotted (optional) edges are ignored, so extra
edges between the six nodes are allowed.  Roles are distinct nodes; targets G are Gene-typed nodes.
  3-node  M1-T1 link, M1->G1, T1->G1
  4-node  M1-T1 link, M1->G1, T1->G2, G1-G2
  5-node  M1-T1 link, M1->G1, M2->G1, M1-M2, T1->G2, G1-G2
  6-node  M1-T1 link, M1->G1, M2->G1, M1-M2, T1-T2 link, T2->G2, G1-G2
Class-specific links (Fig. 1 arrowheads; composite = double-headed arrows):
  miRNA-FFL   M1->T1                 6-node also T1->T2   (chain M1 -> T1 -> T2 -> G2)
  TF-FFL      T1->M1                 6-node also T2->T1   (chain T2 -> T1 -> M1 -> G1)
  composite   M1<->T1                6-node also T1<->T2
so composite = miRNA-FFL AND TF-FFL for the same role assignment, and unique = miRNA-FFL OR TF-FFL.
Undirected links: G1-G2 (co-expression line), M1-M2 (synergy, ball-headed).  Edge sources in this graph:
M->T and M->G = miRNA_target; T->M = TF_miRNA; T->T and T->G = TF_target; G-G = gene_gene layer between
two Gene nodes; M-M = 10 kb polycistronic pairs.
Reconciliation row: '6-node, TF-TF either direction' is the B5 definition sent on 2026-09-26, except that
TF1 and TF2 must be TF-typed nodes here (B5 took the TF-TF edge class, which admits 21 arcs from Gene-typed
regulators).
Counts: instances = role assignments (a node set can satisfy the pattern in more than one way);
distinct node sets also given.  Table 3 columns = distinct nodes and distinct compulsory pairs used.
"""
import sys, os, gzip, json, itertools, collections
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
GRAPHS = os.path.join(HERE, '..', '..', 'analyses/census_and_motif_nulls', 'graphs')
CLS = ['TF->Gene', 'TF->TF', 'TF->miRNA', 'miRNA->Gene', 'miRNA->TF', 'Gene->Gene', 'miRNA-miRNA']
ROLES = {3: ('M1', 'T1', 'G1'), 4: ('M1', 'T1', 'G1', 'G2'), 5: ('M1', 'M2', 'T1', 'G1', 'G2'),
         6: ('M1', 'M2', 'T1', 'T2', 'G1', 'G2')}
CLASSES = ('miRNA_FFL', 'TF_FFL', 'Composite_FFL')
COLS = ['miRNA', 'TF', 'Gene', 'Gene-Gene', 'miRNA-miRNA', 'TF-TF', 'miRNA-Gene', 'miRNA-TF', 'TF-miRNA', 'TF-Gene']


def load(tag):
    L = open(os.path.join(GRAPHS, f'census_{tag}_R4.txt')).read().split('\n')
    nv, ne = map(int, L[0].split())
    names = [l.split('\t')[0] for l in L[1:1 + nv]]; types = [l.split('\t')[1] for l in L[1:1 + nv]]
    cls = {}
    for l in L[1 + nv:1 + nv + ne]:
        u, v, c = map(int, l.split()); cls[(u, v)] = CLS[c]
    O = open(os.path.join(GRAPHS, f'census_{tag}_orig.txt')).read().split('\n')
    no = int(O[0].split()[1])
    for l in O[2:2 + no]:                 # arcs removed by the TF<->miRNA contraction are miRNA->TF
        u, v = map(int, l.split()[:2])
        if (u, v) not in cls:
            assert types[u] == 'miRNA' and types[v] == 'TF'
            cls[(u, v)] = 'miRNA->TF'
    assert len(cls) == no
    return names, types, cls


class G:
    def __init__(self, tag):
        self.tag = tag
        self.names, self.types, cls = load(tag)
        ty = self.types
        self.MG = collections.defaultdict(set); self.TG = collections.defaultdict(set)
        self.MT = set(); self.TM = set(); self.TT = set()
        self.GGnb = collections.defaultdict(set); self.MMnb = collections.defaultdict(set)
        # TF roles take TF-typed nodes only: 36 TF_target arcs start at Gene-typed regulators (e.g. NR4A1,
        # MYBL2, STAT5A), which the census edge classes label TF->TF / TF->Gene
        for (u, v), c in cls.items():
            if c == 'miRNA->Gene' and ty[v] == 'Gene': self.MG[u].add(v)
            elif c == 'TF->Gene' and ty[u] == 'TF' and ty[v] == 'Gene': self.TG[u].add(v)
            elif c == 'miRNA->TF': self.MT.add((u, v))
            elif c == 'TF->miRNA': self.TM.add((u, v))
            elif c == 'TF->TF' and ty[u] == 'TF' and ty[v] == 'TF' and u != v: self.TT.add((u, v))
            elif c == 'Gene->Gene' and ty[u] == 'Gene' and ty[v] == 'Gene' and u != v:
                self.GGnb[u].add(v); self.GGnb[v].add(u)
            elif c == 'miRNA-miRNA' and u != v:
                self.MMnb[u].add(v); self.MMnb[v].add(u)
        self.TTnb_out = collections.defaultdict(set); self.TTnb_in = collections.defaultdict(set)
        for (a, b) in self.TT: self.TTnb_out[a].add(b); self.TTnb_in[b].add(a)
        self.link_pairs = {  # (M, T) pairs carrying the class-specific M1-T1 link
            'miRNA_FFL': sorted(self.MT),
            'TF_FFL': sorted((m, t) for (t, m) in self.TM),
            'Composite_FFL': sorted(p for p in self.MT if (p[1], p[0]) in self.TM)}

    def tf2(self, klass, t1, relaxed=False):
        if relaxed: return self.TTnb_out[t1] | self.TTnb_in[t1]
        if klass == 'miRNA_FFL': return self.TTnb_out[t1]              # T1 -> T2
        if klass == 'TF_FFL': return self.TTnb_in[t1]                  # T2 -> T1
        return self.TTnb_out[t1] & self.TTnb_in[t1]                    # T1 <-> T2

    def enumerate(self, n, klass, relaxed=False):
        out = []
        for (m1, t1) in self.link_pairs[klass]:
            if n == 3:
                for g1 in sorted(self.MG[m1] & self.TG[t1]): out.append((m1, t1, g1))
                continue
            for g1 in sorted(self.MG[m1]):
                if n == 4:
                    for g2 in sorted(self.GGnb[g1] & self.TG[t1]): out.append((m1, t1, g1, g2))
                    continue
                m2s = sorted(m2 for m2 in self.MMnb[m1] if g1 in self.MG[m2])
                if not m2s: continue
                if n == 5:
                    g2s = sorted(self.GGnb[g1] & self.TG[t1])
                    for m2 in m2s:
                        for g2 in g2s: out.append((m1, m2, t1, g1, g2))
                    continue
                for t2 in sorted(self.tf2(klass, t1, relaxed)):
                    g2s = sorted(self.GGnb[g1] & self.TG[t2])
                    for m2 in m2s:
                        for g2 in g2s: out.append((m1, m2, t1, t2, g1, g2))
        return out

    def edges_of(self, n, klass, inst, relaxed=False):
        """compulsory edges of one instance as (source, label, target); undirected links sorted by name."""
        r = dict(zip(ROLES[n], inst)); nm = self.names; E = []

        def und(a, b, lab):
            x, y = sorted((nm[a], nm[b])); E.append((x, lab, y))
        if klass in ('miRNA_FFL', 'Composite_FFL'): E.append((nm[r['M1']], 'miRNA-TF', nm[r['T1']]))
        if klass in ('TF_FFL', 'Composite_FFL'): E.append((nm[r['T1']], 'TF-miRNA', nm[r['M1']]))
        E.append((nm[r['M1']], 'miRNA-Gene', nm[r['G1']]))
        if n == 3: E.append((nm[r['T1']], 'TF-Gene', nm[r['G1']])); return E
        und(r['G1'], r['G2'], 'Gene-Gene')
        if n >= 5:
            E.append((nm[r['M2']], 'miRNA-Gene', nm[r['G1']])); und(r['M1'], r['M2'], 'miRNA-miRNA')
        if n <= 5: E.append((nm[r['T1']], 'TF-Gene', nm[r['G2']])); return E
        E.append((nm[r['T2']], 'TF-Gene', nm[r['G2']]))
        t1, t2 = r['T1'], r['T2']
        if relaxed:
            if (t1, t2) in self.TT: E.append((nm[t1], 'TF-TF', nm[t2]))
            if (t2, t1) in self.TT: E.append((nm[t2], 'TF-TF', nm[t1]))
        else:
            if klass in ('miRNA_FFL', 'Composite_FFL'): E.append((nm[t1], 'TF-TF', nm[t2]))
            if klass in ('TF_FFL', 'Composite_FFL'): E.append((nm[t2], 'TF-TF', nm[t1]))
        return E


def table_row(g, n, klass, insts, relaxed=False):
    nodes = {'miRNA': set(), 'TF': set(), 'Gene': set()}; pairs = collections.defaultdict(set)
    for inst in insts:
        r = dict(zip(ROLES[n], inst))
        for k, v in r.items(): nodes[{'M': 'miRNA', 'T': 'TF', 'G': 'Gene'}[k[0]]].add(v)
        for (s, lab, t) in g.edges_of(n, klass, inst, relaxed):
            key = frozenset((s, t)) if lab in ('Gene-Gene', 'miRNA-miRNA', 'TF-TF') else (s, t)
            pairs[lab].add(key)
    row = {'miRNA': len(nodes['miRNA']), 'TF': len(nodes['TF']), 'Gene': len(nodes['Gene'])}
    for c in COLS[3:]:
        applicable = {'Gene-Gene': n >= 4, 'miRNA-miRNA': n >= 5, 'TF-TF': n == 6,
                      'miRNA-TF': klass != 'TF_FFL', 'TF-miRNA': klass != 'miRNA_FFL'}.get(c, True)
        row[c] = len(pairs[c]) if applicable else 'NA'
    row['No_of_FFL_instances'] = len(insts)
    row['distinct_node_sets'] = len({frozenset(i) for i in insts})
    return row


def main():
    os.makedirs(os.path.join(HERE, 'sif'), exist_ok=True); os.makedirs(os.path.join(HERE, 'instances'), exist_ok=True)
    summary = []
    for tag in ('nolegacy', 'nolegacy_nostring'):
        g = G(tag); nm = g.names
        with open(os.path.join(HERE, 'sif', f'{tag}_node_attributes.tsv'), 'w') as f:
            f.write('node\ttype\n')
            for i in sorted(range(len(nm)), key=lambda i: nm[i]): f.write(f'{nm[i]}\t{g.types[i]}\n')
        print(f'\n== {tag}: M->T {len(g.MT)}, T->M {len(g.TM)}, reciprocal M<->T {len(g.link_pairs["Composite_FFL"])}, '
              f'T->T {len(g.TT)}, reciprocal T<->T pairs {sum(1 for (a, b) in g.TT if (b, a) in g.TT) // 2}, '
              f'Gene-Gene links {sum(len(s) for s in g.GGnb.values()) // 2}, miRNA-miRNA links {sum(len(s) for s in g.MMnb.values()) // 2}')
        merged_comp = set()
        for n in (3, 4, 5, 6):
            sets = {}
            for klass in CLASSES:
                insts = g.enumerate(n, klass)
                sets[klass] = set(insts)
                assert len(sets[klass]) == len(insts)
                row = table_row(g, n, klass, insts)
                summary.append(dict(graph=tag, nodes=n, FFL_type=klass, TF_TF_rule='Fig. 1d' if n == 6 else '', **row))
                E = sorted({e for inst in insts for e in g.edges_of(n, klass, inst)})
                with open(os.path.join(HERE, 'sif', f'{tag}_{n}node_{klass}.sif'), 'w') as f:
                    for (s, lab, t) in E: f.write(f'{s}\t{lab}\t{t}\n')
                if klass == 'Composite_FFL': merged_comp |= set(E)
                with gzip.open(os.path.join(HERE, 'instances', f'{tag}_{n}node_{klass}.tsv.gz'), 'wt') as f:
                    f.write('\t'.join(ROLES[n]) + '\n')
                    for inst in insts: f.write('\t'.join(nm[i] for i in inst) + '\n')
            assert sets['Composite_FFL'] == sets['miRNA_FFL'] & sets['TF_FFL']
            uni = sets['miRNA_FFL'] | sets['TF_FFL']
            summary.append(dict(graph=tag, nodes=n, FFL_type='Unique_FFL', TF_TF_rule='Fig. 1d' if n == 6 else '',
                                **{c: 'NA' for c in COLS}, No_of_FFL_instances=len(uni),
                                distinct_node_sets=len({frozenset(i) for i in uni})))
            for r in summary[-4:]:
                print(f"  {n}-node {r['FFL_type']:14s} " + ' '.join(f"{c}={r[c]}" for c in COLS)
                      + f"  FFL={r['No_of_FFL_instances']:,}  sets={r['distinct_node_sets']:,}")
        # reconciliation with B5 (2026-09-26): TF-TF link in either direction
        rs = {}
        for klass in CLASSES:
            insts = g.enumerate(6, klass, relaxed=True); rs[klass] = set(insts)
            row = table_row(g, 6, klass, insts, relaxed=True)
            summary.append(dict(graph=tag, nodes=6, FFL_type=klass, TF_TF_rule='either direction (B5)', **row))
        uni = rs['miRNA_FFL'] | rs['TF_FFL']
        summary.append(dict(graph=tag, nodes=6, FFL_type='Unique_FFL', TF_TF_rule='either direction (B5)',
                            **{c: 'NA' for c in COLS}, No_of_FFL_instances=len(uni),
                            distinct_node_sets=len({frozenset(i) for i in uni})))
        for r in summary[-4:]:
            print(f"  6-node relaxed {r['FFL_type']:14s} FFL={r['No_of_FFL_instances']:,}  sets={r['distinct_node_sets']:,}")
        with open(os.path.join(HERE, 'sif', f'{tag}_composite_3to6_merged.sif'), 'w') as f:
            for (s, lab, t) in sorted(merged_comp): f.write(f'{s}\t{lab}\t{t}\n')
        print(f'  merged 3-to-6 composite network: {len(merged_comp)} SIF edges, '
              f'{len({x for (s, _, t) in merged_comp for x in (s, t)})} nodes')
    hdr = ['graph', 'nodes', 'FFL_type', 'TF_TF_rule'] + COLS + ['No_of_FFL_instances', 'distinct_node_sets']
    with open(os.path.join(HERE, 'bhat_table3_counts.csv'), 'w') as f:
        f.write(','.join(hdr) + '\n')
        for r in summary: f.write(','.join(str(r[h]) for h in hdr) + '\n')


if __name__ == '__main__':
    main()
