#!/usr/bin/env python3
"""One network file set per feed-forward-loop (FFL) class and size: 3node_, 4node_, 5node_ and 6node_ x miRNA_FFL,
TF_FFL and composite_FFL (12 networks), with README.md and MANIFEST.tsv.

Every network follows the construction of Bhat et al. (2024, Curr Bioinform 19:73-90) on the census graph: the analysed
network (the interaction table without the 30 miRNA-miRNA edges of the original circuits) plus the census layers
(TRRUST TF-target arcs, STRING gene-gene arcs, miRNA pairs within 10 kb).  Patterns are the compulsory edges of their
Fig. 1; roles are distinct nodes, M = miRNA, T = TF, G = gene:
  3-node  M1-T1 link, M1->G1, T1->G1
  4-node  M1-T1 link, M1->G1, T1->G2, G1-G2
  5-node  M1-T1 link, M1->G1, M2->G1, M1-M2, T1->G2, G1-G2
  6-node  M1-T1 link, M1->G1, M2->G1, M1-M2, T1-T2 link, T2->G2, G1-G2
  miRNA-FFL: M1->T1 (6-node also T1->T2); TF-FFL: T1->M1 (6-node also T2->T1); composite: both directions.
G1-G2 (gene-gene) and M1-M2 (miRNA pair) are undirected links.  An instance is one role assignment; the classes overlap
(every composite instance is also a miRNA-FFL and a TF-FFL instance).  This is the enumeration of
analyses/six_node_pattern_networks/build_bhat_networks.py (graph 'nolegacy', TF-TF rule of Fig. 1d).

The script stops unless (a) the six-node composite count is the 3,898 of Note S2 and (b) the 3-node instances are
exactly the typed cores of results/v5/tables/TableS2_ffl_cores.csv whose target is a gene, with the overlapping
classes: composite = the composite cores; TF-FFL = composite + TF-FFL cores; miRNA-FFL = composite + miRNA-FFL cores.

Files per network <name>:
  <name>.sif        source <tab> interacts with <tab> target, one line per edge (the original submission's format);
                    undirected links appear once, the two names sorted
  <name>_edges.tsv  source, target, edge_type, sign, sign_source, evidence_tier, layer, n_ffl
  <name>_nodes.tsv  node, type, n_ffl
n_ffl = the number of instances of that network that contain the edge or node.
sign and sign_source follow the published interaction table: miRNA -> target -1 (mechanistic); TF edges +1 or -1 from
the TRRUST or TransmiR mode (Activation, Repression); otherwise blank, 'unannotated (excluded from sign-dependent
analyses)'.  evidence_tier is the table's tier (miRNA -> target edges).  Census-layer links carry their layer's source.

usage: python3 make_ffl_networks.py <analysis root> <folder of the original SIF files> <output folder>"""
import os, sys, csv, hashlib, collections
sys.dont_write_bytecode = True
UNANN = 'unannotated (excluded from sign-dependent analyses)'
CLASSES = ('miRNA_FFL', 'TF_FFL', 'composite_FFL')
ROLES = {3: ('M1', 'T1', 'G1'), 4: ('M1', 'T1', 'G1', 'G2'), 5: ('M1', 'M2', 'T1', 'G1', 'G2'), 6: ('M1', 'M2', 'T1', 'T2', 'G1', 'G2')}
NOTE_S2 = 3898                                        # six-node composite instances reported in Note S2
LAYER = {'canonical': 'analysed network', 'layer_TF_target': 'census layer: TRRUST TF-target',
         'layer_gene_gene': 'census layer: STRING gene-gene', 'layer_miRNA_miRNA': 'census layer: miRNA pair within 10 kb'}


def tsv(path): return list(csv.DictReader(open(path, newline=''), delimiter='\t'))


def load(root):
    nodes = {r['name']: r['type'] for r in tsv(os.path.join(root, 'data/canonical_nodes.tsv'))}
    analysed = {}                                    # the interaction table without the 30 legacy miRNA-miRNA edges
    for r in tsv(os.path.join(root, 'data/canonical_edges.tsv')):
        if r['edge_type'] != 'miRNA_miRNA': analysed[(r['source'], r['target'])] = r['edge_type']
    census = {e: (t, 'canonical', None) for e, t in analysed.items()}

    def add(fn, sc, tc, et, lab, undirected=False, keep=lambda r: True):
        for r in tsv(os.path.join(root, 'data', fn)):
            if not keep(r): continue
            a, b = r[sc], r[tc]
            if a not in nodes or b not in nodes or a == b: continue
            for e in ([(a, b), (b, a)] if undirected else [(a, b)]):
                if e not in census: census[e] = (et, lab, r)
    # the census graph of scripts/03_ffl_census/03_ffl_census_nolegacymirna.py build_graph(), before its TF<->miRNA contraction
    add('layer_TF_target.tsv', 'source', 'target', 'TF_target', 'layer_TF_target')
    add('layer_gene_gene.tsv', 'source', 'target', 'gene_gene', 'layer_gene_gene')
    add('layer_miRNA_miRNA.tsv', 'miRNA_1', 'miRNA_2', 'miRNA_miRNA', 'layer_miRNA_miRNA', undirected=True,
        keep=lambda r: str(r.get('threshold_10kb', '')).upper() in ('TRUE', '1'))
    table = {(r['source'], r['target']): r for r in csv.DictReader(open(os.path.join(root, 'results/v5/tables/TableS1_all_interactions.csv'), newline=''))}
    assert len(nodes) == 587 and len(analysed) == 6829 and len(census) == 9226, (len(nodes), len(analysed), len(census))
    assert all(e in table and table[e]['edge_type'] == t for e, t in analysed.items())
    return nodes, analysed, census, table


def sign_of(edge_type, mode):
    if edge_type == 'miRNA_target': return '-1', 'mechanistic (miRNA repression)'
    src = {'TF_target': 'TRRUST', 'TF_miRNA': 'TransmiR'}.get(edge_type)
    if src and mode == 'Activation': return '1', src
    if src and mode == 'Repression': return '-1', src
    return '', UNANN


def attributes(e, census, table):
    """edge_type, sign, sign_source, evidence_tier, layer of one census arc"""
    et, lab, r = census[e]
    if lab == 'canonical':
        t = table[e]; mode = t['trrust_mode'] if et == 'TF_target' else t['transmir_mode'] if et == 'TF_miRNA' else ''
        return (et,) + sign_of(et, mode) + (t['tier'], LAYER[lab])
    if lab == 'layer_TF_target': return (et,) + sign_of(et, r['mode']) + ('', LAYER[lab])
    if lab == 'layer_gene_gene': return (et, '', 'unsigned (STRING association)', '', LAYER[lab])
    return (et, '', 'unsigned (co-transcription within 10 kb)', '', LAYER[lab])


def check_published_signs(analysed, table):
    """the published interaction table carries the same sign and sign_source; older copies carry a placeholder sign"""
    if not any(r['sign_source'] == UNANN for r in table.values()): return 'not checked (this copy of the table predates the published sign columns)'
    bad = [e for e, et in analysed.items() if attributes(e, {e: (et, 'canonical', None)}, table)[1:3] != (table[e]['sign'], table[e]['sign_source'])]
    assert not bad, ('sign differs from the published table', bad[:5])
    return f'checked: equal to the table for all {len(analysed):,} edges'


class Bhat:
    """the enumeration of analyses/six_node_pattern_networks build_bhat_networks.py, on the census graph"""
    def __init__(self, nodes, census):
        def cls(a, b):
            t = census[(a, b)][0]
            if t == 'TF_target': return 'TF->TF' if nodes[b] == 'TF' else 'TF->Gene'
            if t == 'miRNA_target': return 'miRNA->TF' if nodes[b] == 'TF' else 'miRNA->Gene'
            return {'TF_miRNA': 'TF->miRNA', 'gene_gene': 'Gene->Gene', 'miRNA_miRNA': 'miRNA-miRNA'}[t]
        ty = nodes
        self.MG = collections.defaultdict(set); self.TG = collections.defaultdict(set)
        self.MT, self.TM, self.TT = set(), set(), set()
        self.GG = collections.defaultdict(set); self.MM = collections.defaultdict(set)
        for (u, v) in census:
            c = cls(u, v)
            if c == 'miRNA->Gene' and ty[v] == 'Gene': self.MG[u].add(v)
            elif c == 'TF->Gene' and ty[u] == 'TF' and ty[v] == 'Gene': self.TG[u].add(v)
            elif c == 'miRNA->TF': self.MT.add((u, v))
            elif c == 'TF->miRNA': self.TM.add((u, v))
            elif c == 'TF->TF' and ty[u] == 'TF' and ty[v] == 'TF' and u != v: self.TT.add((u, v))
            elif c == 'Gene->Gene' and ty[u] == 'Gene' and ty[v] == 'Gene' and u != v: self.GG[u].add(v); self.GG[v].add(u)
            elif c == 'miRNA-miRNA' and u != v: self.MM[u].add(v); self.MM[v].add(u)
        self.Tout = collections.defaultdict(set); self.Tin = collections.defaultdict(set)
        for (a, b) in self.TT: self.Tout[a].add(b); self.Tin[b].add(a)
        self.links = {'miRNA_FFL': sorted(self.MT), 'TF_FFL': sorted((m, t) for (t, m) in self.TM),
                      'composite_FFL': sorted(p for p in self.MT if (p[1], p[0]) in self.TM)}

    def tf2(self, klass, t1):
        return self.Tout[t1] if klass == 'miRNA_FFL' else self.Tin[t1] if klass == 'TF_FFL' else self.Tout[t1] & self.Tin[t1]

    def enumerate(self, n, klass):
        out = []
        for (m1, t1) in self.links[klass]:
            if n == 3:
                out += [(m1, t1, g1) for g1 in sorted(self.MG[m1] & self.TG[t1])]; continue
            for g1 in sorted(self.MG[m1]):
                if n == 4:
                    out += [(m1, t1, g1, g2) for g2 in sorted(self.GG[g1] & self.TG[t1])]; continue
                m2s = sorted(m2 for m2 in self.MM[m1] if g1 in self.MG[m2])
                if not m2s: continue
                if n == 5:
                    out += [(m1, m2, t1, g1, g2) for m2 in m2s for g2 in sorted(self.GG[g1] & self.TG[t1])]; continue
                for t2 in sorted(self.tf2(klass, t1)):
                    out += [(m1, m2, t1, t2, g1, g2) for m2 in m2s for g2 in sorted(self.GG[g1] & self.TG[t2])]
        return out

    @staticmethod
    def edges_of(n, klass, inst):
        """compulsory edges of one instance: directed arcs (s, t, False) and undirected links (a, b, True), a < b"""
        r = dict(zip(ROLES[n], inst)); E = []
        und = lambda a, b: E.append(tuple(sorted((a, b))) + (True,))
        if klass in ('miRNA_FFL', 'composite_FFL'): E.append((r['M1'], r['T1'], False))
        if klass in ('TF_FFL', 'composite_FFL'): E.append((r['T1'], r['M1'], False))
        E.append((r['M1'], r['G1'], False))
        if n == 3: E.append((r['T1'], r['G1'], False)); return E
        und(r['G1'], r['G2'])
        if n >= 5: E.append((r['M2'], r['G1'], False)); und(r['M1'], r['M2'])
        if n <= 5: E.append((r['T1'], r['G2'], False)); return E
        E.append((r['T2'], r['G2'], False))
        if klass in ('miRNA_FFL', 'composite_FFL'): E.append((r['T1'], r['T2'], False))
        if klass in ('TF_FFL', 'composite_FFL'): E.append((r['T2'], r['T1'], False))
        return E


def check_typed_cores(root, nodes, three):
    """3-node instances = the typed cores of the paper whose target is a gene, with overlapping classes"""
    s2 = [r for r in csv.DictReader(open(os.path.join(root, 'results/v5/tables/TableS2_ffl_cores.csv'), newline=''))]
    assert len(s2) == 1649
    g = collections.defaultdict(set)                  # as (miRNA, TF, gene)
    for r in s2:
        if nodes[r['target']] != 'Gene': continue
        mir, tf = (r['intermediate'], r['regulator']) if r['type_R'] == 'TF' else (r['regulator'], r['intermediate'])
        g[r['class']].add((mir, tf, r['target']))
    want = {'composite_FFL': g['Composite-FFL'], 'TF_FFL': g['Composite-FFL'] | g['TF-FFL'],
            'miRNA_FFL': g['Composite-FFL'] | g['miRNA-FFL']}
    for k in CLASSES: assert set(three[k]) == want[k], k
    n_tf_target = sum(nodes[r['target']] == 'TF' for r in s2)
    return {k: len(v) for k, v in g.items()}, n_tf_target


def original_counts(folder):
    rows = []
    for f in ('3-miR', '3-TF', '3-Comp', '4-TF', '5-TF', '6-TF'):
        L = [l.rstrip('\r\n').split('\t') for l in open(os.path.join(folder, f + '.sif'), newline='') if l.strip()]
        names = {x for p in L for x in (p[0], p[2])}
        rows.append(dict(file=f + '.sif', lines=len(L), distinct_edges=len({(p[0], p[2]) for p in L}), nodes=len(names),
                         nodes_case_merged=len({x.lower() for x in names})))
    return rows


def main(root, orig_folder, out):
    os.makedirs(out, exist_ok=True)
    nodes, analysed, census, table = load(root)
    sign_check = check_published_signs(analysed, table)
    bh = Bhat(nodes, census)
    nets = collections.OrderedDict()                   # name -> (edges {key: n}, nodes {node: n}, instances, node sets)
    three = {}
    for size in (3, 4, 5, 6):
        for klass in CLASSES:
            insts = bh.enumerate(size, klass); assert len(set(insts)) == len(insts)
            if size == 3: three[klass] = insts
            E, V = collections.Counter(), collections.Counter()
            for inst in insts:
                V.update(set(inst)); E.update(set(bh.edges_of(size, klass, inst)))
            nets[f'{size}node_{klass}'] = (E, V, len(insts), len({frozenset(i) for i in insts}))
    assert nets['6node_composite_FFL'][2] == NOTE_S2, nets['6node_composite_FFL'][2]
    core_counts, n_tf_target = check_typed_cores(root, nodes, three)
    summary = []
    for name, (E, V, n, nsets) in nets.items():
        rows = []
        for (a, b, undirected), k in E.items():
            if undirected:
                arcs = sorted((e for e in ((a, b), (b, a)) if e in census), key=lambda e: census[e][1] != 'canonical')
                at = attributes(arcs[0], census, table)
            else: at = attributes((a, b), census, table)
            rows.append((a, b) + at + (k,))
        rows.sort(key=lambda r: (r[0], r[1]))
        with open(os.path.join(out, name + '.sif'), 'w', newline='') as f:
            for r in rows: f.write(f'{r[0]}\tinteracts with\t{r[1]}\n')
        with open(os.path.join(out, name + '_edges.tsv'), 'w', newline='') as f:
            f.write('source\ttarget\tedge_type\tsign\tsign_source\tevidence_tier\tlayer\tn_ffl\n')
            for r in rows: f.write('\t'.join(map(str, r)) + '\n')
        with open(os.path.join(out, name + '_nodes.tsv'), 'w', newline='') as f:
            f.write('node\ttype\tn_ffl\n')
            for v in sorted(V): f.write(f'{v}\t{nodes[v]}\t{V[v]}\n')
        et = collections.Counter(r[2] for r in rows); ly = collections.Counter(r[6] for r in rows)
        vt = collections.Counter(nodes[v] for v in V)
        summary.append(dict(network=name, n_ffl=n, node_sets=nsets, nodes=len(V), TF=vt['TF'], miRNA=vt['miRNA'], Gene=vt['Gene'],
                            edges=len(rows), by_type=dict(sorted(et.items())), by_layer=dict(sorted(ly.items()))))
    write_readme(out, summary, original_counts(orig_folder), core_counts, n_tf_target)
    me = os.path.abspath(__file__)                      # the output folder carries the script that wrote it
    if os.path.dirname(me) != os.path.abspath(out): open(os.path.join(out, os.path.basename(me)), 'wb').write(open(me, 'rb').read())
    write_manifest(out, summary)
    for s in summary: print(f"{s['network']:22s} FFLs {s['n_ffl']:6d}  node sets {s['node_sets']:6d}  nodes {s['nodes']:4d}  edges {s['edges']:5d}")
    print('signs:', sign_check)


def write_readme(out, S, orig, core_counts, n_tf_target):
    g = {s['network']: s for s in S}
    W = []; w = W.append
    w('# Feed-forward-loop networks, one file set per class and size\n')
    w('Twelve networks: three FFL classes (miRNA-FFL, TF-FFL, composite FFL) at three to six nodes. They are the revised')
    w('counterparts of the six SIF files of the original submission (`analyses/original_submission_code/SIF_files/`).')
    w('`make_ffl_networks.py` writes every file here, and `MANIFEST.tsv` lists them with their md5.\n')
    w('## Files\n')
    w('For each network `<name>` (for example `6node_composite_FFL`):\n')
    w('- `<name>.sif`: one line per edge, `source <tab> interacts with <tab> target`, the format of the original files.')
    w('  Undirected links (gene-gene, miRNA pair) appear once, with the two names sorted.')
    w('- `<name>_edges.tsv`: `source`, `target`, `edge_type`, `sign`, `sign_source`, `evidence_tier`, `layer`, `n_ffl`.')
    w('- `<name>_nodes.tsv`: `node`, `type` (TF, miRNA or Gene), `n_ffl`.\n')
    w('`n_ffl` is the number of FFL instances of that network that contain the edge or node (defined below).\n')
    w('## Graph\n')
    w('All twelve networks are built on the census graph (the graph of Table 3), 9,226 arcs:')
    w('- the analysed network: the interaction table (Table S2) without the 30 miRNA-miRNA edges of the original circuits,')
    w('  587 nodes and 6,829 edges;')
    w('- 604 TF-to-target arcs from TRRUST;')
    w('- 833 gene-gene arcs from STRING;')
    w('- 960 arcs joining 480 pairs of miRNAs that lie within 10 kb of each other.\n')
    w('The `layer` column of every edge file names the source of each edge. The 3-node networks use only edges of the')
    w('analysed network.\n')
    w('## Definition: Bhat et al. (2024), at every size\n')
    w('Each size follows a fixed pattern of compulsory edges. M1 and M2 are miRNAs, T1 and T2 are TFs, G1 and G2 are genes,')
    w('and all are distinct nodes:\n')
    w('| size | compulsory edges |'); w('|---|---|')
    w('| 3 | M1-T1 link, M1 -> G1, T1 -> G1 |')
    w('| 4 | M1-T1 link, M1 -> G1, T1 -> G2, G1 - G2 |')
    w('| 5 | M1-T1 link, M1 -> G1, M2 -> G1, M1 - M2, T1 -> G2, G1 - G2 |')
    w('| 6 | M1-T1 link, M1 -> G1, M2 -> G1, M1 - M2, T1-T2 link, T2 -> G2, G1 - G2 |\n')
    w('The class sets the direction of the links:')
    w('- **miRNA-FFL**: M1 -> T1, and at six nodes also T1 -> T2.')
    w('- **TF-FFL**: T1 -> M1, and at six nodes also T2 -> T1.')
    w('- **Composite FFL**: both directions.\n')
    w('G1 - G2 (gene-gene) and M1 - M2 (miRNA pair) are undirected. Extra edges among the nodes are allowed. An instance is')
    w('one assignment of nodes to the roles.\n')
    w('The classes overlap: every composite instance is also a miRNA-FFL and a TF-FFL instance. In this network almost every')
    w('TF -> miRNA edge is matched by a miRNA -> TF edge (1,223 of 1,239), so most TF-FFL and many miRNA-FFL instances are')
    w('composite. The six-node composite FFL is the pattern of Note S2.\n')
    w('## Counts\n')
    w('| network | FFL instances | distinct node sets | nodes (TF / miRNA / gene) | edges | edges by type |')
    w('|---|---|---|---|---|---|')
    for s in S:
        w(f"| {s['network']} | {s['n_ffl']:,} | {s['node_sets']:,} | {s['nodes']:,} ({s['TF']} / {s['miRNA']} / {s['Gene']}) | {s['edges']:,} | "
          + ', '.join(f'{k} {v:,}' for k, v in s['by_type'].items()) + ' |')
    w('\nEdges by source (the `layer` column):\n')
    w('| network | analysed network | TRRUST layer | STRING layer | 10 kb miRNA-pair layer |'); w('|---|---|---|---|---|')
    for s in S:
        L = s['by_layer']; w(f"| {s['network']} | " + ' | '.join(f"{L.get(LAYER[k], 0):,}" for k in ('canonical', 'layer_TF_target', 'layer_gene_gene', 'layer_miRNA_miRNA')) + ' |')
    comp, tfo, mio = core_counts.get('Composite-FFL', 0), core_counts.get('TF-FFL', 0), core_counts.get('miRNA-FFL', 0)
    w('\n## How these counts relate to the paper\n')
    w(f"- **Note S2:** {g['6node_composite_FFL']['n_ffl']:,} six-node composite instances. The file reproduces this, and the")
    w('  script stops otherwise.')
    w('- **Table S4 and Table 3, three-node row:** the paper counts 1,649 typed three-node cores with exclusive classes:')
    w('  1,434 composite, 206 miRNA-FFL and 9 TF-FFL. There, a composite core counts only as composite.')
    w(f'  - {n_tf_target} of those cores have a TF as target; the Bhat pattern requires a gene target. The other '
      f'{comp + tfo + mio:,} have a gene target:')
    w(f'    {comp:,} composite, {mio:,} miRNA-FFL and {tfo:,} TF-FFL.')
    w(f"  - The 3-node files are exactly these cores, with overlapping classes:")
    w(f"    - composite: {g['3node_composite_FFL']['n_ffl']:,} = {comp:,};")
    w(f"    - TF-FFL: {g['3node_TF_FFL']['n_ffl']:,} = {comp:,} + {tfo:,};")
    w(f"    - miRNA-FFL: {g['3node_miRNA_FFL']['n_ffl']:,} = {comp:,} + {mio:,}.")
    w('  - The script checks this correspondence core by core, and stops otherwise.')
    w('- **Table 3, other counts:** 5,833 at n = 3 and the modules at four to seven nodes are D1-D4 modules counted')
    w('  irrespective of node type. They are a different object from these typed networks.\n')
    w('## The original submission\'s files\n')
    w('Counts are taken from the files as deposited. Many miRNAs occur twice with different capitalisation (`hsa-miR-` and')
    w('`hsa-mir-`), so nodes are given both as written and with capitalisation merged.\n')
    w('| file | lines | distinct edges | nodes as written | nodes, capitalisation merged | revised counterpart |')
    w('|---|---|---|---|---|---|')
    cp = {'3-miR.sif': '3node_miRNA_FFL', '3-TF.sif': '3node_TF_FFL', '3-Comp.sif': '3node_composite_FFL',
          '4-TF.sif': '4node_* (exemplar circuit)', '5-TF.sif': '5node_* (exemplar circuit)', '6-TF.sif': '6node_* (exemplar circuit)'}
    for o in orig:
        w(f"| {o['file']} | {o['lines']:,} | {o['distinct_edges']:,} | {o['nodes']:,} | {o['nodes_case_merged']:,} | {cp[o['file']]} |")
    w('\nThe original files are of a different kind from the revised ones:')
    w('- The original 3-node files are the full networks of the original analysis (5,360 to 6,526 distinct edges each).')
    w('- The original 4-, 5- and 6-node files are single exemplar circuits. The 5- and 6-node files hold 20 and 10')
    w('  miRNA-miRNA edges: the 30 edges the revision excludes.')
    w('- Each revised file contains only the edges that the FFL instances of its class use.\n')
    w('## Signs\n')
    w('`sign` and `sign_source` follow the published interaction table:')
    w('- miRNA -> target: -1 (mechanistic repression).')
    w('- TF edges: +1 or -1 from the TRRUST or TransmiR mode (Activation or Repression).')
    w(f'- All others: blank, with sign_source `{UNANN}`.\n')
    w('Census-layer arcs: TRRUST arcs follow the same rule; STRING and miRNA-pair links are unsigned. `evidence_tier` is')
    w('the table\'s tier for miRNA -> target edges. When the table read is the published one, the script also checks that')
    w('these values equal its `sign` and `sign_source` columns for every edge, and stops otherwise.\n')
    w('## Rebuild\n')
    w('`python3 make_ffl_networks.py <analysis root> <folder of the original SIF files> <output folder>`. The script reads:')
    w('- `data/canonical_nodes.tsv` and `data/canonical_edges.tsv`;')
    w('- the three layer files `data/layer_TF_target.tsv`, `data/layer_gene_gene.tsv` and `data/layer_miRNA_miRNA.tsv`;')
    w('- the interaction table `results/v5/tables/TableS1_all_interactions.csv`;')
    w('- the typed-core table `results/v5/tables/TableS2_ffl_cores.csv`, for the 3-node check.')
    open(os.path.join(out, 'README.md'), 'w').write('\n'.join(W) + '\n')


def write_manifest(out, S):
    backs = {s['network']: f"{s['network']}: {s['n_ffl']:,} instances, {s['nodes']:,} nodes, {s['edges']:,} edges" for s in S}
    rep = {'6node_composite_FFL': 'Note S2 (3,898 instances)',
           '3node_composite_FFL': 'Table S4 composite cores with a gene target',
           '3node_TF_FFL': 'Table S4 composite + TF-FFL cores with a gene target',
           '3node_miRNA_FFL': 'Table S4 composite + miRNA-FFL cores with a gene target'}
    md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
    rows = []
    for fn in sorted(os.listdir(out)):
        if fn == 'MANIFEST.tsv' or not os.path.isfile(os.path.join(out, fn)): continue
        net = fn.replace('_edges.tsv', '').replace('_nodes.tsv', '').replace('.sif', '')
        b = (backs[net] + ('; ' + rep[net] if net in rep else '')) if net in backs else \
            ('the definitions and counts of every network' if fn == 'README.md' else 'writes every file here')
        p = os.path.join(out, fn); rows.append((fn, md5(p), os.path.getsize(p), b))
    with open(os.path.join(out, 'MANIFEST.tsv'), 'w') as f:
        f.write('path\tmd5\tsize_bytes\tbacks\n')
        for r in rows: f.write('\t'.join(map(str, r)) + '\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
