#!/usr/bin/env python3
"""
S0(b) overlaps and S0(c) MODEL6 in network terms (analyses/six_node_pattern).  Read-only on the project.

Graph: the nolegacy census graph (S1/graph_nolegacy_null.txt + S1/graph_nolegacy_labels.txt, arc order
and labels as documented in S1/make_s1_inputs.py).  Contracted graph for D1-D4 = the uncontracted one
minus every miRNA->TF arc whose TF->miRNA reverse exists (the census builder's rule).

(c) MODEL6-architecture = the COMP_C2_toggle six-node module of scripts/11_dynamics/dyn_models.py as a
    sign-agnostic, non-induced pattern (definition in S1/s1_null6.c; FULL = every interaction in the n6
    equations, DEF = the layer-defining edges only).  Counted here independently of the C counter
    (S1/s1_null6, OBS mode) and compared with it.  The signed edge list of every six-node module comes
    from dyn_models.py and is checked against the Jacobian signs stored in
    analyses/dynamics_checks/c2_interaction_signs.csv (COMP_C2_toggle and I1_miRNA_FFL);
    COMP_I1 differs from COMP_C2_toggle only in the sign of TF1 -> miRNA (s_TM = +1).
(b) overlaps between EA6 circuits, BHAT6 sets, CENSUS6 (c)-with-core modules and MODEL6 sets.
"""
import os, csv, gzip, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
S1 = os.path.join(HERE, '..', 'S1')
names = open(f'{S1}/node_names.txt').read().split('\n')[:587]
L = open(f'{S1}/graph_nolegacy_null.txt').read().rstrip('\n').split('\n')
nv, ne, _ = map(int, L[0].split()); ty = list(map(int, L[1].split()))
labs = list(map(int, open(f'{S1}/graph_nolegacy_labels.txt').read().split()))
arcs = [tuple(map(int, l.split()[:2])) for l in L[2:2 + ne]]
MIR, TF, GEN = 0, 1, 2
idx = {n: i for i, n in enumerate(names)}
out = []
P = lambda *a: out.append(' '.join(str(x) for x in a))


def build(nostring):
    A = [(u, v, l) for (u, v), l in zip(arcs, labs) if not (nostring and l == 5)]
    S = {(u, v) for u, v, _ in A}
    CO = collections.defaultdict(set)
    for u, v, _ in A:
        if ty[u] == MIR and ty[v] == TF and (v, u) in S: continue
        CO[u].add(v)
    by = collections.defaultdict(lambda: collections.defaultdict(set))
    for u, v, l in A:
        if l == 0: by['TM'][u].add(v)
        elif l == 1: by['MX'][u].add(v)
        elif l in (2, 3): by['TX'][u].add(v)
        elif l in (4, 5): by['GG'][u].add(v); by['GG'][v].add(u)
        elif l == 6: by['MM'][u].add(v); by['MM'][v].add(u)
    return S, CO, by


def model6(by):
    TM, MX, TX, GG, MM = by['TM'], by['MX'], by['TX'], by['GG'], by['MM']
    full, dfn = [], []
    for t1 in range(nv):
        if ty[t1] != TF: continue
        for m1 in TM[t1]:
            if t1 not in MX[m1]: continue                                   # TF1 <-> miR1
            for g1 in MX[m1] & TX[t1]:
                if g1 == t1 or ty[g1] == MIR: continue
                for g2 in GG[g1]:
                    if g2 in (t1, g1) or ty[g2] == MIR: continue
                    for m2 in MM[m1]:
                        if m2 == m1: continue
                        for t2 in TX[t1]:
                            if ty[t2] != TF or t2 in (t1, g1, g2): continue
                            if not (g1 in TX[t2] and m1 in TM[t2]): continue
                            inst = (t1, t2, m1, m2, g1, g2); dfn.append(inst)
                            if (g2 in TX[t1] and g2 in MX[m1] and m2 in TM[t1] and {g1, g2, t1} <= MX[m2]
                                    and g2 in TX[t2] and m2 in TM[t2]):
                                full.append(inst)
    return full, dfn


# ---------------- D1-D4 (port of acyclic() + check_module() of v2_census2.c) ----------------
def is_module(Sset, CO):
    V = list(Sset); n = len(V)
    adj = [[(i != j) and (V[j] in CO[V[i]]) for j in range(n)] for i in range(n)]
    indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
    gone = [False] * n; removed = 0; ch = True; d = indeg[:]
    while ch:
        ch = False
        for i in range(n):
            if not gone[i] and d[i] == 0:
                gone[i] = True; removed += 1; ch = True
                for j in range(n):
                    if adj[i][j]: d[j] -= 1
    if removed != n: return False
    outdeg = [sum(adj[i]) for i in range(n)]
    src = [i for i in range(n) if indeg[i] == 0]; snk = [i for i in range(n) if outdeg[i] == 0]
    if len(src) != 1 or len(snk) != 1 or src[0] == snk[0]: return False
    s, t = src[0], snk[0]
    def reach(start, fwd):
        seen = {start}; st = [start]
        while st:
            u = st.pop()
            for v in range(n):
                if (adj[u][v] if fwd else adj[v][u]) and v not in seen: seen.add(v); st.append(v)
        return seen
    if len(reach(s, True)) != n or len(reach(t, False)) != n: return False
    # node-split max flow >= 2
    NN = 2 * n; cap = [[0] * NN for _ in range(NN)]
    for i in range(n): cap[2 * i][2 * i + 1] = n if i in (s, t) else 1
    for i in range(n):
        for j in range(n):
            if adj[i][j]: cap[2 * i + 1][2 * j] = 1
    S_, T_ = 2 * s + 1, 2 * t; flow = 0
    for _ in range(2):
        prev = [-1] * NN; prev[S_] = S_; q = [S_]
        while q:
            u = q.pop(0)
            for v in range(NN):
                if prev[v] < 0 and cap[u][v] > 0: prev[v] = u; q.append(v)
        if prev[T_] < 0: break
        v = T_
        while v != S_: u = prev[v]; cap[u][v] -= 1; cap[v][u] += 1; v = u
        flow += 1
    return flow >= 2


rows_c = []
for tag, nost in (('with STRING', False), ('without STRING', True)):
    S, CO, by = build(nost)
    full, dfn = model6(by)
    fs = {frozenset(i) for i in full}; ds = {frozenset(i) for i in dfn}
    # compare with the C counter's listings (OBS mode)
    pre = 'nolegacy' if not nost else 'nolegacy_nostring'
    cf = [tuple(map(int, l.split())) for l in open(f'{S1}/obs/{pre}_model6_full_instances.tsv') if l.strip()]
    cd = [tuple(map(int, l.split())) for l in open(f'{S1}/obs/{pre}_model6_def_instances.tsv') if l.strip()]
    assert sorted(cf) == sorted(full) and sorted(cd) == sorted(dfn), 'C and Python MODEL6 counts differ'
    P(f'MODEL6 {tag}: FULL {len(full)} role assignments, {len(fs)} node sets; DEF {len(dfn)} role assignments, '
      f'{len(ds)} node sets  (identical to the C counter)')
    rows_c.append(dict(graph=f'nolegacy {tag}', model6_full_role_assignments=len(full), model6_full_node_sets=len(fs),
                       model6_def_role_assignments=len(dfn), model6_def_node_sets=len(ds)))
    if not nost:
        FULL, DEF, FS, DS, CO0, BY0 = full, dfn, fs, ds, CO, by
        with open(f'{HERE}/model6_full_instances_nolegacy.tsv', 'w') as f:
            f.write('TF1\tTF2\tmiR1\tmiR2\tG1\tG2\n')
            for i in full: f.write('\t'.join(names[x] for x in i) + '\n')
        with gzip.open(f'{HERE}/model6_def_instances_nolegacy.tsv.gz', 'wt') as f:
            f.write('TF1\tTF2\tmiR1\tmiR2\tG1\tG2\n')
            for i in dfn: f.write('\t'.join(names[x] for x in i) + '\n')
with open(f'{HERE}/s0c_model6_counts.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_c[0])); w.writeheader(); w.writerows(rows_c)
# node / pair usage of the FULL instances
P('MODEL6 FULL (with STRING) nodes: TF1', sorted({names[i[0]] for i in FULL}), '| TF2', sorted({names[i[1]] for i in FULL}),
  '| miR1', sorted({names[i[2]] for i in FULL}), '| miR2', sorted({names[i[3]] for i in FULL}),
  '| G1', sorted({names[i[4]] for i in FULL}), '| G2', sorted({names[i[5]] for i in FULL}))

# ---------------- signed edge lists (S0c) ----------------
J = list(csv.DictReader(open(f'{REV}/analyses/dynamics_checks/c2_interaction_signs.csv')))
jac = {(r['family'], r['source'], r['target']): r['sign'] for r in J if r['module'] == 'n6'}
node = lambda s: s.split('_')[0]
E = [  # (source, target, sign C2 / I1-composite / I1 miRNA-FFL, mechanism, layer, network arc, FULL, DEF)
 ('TF1', 'miR1', '-', '+', '+', 'TF1 promoter input u_m (Hill); gated with TF2 (AND = product)', 'core', 'TF_miRNA', 1, 1),
 ('miR1', 'G1', '-', '-', '-', 'RISC: accelerated G1 mRNA decay (lam_y), titration (theta), translational block of G1 protein (AND gate)', 'core', 'miRNA_target', 1, 1),
 ('TF1', 'G1', '+', '+', '+', 'TF1 transcription of G1 (u_TY); gated with TF2', 'core', 'TF_target', 1, 1),
 ('miR1', 'TF1', '-', '-', 'absent', 'RISC on the TF1 transcript: translational block + accelerated decay (lam_x) [composite only]', 'core', 'miRNA_target (reciprocal)', 1, 1),
 ('G1', 'G2', '+', '+', '+', 'mutual protein stabilisation, saturable kappa*f+(p2) (and G2 -> G1)', 'gene-gene', 'gene_gene (STRING/deposit)', 1, 1),
 ('G2', 'G1', '+', '+', '+', 'mutual protein stabilisation (reciprocal of the above)', 'gene-gene', 'gene_gene (same link)', 1, 1),
 ('TF1', 'G2', '+', '+', '+', 'G2 is co-regulated: transcription rho_g2*u_TY (same TF input)', 'gene-gene', 'TF_target', 1, 0),
 ('miR1', 'G2', '-', '-', '-', 'RISC on G2 mRNA (decay, titration) and G2 translation', 'gene-gene', 'miRNA_target', 1, 0),
 ('miR1', 'miR2', 'coupled', 'coupled', 'coupled', 'co-transcription from the shared cluster promoter (miR2 production = rho*u_m) and shared saturable RISC pool Meff; no direct term', 'miRNA-miRNA', 'miRNA_miRNA (10 kb)', 1, 1),
 ('TF1', 'miR2', '-', '+', '+', 'shared promoter: miR2 production = rho*u_m', 'miRNA-miRNA', 'TF_miRNA', 1, 0),
 ('miR2', 'G1', '-', '-', '-', 'shared RISC pool: Meff = f(m1 + phi*m2) acts on G1', 'miRNA-miRNA', 'miRNA_target', 1, 0),
 ('miR2', 'G2', '-', '-', '-', 'shared RISC pool acts on G2', 'miRNA-miRNA', 'miRNA_target', 1, 0),
 ('miR2', 'TF1', '-', '-', 'absent', 'shared RISC pool acts on the TF1 transcript [composite only]', 'miRNA-miRNA', 'miRNA_target', 1, 0),
 ('TF1', 'TF2', '+', '+', '+', 'TF1 transcription of TF2 (s_T12, Hill)', 'TF-TF', 'TF_target (TF->TF)', 1, 1),
 ('TF2', 'G1', '+', '+', '+', 'TF2 input to G1 transcription, AND-gated with TF1 (s_T2Y)', 'TF-TF', 'TF_target', 1, 1),
 ('TF2', 'miR1', '+', '+', '+', 'TF2 input to the miRNA promoter, AND-gated with TF1 (s_T2M)', 'TF-TF', 'TF_miRNA', 1, 1),
 ('TF2', 'G2', '+', '+', '+', 'G2 transcription uses the same gated TF input as G1', 'TF-TF', 'TF_target', 1, 0),
 ('TF2', 'miR2', '+', '+', '+', 'miR2 production uses the same gated promoter input as miR1', 'TF-TF', 'TF_miRNA', 1, 0),
 ('G1', 'miR1', '-', '-', '-', 'titration: G1 mRNA sequesters miRNA (theta); reverse action of the miR1 -> G1 arc', 'core', '(same arc as miR1 -> G1)', 0, 0),
 ('G1', 'miR2', '-', '-', '-', 'titration by G1 mRNA', 'miRNA-miRNA', '(same arc as miR2 -> G1)', 0, 0),
 ('G2', 'miR1', '-', '-', '-', 'titration by G2 mRNA (weight w2)', 'gene-gene', '(same arc as miR1 -> G2)', 0, 0),
 ('G2', 'miR2', '-', '-', '-', 'titration by G2 mRNA', 'gene-gene', '(same arc as miR2 -> G2)', 0, 0),
]
chk = []
for s, t, c2, ci1, i1, mech, layer, arc, f_, d_ in E:
    for fam, sg in (('COMP_C2_toggle', c2), ('I1_miRNA_FFL', i1)):
        if sg in ('absent', 'coupled'): continue
        tgt_states = {'G1': ['G1_mRNA', 'G1_protein'], 'G2': ['G2_mRNA', 'G2_protein']}.get(t, [t])
        src_states = {'G1': ['G1_mRNA', 'G1_protein'], 'G2': ['G2_mRNA', 'G2_protein']}.get(s, [s])
        found = {jac[(fam, a, b)] for a in src_states for b in tgt_states if (fam, a, b) in jac}
        chk.append((fam, s, t, sg, found))
bad = [c for c in chk if c[4] != {c[3]}]
P('signed edge list checked against the stored n6 Jacobian signs:', len(chk) - len(bad), 'of', len(chk), 'agree;',
  'disagreements:', bad)
with open(f'{HERE}/s0c_model6_signed_edge_list.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['source', 'target', 'sign_COMP_C2_toggle_n6', 'sign_COMP_I1_negfeedback_n6', 'sign_I1_miRNA_FFL_n6',
                'mechanism_in_dyn_models', 'layer', 'network_arc_class', 'required_MODEL6_FULL', 'required_MODEL6_DEF'])
    for r in E: w.writerow(r)

# ---------------- S0(b) overlaps ----------------
def bhat_sets(klass):
    p = f'{REV}/analyses/six_node_pattern_networks/instances/nolegacy_6node_{klass}.tsv.gz'
    with gzip.open(p, 'rt') as f:
        next(f); return {frozenset(l.split()) for l in f if l.strip()}
B = {k: bhat_sets(k) for k in ('Composite_FFL', 'miRNA_FFL', 'TF_FFL')}
P('BHAT6 node sets (nolegacy): composite', len(B['Composite_FFL']), 'miRNA-FFL', len(B['miRNA_FFL']), 'TF-FFL', len(B['TF_FFL']))
ea6 = {frozenset(r['members'].split(';')) for r in csv.DictReader(open(f'{REV}/results/circuit_level_3node_vs_6node.csv'))
       if r['circuit_class'] == '6-node'}
anyB = B['Composite_FFL'] | B['miRNA_FFL'] | B['TF_FFL']
FSn = {frozenset(names[x] for x in s) for s in FS}; DSn = {frozenset(names[x] for x in s) for s in DS}
P('EA6 circuits:', len(ea6), '| that are BHAT6 composite sets:', len(ea6 & B['Composite_FFL']),
  '| BHAT6 sets of any class:', len(ea6 & anyB), '| MODEL6 FULL sets:', len(ea6 & FSn), '| MODEL6 DEF sets:', len(ea6 & DSn))
for k in B:
    n_pass = sum(is_module({idx[x] for x in s}, CO0) for s in B[k])
    P(f'BHAT6 {k} node sets passing D1-D4 in the contracted nolegacy census graph: {n_pass} of {len(B[k])}')
c6 = [frozenset(map(int, l.split())) for l in open(f'{S1}/obs/nolegacy_census6c_modules.tsv') if l.strip()]
P('CENSUS6 (c)-with-core modules:', len(c6), '| equal to a MODEL6 FULL node set:', len(set(c6) & FS),
  '| equal to a MODEL6 DEF node set:', len(set(c6) & DS), '| equal to a BHAT6 composite set:',
  len({frozenset(names[x] for x in s) for s in c6} & B['Composite_FFL']))
nf = sum(is_module(s, CO0) for s in FS); nd = sum(is_module(s, CO0) for s in DS)
P('MODEL6 node sets passing D1-D4 themselves: FULL', nf, 'of', len(FS), '; DEF', nd, 'of', len(DS))
P('MODEL6 FULL sets that are BHAT6 composite sets:', len(FSn & B['Composite_FFL']), '; BHAT6 any class:', len(FSn & anyB))
open(f'{HERE}/s0bc_overlaps_and_model6.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
