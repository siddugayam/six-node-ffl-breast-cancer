#!/usr/bin/env python3
"""
S5 (a)-(c): which real six-node circuits carry the modelled loops?  (analyses/six_node_pattern)  Read-only.

Graph: nolegacy census graph with STRING (S1/graph_nolegacy_null.txt + labels; 9,226 arcs).
Signs (fixed before any count was made):
  miRNA_target (label 1)            -1, mechanistic repression
  TF_target, analysed (2)           TRRUST mode in Table S1 (published as Table S2): Activation +1,
                                    Repression -1, Unknown / Ambiguous unsigned ('assumed' +1 NOT used)
  TF_target, TRRUST layer (3)       mode column of data/layer_TF_target.tsv, same rule
  TF_miRNA (0)                      TransmiR mode in Table S1: Activation +1, Repression -1,
                                    Regulation / Ambiguous unsigned
  gene_gene deposit (4), STRING (5), 10 kb (6)   unsigned
Evidence per edge: miRNA_target -> Table S1 tier (strong / weak / predicted_only) and databases;
TF_target -> TRRUST with its PMID count; TF_miRNA -> TransmiR (literature) with its PMID count;
STRING -> combined score and network type; 10 kb -> distance (bp).  Tiers exist only for
miRNA-target edges, so circuits are sorted by their weakest miRNA-target tier (strong < weak <
predicted_only), then by the smallest TRRUST PMID count among their TF_target arcs (larger first).

Objects and the model roles (TF1, TF2, miR1, miR2, G1, G2) mapped onto them:
  BHAT6 composite   T1->TF1, T2->TF2, M1->miR1, M2->miR2, G1->G1, G2->G2 (3,898 role assignments)
  MODEL6 FULL/DEF   as enumerated (S1/obs listings; 336 / 45,666 role assignments)
  CENSUS6 (c)       every seed (TF1, TF2, m, t) inside every one of the 210,115 modules: miR1 = m,
                    G1 = t; no miR2 / G2 roles (a D1-D4 module cannot hold a co-transcribed pair)
Model edges (dyn_models.py, S0/s0c_model6_signed_edge_list.csv): core TF1->miR1, miR1->TF1,
miR1->G1, TF1->G1; GG G1-G2, TF1->G2, miR1->G2; MM miR1-miR2, TF1->miR2, miR2->G1, miR2->G2,
miR2->TF1; TT TF1->TF2, TF2->G1, TF2->G2, TF2->miR1, TF2->miR2.  Only edges between roles the object
has are evaluated.
Definitions
  sign-resolved core   TF1->miR1 and TF1->G1 present and signed -> configuration (dyn_models.py
                       identifiers): COMP_I1 (+,+), COMP_C2 (-,+), COMP_C4 (+,-), COMP_I3 (-,-)
  loop                 TF1 -> TF2 -> miR1 -| TF1; fully signed when TF1->TF2 and TF2->miR1 are
                       present and signed; sign = s(TF1->TF2) * s(TF2->miR1) * (-1)
  exact match          every model edge between the object's roles present, with the family's sign on
                       every directed arc (the unsigned G1-G2 and miR1-miR2 links need only be present):
                       COMP_C2_toggle n6; COMP_I1 n6 (TF1->miR1, TF1->miR2 positive); S3a (TF2 -| miR1,
                       miR2); S3b (TF2 -| G1, G2); S3c (TF1 -| TF2)
  fully signed config  every compulsory TF-sourced arc of the object signed (compulsory: BHAT6 T1->M1,
                       T1->T2, T2->T1, T2->G2; MODEL6 FULL the nine TF arcs; MODEL6 DEF TF1->miR1,
                       TF1->G1, TF1->TF2, TF2->G1, TF2->miR1; CENSUS6 TF1->m, TF1->t, TF1->TF2,
                       TF2->m, TF2->t); 'at least one unsigned' counts the same arcs
(b) sign permutation: the known signs are shuffled within each paper edge class (TF->Gene, TF->TF,
    TF->miRNA; Gene-typed TF_target regulators included in their class), 1,000 times (numpy seed
    20250908); miRNA repression fixed; unsigned arcs stay unsigned.
"""
import os, csv, gzip, itertools, collections, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
S1 = os.path.join(HERE, '..', 'S1')
names = open(f'{S1}/node_names.txt').read().split('\n')[:587]
L = open(f'{S1}/graph_nolegacy_null.txt').read().rstrip('\n').split('\n')
nv, ne, _ = map(int, L[0].split()); ty = list(map(int, L[1].split()))
labs = list(map(int, open(f'{S1}/graph_nolegacy_labels.txt').read().split()))
arcs = [tuple(map(int, l.split()[:2])) for l in L[2:2 + ne]]
MIR, TF, GEN = 0, 1, 2
TYN = {0: 'miRNA', 1: 'TF', 2: 'Gene'}
AI = {a: i for i, a in enumerate(arcs)}

# ---------------- signs and evidence ----------------
S1tab = {(r['source'], r['target']): r for r in csv.DictReader(open(
    f'{REV}/results/v5/tables/TableS1_all_interactions.csv'))}
TRL = {(r['source'], r['target']): r for r in csv.DictReader(open(f'{REV}/data/layer_TF_target.tsv'), delimiter='\t')}
TML = {(r['source'], r['target']): r for r in csv.DictReader(open(f'{REV}/data/layer_TF_miRNA.tsv'), delimiter='\t')}
GGL = {}
for r in csv.DictReader(open(f'{REV}/data/layer_gene_gene.tsv'), delimiter='\t'):
    if r['evidence'] == 'STRING': GGL[frozenset((r['source'], r['target']))] = r
MML = {}
for r in csv.DictReader(open(f'{REV}/data/layer_miRNA_miRNA.tsv'), delimiter='\t'):
    MML[frozenset((r['miRNA_1'], r['miRNA_2']))] = r
MODE = {'Activation': 1, 'Repression': -1}
sign = np.zeros(ne, dtype=np.int8); prov = [''] * ne; evid = [''] * ne; tier = [''] * ne; pcls = [''] * ne
npmid = np.full(ne, -1)
for i, ((u, v), lb) in enumerate(zip(arcs, labs)):
    a, b = names[u], names[v]
    if lb == 1:
        sign[i] = -1; t = S1tab[(a, b)]; prov[i] = 'mechanistic (miRNA repression)'
        tier[i] = t['tier']; evid[i] = (t['databases'] or 'no validated record') + f"; tier {t['tier']}"
        pcls[i] = 'miRNA->TF' if ty[v] == TF else 'miRNA->Gene'
    elif lb == 0:
        m = S1tab[(a, b)]['transmir_mode']; sign[i] = MODE.get(m, 0)
        prov[i] = f'TransmiR mode {m}' if m else 'TransmiR (no mode)'
        if (a, b) in TML:
            npmid[i] = len([x for x in TML[(a, b)]['pmid'].split(';') if x]); evid[i] = f'TransmiR literature, {npmid[i]} PMID(s)'
        else:
            evid[i] = 'TransmiR literature (canonical miRNA name not in data/layer_TF_miRNA.tsv; PMIDs not traced)'
        pcls[i] = 'TF->miRNA'
    elif lb in (2, 3):
        m = (S1tab[(a, b)]['trrust_mode'] if lb == 2 else TRL[(a, b)]['mode'])
        sign[i] = MODE.get(m, 0); prov[i] = f'TRRUST mode {m}'
        if (a, b) in TRL:
            npmid[i] = int(TRL[(a, b)]['n_pmid']); evid[i] = f"TRRUST, {npmid[i]} PMID(s)" + ('' if lb == 2 else ' (census layer)')
        else:
            evid[i] = 'TRRUST (pair not in data/layer_TF_target.tsv; PMIDs not traced)'
        pcls[i] = 'TF->TF' if ty[v] == TF else 'TF->Gene'
    elif lb == 4:
        prov[i] = 'unsigned (gene-gene)'; evid[i] = 'deposit gene_gene edge'; pcls[i] = 'Gene->Gene'
    elif lb == 5:
        g = GGL[frozenset((a, b))]; prov[i] = 'unsigned (STRING)'
        evid[i] = f"STRING {g['score_or_mode']} ({g['string_network_type']})"; pcls[i] = 'Gene->Gene'
    elif lb == 6:
        g = MML[frozenset((a, b))]; prov[i] = 'unsigned (co-transcription)'
        evid[i] = f"10 kb cluster pair, {g['distance_bp']} bp"; pcls[i] = 'miRNA-miRNA'
arc_of = lambda u, v: AI.get((u, v), -1)
print('arcs whose PMIDs could not be traced:', collections.Counter(labs[i] for i in range(ne) if 'not traced' in evid[i]), flush=True)


def link(u, v):
    """undirected link (G1-G2 gene_gene or miR1-miR2 10 kb): an arc index in either direction."""
    for x in (arc_of(u, v), arc_of(v, u)):
        if x >= 0 and labs[x] in (4, 5, 6): return x
    return -1


# model edges: (name, src role, tgt role, kind) ; kind 'd' directed arc, 'l' undirected link
ME = [('TF1->miR1', 'TF1', 'miR1', 'd'), ('miR1->TF1', 'miR1', 'TF1', 'd'), ('miR1->G1', 'miR1', 'G1', 'd'),
      ('TF1->G1', 'TF1', 'G1', 'd'), ('G1-G2', 'G1', 'G2', 'l'), ('TF1->G2', 'TF1', 'G2', 'd'),
      ('miR1->G2', 'miR1', 'G2', 'd'), ('miR1-miR2', 'miR1', 'miR2', 'l'), ('TF1->miR2', 'TF1', 'miR2', 'd'),
      ('miR2->G1', 'miR2', 'G1', 'd'), ('miR2->G2', 'miR2', 'G2', 'd'), ('miR2->TF1', 'miR2', 'TF1', 'd'),
      ('TF1->TF2', 'TF1', 'TF2', 'd'), ('TF2->G1', 'TF2', 'G1', 'd'), ('TF2->G2', 'TF2', 'G2', 'd'),
      ('TF2->miR1', 'TF2', 'miR1', 'd'), ('TF2->miR2', 'TF2', 'miR2', 'd')]
MEI = {e[0]: k for k, e in enumerate(ME)}
C2 = {'TF1->miR1': -1, 'TF1->miR2': -1, 'TF1->G1': 1, 'TF1->G2': 1, 'TF1->TF2': 1, 'TF2->G1': 1, 'TF2->G2': 1,
      'TF2->miR1': 1, 'TF2->miR2': 1, 'miR1->TF1': -1, 'miR1->G1': -1, 'miR1->G2': -1, 'miR2->G1': -1,
      'miR2->G2': -1, 'miR2->TF1': -1}
FAM = {'COMP_C2_toggle_n6': C2,
       'COMP_I1_n6': {**C2, 'TF1->miR1': 1, 'TF1->miR2': 1},
       'S3a_TF2-|miR': {**C2, 'TF2->miR1': -1, 'TF2->miR2': -1},
       'S3b_TF2-|G1': {**C2, 'TF2->G1': -1, 'TF2->G2': -1},
       'S3c_TF1-|TF2': {**C2, 'TF1->TF2': -1}}
COMPULSORY = {'BHAT6_composite': None,        # handled separately (Bhat arcs are not all model edges)
              'MODEL6_FULL': ['TF1->miR1', 'TF1->miR2', 'TF1->G1', 'TF1->G2', 'TF1->TF2', 'TF2->G1', 'TF2->G2', 'TF2->miR1', 'TF2->miR2'],
              'MODEL6_DEF': ['TF1->miR1', 'TF1->G1', 'TF1->TF2', 'TF2->G1', 'TF2->miR1'],
              'CENSUS6_c': ['TF1->miR1', 'TF1->G1', 'TF1->TF2', 'TF2->miR1', 'TF2->G1']}


def edge_arcs(roles):
    """arc index (or -1 absent, -2 role missing) for every model edge, given a role dict."""
    out = np.full(len(ME), -2, dtype=np.int64)
    for k, (nm, a, b, kind) in enumerate(ME):
        if roles.get(a) is None or roles.get(b) is None: continue
        out[k] = arc_of(roles[a], roles[b]) if kind == 'd' else link(roles[a], roles[b])
    return out


# ---------------- instances ----------------
def rd(p, cols):
    return [dict(zip(cols, map(int, l.split()))) for l in open(p) if l.strip()]
objs = {}
bh = rd(f'{S1}/obs/nolegacy_bhat6_composite_instances.tsv', ['miR1', 'miR2', 'TF1', 'TF2', 'G1', 'G2'])
objs['BHAT6_composite'] = bh
objs['MODEL6_FULL'] = rd(f'{S1}/obs/nolegacy_model6_full_instances.tsv', ['TF1', 'TF2', 'miR1', 'miR2', 'G1', 'G2'])
objs['MODEL6_DEF'] = rd(f'{S1}/obs/nolegacy_model6_def_instances.tsv', ['TF1', 'TF2', 'miR1', 'miR2', 'G1', 'G2'])
# CENSUS6: seeds inside every module (contracted adjacency as the census builder: drop miRNA->TF of reciprocal pairs)
A = set(arcs)
CO = {(u, v) for (u, v) in arcs if not (ty[u] == MIR and ty[v] == TF and (v, u) in A)}
seeds = collections.defaultdict(list)
for (a, m) in CO:
    if ty[a] != TF or ty[m] != MIR or (m, a) not in A: continue
    for b in range(nv):
        if b == a or ty[b] != TF or (a, b) not in CO or (b, m) not in CO: continue
        for t in range(nv):
            if t in (a, b, m) or ty[t] == MIR: continue
            if (a, t) in CO and (b, t) in CO and (m, t) in CO:
                seeds[frozenset((a, b, m, t))].append(dict(TF1=a, TF2=b, miR1=m, G1=t))
assert sum(len(v) for v in seeds.values()) == 1339
mods = [tuple(map(int, l.split())) for l in open(f'{S1}/obs/nolegacy_census6c_modules.tsv') if l.strip()]
cen = []; cen_mod = []
for mi, M in enumerate(mods):
    for sub in itertools.combinations(M, 4):
        for s in seeds.get(frozenset(sub), ()):
            cen.append(s); cen_mod.append(mi)
objs['CENSUS6_c'] = cen
cen_mod = np.array(cen_mod)
print({k: len(v) for k, v in objs.items()}, 'CENSUS6 modules', len(mods), flush=True)

EA = {k: np.array([edge_arcs(r) for r in v]) for k, v in objs.items()}
# BHAT6 compulsory TF-sourced arcs: T1->M1, T1->T2, T2->T1, T2->G2
BH_ARCS = np.array([[arc_of(r['TF1'], r['miR1']), arc_of(r['TF1'], r['TF2']), arc_of(r['TF2'], r['TF1']),
                     arc_of(r['TF2'], r['G2'])] for r in bh])
assert (BH_ARCS >= 0).all()
BH_NAMES = ['T1->M1', 'T1->T2', 'T2->T1', 'T2->G2']


def evaluate(S, with_conf=False):
    """S: sign vector over arcs.  Returns dict of per-object boolean/str arrays."""
    res = {}
    Sx = np.append(S, 0).astype(np.int64)          # index -1 -> 0 (absent is handled via masks)
    for k, E in EA.items():
        pres = E >= 0
        sg = np.where(pres, Sx[np.where(pres, E, len(S))], 0)          # 0 = absent or unsigned
        col = lambda nm: MEI[nm]
        tm, ty_ = sg[:, col('TF1->miR1')], sg[:, col('TF1->G1')]
        core_ok = (tm != 0) & (ty_ != 0)
        cfg = np.where(~core_ok, 'core not sign-resolved',
                       np.where(tm > 0, np.where(ty_ > 0, 'COMP_I1', 'COMP_C4'), np.where(ty_ > 0, 'COMP_C2', 'COMP_I3')))
        s12, s2m = sg[:, col('TF1->TF2')], sg[:, col('TF2->miR1')]
        loop_ok = (s12 != 0) & (s2m != 0)
        loop = np.where(~loop_ok, 'not fully signed', np.where(s12 * s2m > 0, 'negative', 'positive'))
        ex = {}
        for fam, sgn in FAM.items():
            ok = np.ones(len(E), bool)
            for kk, (nm, a, b, kind) in enumerate(ME):
                role_missing = E[:, kk] == -2
                if kind == 'l':
                    ok &= role_missing | (E[:, kk] >= 0)
                else:
                    ok &= role_missing | (sg[:, kk] == sgn[nm])
            ex[fam] = ok
        if k == 'BHAT6_composite':
            cs = Sx[BH_ARCS]
            fully = (cs != 0).all(1); anyun = (cs == 0).any(1)
            conf = np.array([' '.join(f'{n}:{"+" if x > 0 else "-"}' for n, x in zip(BH_NAMES, row)) if f else ''
                             for row, f in zip(cs, fully)]) if with_conf else None
        else:
            cols = [MEI[n] for n in COMPULSORY[k]]
            cs = sg[:, cols]; assert (E[:, cols] >= 0).all()
            fully = (cs != 0).all(1); anyun = (cs == 0).any(1)
            conf = np.array([' '.join(f'{n}:{"+" if x > 0 else "-"}' for n, x in zip(COMPULSORY[k], row)) if f else ''
                             for row, f in zip(cs, fully)]) if with_conf else None
        res[k] = dict(cfg=cfg, loop=loop, ex=ex, fully=fully, anyun=anyun, conf=conf)
    return res


R0 = evaluate(sign.astype(np.int64), with_conf=True)
rows = []; top = []
for k, r in R0.items():
    n = len(r['cfg'])
    add = lambda cat, val, mask: rows.append(dict(object=k, unit='role assignments', category=cat, value=val,
                                                  count=int(mask.sum()), n=n, pct=100 * mask.mean(),
                                                  modules=(len(set(cen_mod[mask])) if k == 'CENSUS6_c' else '')))
    for c in ('COMP_C2', 'COMP_I1', 'COMP_C4', 'COMP_I3', 'core not sign-resolved'):
        add('core configuration', c, r['cfg'] == c)
    for c in ('negative', 'positive', 'not fully signed'):
        add('TF1->TF2->miR1-|TF1 loop', c, r['loop'] == c)
    for fam in FAM: add('exact match', fam, r['ex'][fam])
    add('compulsory TF arcs', 'all signed', r['fully'])
    add('compulsory TF arcs', 'at least one unsigned', r['anyun'])
    cc = collections.Counter(r['conf'][r['fully']])
    for i, (c, v) in enumerate(cc.most_common(10)):
        top.append(dict(object=k, rank=i + 1, configuration=c, count=v, of_fully_signed=int(r['fully'].sum())))
with open(f'{HERE}/s5a_counts.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(f'{HERE}/s5a_top10_fully_signed_configurations.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(top[0])); w.writeheader(); w.writerows(top)
for x in rows: print(x)

# ---------------- (b) sign permutation null ----------------
rng = np.random.default_rng(20250908)
cls_idx = collections.defaultdict(list)
for i in range(ne):
    if sign[i] != 0 and labs[i] != 1: cls_idx[pcls[i]].append(i)
print('signed arcs permuted per class:', {k: len(v) for k, v in cls_idx.items()}, flush=True)
stat = lambda R: {k: (int((r['loop'] == 'negative').sum()), int(r['ex']['COMP_C2_toggle_n6'].sum()),
                      int(r['ex']['COMP_I1_n6'].sum()), int((r['cfg'] == 'COMP_C2').sum()), int((r['cfg'] == 'COMP_I1').sum()))
                  for k, r in R.items()}
obs = stat(R0)
NP = 1000; null = {k: np.zeros((NP, 5)) for k in obs}
for p in range(NP):
    S = sign.astype(np.int64).copy()
    for c, idx in cls_idx.items():
        idx = np.array(idx); S[idx] = S[idx][rng.permutation(len(idx))]
    st = stat(evaluate(S))
    for k in st: null[k][p] = st[k]
    if (p + 1) % 100 == 0: print(f'  permutation {p + 1}/{NP}', flush=True)
nb = []
for k in obs:
    for j, nm in enumerate(['negative TF1->TF2->miR1-|TF1 loop', 'COMP_C2_toggle n6 exact', 'COMP_I1 n6 exact',
                            'core configuration COMP_C2', 'core configuration COMP_I1']):
        x = null[k][:, j]; o = obs[k][j]
        nb.append(dict(object=k, statistic=nm, observed=o, null_mean=x.mean(), null_sd=x.std(ddof=1),
                       p_upper=(np.sum(x >= o) + 1) / (NP + 1), p_lower=(np.sum(x <= o) + 1) / (NP + 1), permutations=NP))
with open(f'{HERE}/s5b_sign_permutation_null.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(nb[0])); w.writeheader(); w.writerows(nb)
for x in nb: print(x)

# ---------------- (c) named circuits ----------------
TR = {'strong': 0, 'weak': 1, 'predicted_only': 2}
k = 'MODEL6_FULL'
sel = np.where(R0[k]['loop'] == 'negative')[0]
listed = k
if len(sel) == 0:
    k = 'BHAT6_composite'; sel = np.where(R0[k]['loop'] == 'negative')[0]; listed = 'BHAT6_composite (no MODEL6 instance qualifies)'
out_rows = []; circ = []
for i in sel:
    roles = objs[k][i]; E = EA[k][i]
    used = [(ME[kk][0], E[kk]) for kk in range(len(ME)) if E[kk] >= 0]
    if k == 'BHAT6_composite':
        extra = [('T2->T1', arc_of(roles['TF2'], roles['TF1'])), ('M2->G1 (Bhat)', arc_of(roles['miR2'], roles['G1'])),
                 ('T2->G2 (Bhat)', arc_of(roles['TF2'], roles['G2']))]
        used += [(n, a) for n, a in extra if a >= 0 and a not in {x for _, x in used}]
    mt = [TR[tier[a]] for _, a in used if labs[a] == 1]
    tp = [npmid[a] for _, a in used if labs[a] in (2, 3)]
    key = (max(mt) if mt else -1, -(min(tp) if tp else 0))
    cid = '|'.join(f"{r}={names[roles[r]]}" for r in ('TF1', 'TF2', 'miR1', 'miR2', 'G1', 'G2'))
    circ.append((key, cid, i, used))
circ.sort(key=lambda x: (x[0], x[1]))
for rank, (key, cid, i, used) in enumerate(circ, 1):
    for nm, a in used:
        u, v = arcs[a]
        out_rows.append(dict(rank=rank, circuit=cid, weakest_miRNA_tier=['strong', 'weak', 'predicted_only'][key[0]] if key[0] >= 0 else '',
                             min_TRRUST_pmids=-key[1], model_edge=nm, source=names[u], target=names[v],
                             layer={0: 'TF_miRNA', 1: 'miRNA_target', 2: 'TF_target', 3: 'TF_target (TRRUST census layer)',
                                    4: 'gene_gene (deposit)', 5: 'gene_gene (STRING)', 6: 'miRNA_miRNA (10 kb)'}[labs[a]],
                             sign={1: '+', -1: '-', 0: 'unsigned'}[int(sign[a])], sign_provenance=prov[a], evidence=evid[a]))
with open(f'{HERE}/s5c_named_circuits_edges.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0]) if out_rows else ['rank']); w.writeheader(); w.writerows(out_rows)
print('S5(c) listed object:', listed, 'circuits:', len(circ))
