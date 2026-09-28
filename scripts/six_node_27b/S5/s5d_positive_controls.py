#!/usr/bin/env python3
"""
S5(d): two experimentally characterised TF-miRNA circuits as positive controls (context, not claims).
Edges are taken from the analysis plan's descriptions only; absent edges are reported as absent, none is added.
  MYC/E2F - miR-17~92 (O'Donnell 2005; Woods 2007; Sylvestre 2007; Aguda 2008):
    MYC, E2F1, E2F3 -> miR-17, -18a, -19a, -20a, -19b, -92a (activation); miR-17, miR-20a -| E2F1;
    MYC -> E2F1 and E2F1 -> MYC (activation)
  ZEB/SNAIL - miR-200 (Bracken 2008; Burk 2008; Lu 2013):
    ZEB1, ZEB2 -| miR-200a, -200b, -200c, -141, -429; those miRNAs -| ZEB1, ZEB2; SNAI1 -> ZEB1
Presence, sign and evidence are read from the nolegacy census graph with the S5 sign rules
(s5_signed_circuits.py).  Containment is checked in the S1 observed listings (BHAT6 composite, MODEL6
FULL / DEF), the CENSUS6 (c)-with-core modules and the typed three-node cores (D1-D4 three-node FFLs of
the analysed network, INBOX_2026-09-26b/graphs/null_dep_nolegacy.txt).  Read-only.
"""
import os, csv, json, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
S1 = os.path.join(HERE, '..', 'S1')
names = open(f'{S1}/node_names.txt').read().split('\n')[:587]; ix = {n: i for i, n in enumerate(names)}
L = open(f'{S1}/graph_nolegacy_null.txt').read().rstrip('\n').split('\n')
nv, ne, _ = map(int, L[0].split()); ty = list(map(int, L[1].split()))
labs = list(map(int, open(f'{S1}/graph_nolegacy_labels.txt').read().split()))
arcs = [tuple(map(int, l.split()[:2])) for l in L[2:2 + ne]]; AI = {a: i for i, a in enumerate(arcs)}
S1tab = {(r['source'], r['target']): r for r in csv.DictReader(open(
    f'{REV}/results/v5/tables/TableS1_all_interactions.csv'))}
TRL = {(r['source'], r['target']): r for r in csv.DictReader(open(f'{REV}/data/layer_TF_target.tsv'), delimiter='\t')}
MODE = {'Activation': '+', 'Repression': '-'}
LAYER = {0: 'TF_miRNA (analysed network)', 1: 'miRNA_target (analysed network)', 2: 'TF_target (analysed network)',
         3: 'TF_target (TRRUST census layer)', 4: 'gene_gene (deposit)', 5: 'gene_gene (STRING census layer)',
         6: 'miRNA_miRNA (10 kb census layer)'}


def edge_info(a, b):
    if a not in ix or b not in ix: return dict(present='no (node not in the network: ' + ', '.join(x for x in (a, b) if x not in ix) + ')')
    k = AI.get((ix[a], ix[b]))
    if k is None: return dict(present='no')
    lb = labs[k]; d = dict(present='yes', layer=LAYER[lb])
    if lb == 1:
        t = S1tab[(a, b)]; d.update(sign='-', sign_provenance='mechanistic (miRNA repression)', evidence=f"{t['databases'] or 'no validated record'}", tier=t['tier'])
    elif lb == 0:
        m = S1tab[(a, b)]['transmir_mode']; d.update(sign=MODE.get(m, 'unsigned'), sign_provenance=f'TransmiR mode {m}', evidence='TransmiR (literature)', tier='')
    elif lb in (2, 3):
        m = S1tab[(a, b)]['trrust_mode'] if lb == 2 else TRL[(a, b)]['mode']
        d.update(sign=MODE.get(m, 'unsigned'), sign_provenance=f'TRRUST mode {m}', evidence=f"TRRUST, {TRL.get((a, b), {}).get('n_pmid', '?')} PMID(s)", tier='')
    else:
        d.update(sign='unsigned', sign_provenance='unsigned layer', evidence=LAYER[lb], tier='')
    return d


M17 = ['hsa-miR-17', 'hsa-miR-18a', 'hsa-miR-19a', 'hsa-miR-20a', 'hsa-miR-19b', 'hsa-miR-92a']
M200 = ['hsa-miR-200a', 'hsa-miR-200b', 'hsa-miR-200c', 'hsa-miR-141', 'hsa-miR-429']
CIRC = {'MYC/E2F-miR-17~92': [(tf, m, '+') for tf in ('MYC', 'E2F1', 'E2F3') for m in M17] +
                             [('hsa-miR-17', 'E2F1', '-'), ('hsa-miR-20a', 'E2F1', '-'), ('MYC', 'E2F1', '+'), ('E2F1', 'MYC', '+')],
        'ZEB/SNAIL-miR-200': [(z, m, '-') for z in ('ZEB1', 'ZEB2') for m in M200] +
                             [(m, z, '-') for m in M200 for z in ('ZEB1', 'ZEB2')] + [('SNAI1', 'ZEB1', '+')]}
rows = []
for c, E in CIRC.items():
    for a, b, lit in E:
        d = edge_info(a, b)
        rows.append(dict(circuit=c, source=a, target=b, literature_sign=lit, **{k: d.get(k, '') for k in
                         ('present', 'layer', 'sign', 'sign_provenance', 'evidence', 'tier')}))
with open(f'{HERE}/s5d_positive_control_edges.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
out = []
P = lambda *a: out.append(' '.join(str(x) for x in a))
for c in CIRC:
    R = [r for r in rows if r['circuit'] == c]
    P(f'{c}: {sum(r["present"] == "yes" for r in R)} of {len(R)} edges present; '
      f'network sign agrees with the literature sign for {sum(r["sign"] == r["literature_sign"] for r in R)}; '
      f'present but unsigned {sum(r["present"] == "yes" and r["sign"] == "unsigned" for r in R)}; '
      f'present with the opposite sign {sum(r["present"] == "yes" and r["sign"] in "+-" and r["sign"] != r["literature_sign"] for r in R)}')

# ---- containment
def rd(p, cols): return [dict(zip(cols, map(int, l.split()))) for l in open(p) if l.strip()]
BH = rd(f'{S1}/obs/nolegacy_bhat6_composite_instances.tsv', ['miR1', 'miR2', 'TF1', 'TF2', 'G1', 'G2'])
MF = rd(f'{S1}/obs/nolegacy_model6_full_instances.tsv', ['TF1', 'TF2', 'miR1', 'miR2', 'G1', 'G2'])
MD = rd(f'{S1}/obs/nolegacy_model6_def_instances.tsv', ['TF1', 'TF2', 'miR1', 'miR2', 'G1', 'G2'])
CM = [set(map(int, l.split())) for l in open(f'{S1}/obs/nolegacy_census6c_modules.tsv') if l.strip()]
# typed three-node cores: D1-D4 three-node FFLs of the analysed network (count3 enumeration)
G = f'{REV}/INBOX_2026-09-26b/graphs/null_dep_nolegacy.txt'
Lg = open(G).read().rstrip('\n').split('\n'); nm2 = json.load(open(G + '.meta.json'))['nodes']; assert nm2 == names
t2 = list(map(int, Lg[1].split())); Ag = {tuple(map(int, l.split()[:2])) for l in Lg[2:]}
red = {(u, v) for (u, v) in Ag if not (t2[u] == 0 and t2[v] == 1 and (v, u) in Ag)}
out_ = collections.defaultdict(set); in_ = collections.defaultdict(set)
for u, v in red: out_[u].add(v); in_[v].add(u)
cores = []
for (u, v) in red:
    if u in out_[v]: continue
    for t in out_[u] & out_[v]:
        if t not in in_[u] and t not in in_[v]: cores.append((u, v, t))
assert len(cores) == 1774
KEYS = {'MYC/E2F-miR-17~92': dict(roles=[dict(TF1='E2F1', TF2='MYC', miR1=m) for m in ('hsa-miR-17', 'hsa-miR-20a')],
                                  sets=[{'E2F1', 'MYC', m} for m in ('hsa-miR-17', 'hsa-miR-20a')],
                                  pairs=[{'E2F1', m} for m in ('hsa-miR-17', 'hsa-miR-20a')]),
        'ZEB/SNAIL-miR-200': dict(roles=[dict(TF1='ZEB2', miR1=m) for m in M200],
                                  sets=[{'ZEB2', 'SNAI1', m} for m in M200], pairs=[{'ZEB2', m} for m in M200])}
for c, K in KEYS.items():
    sets = [{ix[x] for x in s} for s in K['sets']]; pairs = [{ix[x] for x in s} for s in K['pairs']]
    roles = [{k: ix[v] for k, v in r.items()} for r in K['roles']]
    def cnt(inst):
        byrole = sum(any(all(i[k] == v for k, v in r.items()) for r in roles) for i in inst)
        s_ = sum(any(s <= set(i.values()) for s in sets) for i in inst); p_ = sum(any(p <= set(i.values()) for p in pairs) for i in inst)
        return byrole, s_, p_
    P(f'{c}: key roles {K["roles"]}; key node sets {K["sets"]}; key pairs {K["pairs"]}')
    for nm, inst in (('BHAT6 composite', BH), ('MODEL6 FULL', MF), ('MODEL6 DEF', MD)):
        a, b, d = cnt(inst)
        P(f'   {nm}: role assignments with the key roles {a}; containing a key node set {b}; containing a key pair {d} (of {len(inst)})')
    P(f'   CENSUS6 (c)-with-core modules containing a key node set {sum(any(s <= M for s in sets) for M in CM)}; '
      f'containing a key pair {sum(any(p <= M for p in pairs) for M in CM)} (of {len(CM)})')
    tc = [(names[u], names[v], names[t]) for (u, v, t) in cores]
    P(f'   typed three-node cores (analysed network, 1,774) containing a key pair: '
      f'{sum(any(p <= {ix[x] for x in t} for p in pairs) for t in tc)}; containing a key node set: {sum(any(s <= {ix[x] for x in t} for s in sets) for t in tc)}')
    ex = [t for t in tc if any(p <= {ix[x] for x in t} for p in pairs)]
    P(f'   those cores (u -> v -> t): {ex}')
open(f'{HERE}/s5d_positive_controls.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
for r in rows:
    if r['present'] == 'yes': print(r)

# ---- added: sign configuration of the typed cores that contain a key pair, and CENSUS6 role-level counts
def sgn(a, b):
    d = edge_info(a, b); return d.get('sign', 'absent') if d.get('present') == 'yes' else 'absent'
P('sign configuration of the typed three-node cores containing a key pair (TF -> miRNA, TF -> target; miRNA -> target is -):')
CFG = {('+', '+'): 'COMP_I1', ('-', '+'): 'COMP_C2', ('+', '-'): 'COMP_C4', ('-', '-'): 'COMP_I3'}
for (u, v, t) in [(names[a], names[b], names[c]) for (a, b, c) in cores]:
    if (u, v) in (('E2F1', 'hsa-miR-17'), ('E2F1', 'hsa-miR-20a')):
        tm, tt = sgn(u, v), sgn(u, t)
        P(f'   {u} -> {v} ({tm}), {u} -> {t} ({tt}), {v} -| {t}; reciprocal {v} -| {u}: {sgn(v, u)} -> '
          f'{CFG.get((tm, tt), "not sign-resolved")}')
A = set(arcs); CO = {(u, v) for (u, v) in arcs if not (ty[u] == 0 and ty[v] == 1 and (v, u) in A)}
def census_seeds(TF1, TF2, m):
    a, b, mm = ix[TF1], ix[TF2], ix[m]
    if not ((a, mm) in CO and (mm, a) in A and (a, b) in CO and (b, mm) in CO): return []
    return [names[t] for t in range(nv) if t not in (a, b, mm) and ty[t] != 0 and (a, t) in CO and (b, t) in CO and (mm, t) in CO]
for m in ('hsa-miR-17', 'hsa-miR-20a'):
    ts = census_seeds('E2F1', 'MYC', m)
    P(f'CENSUS6 seeds with TF1=E2F1, TF2=MYC, m={m}: {len(ts)} (targets t: {ts}); modules containing such a seed: '
      f'{sum(any({ix["E2F1"], ix["MYC"], ix[m], ix[t]} <= M for t in ts) for M in CM)}')
    P(f'   loop E2F1 -> MYC -> {m} -| E2F1 signs: E2F1->MYC {sgn("E2F1", "MYC")} ({edge_info("E2F1", "MYC").get("layer", "")}), '
      f'MYC->{m} {sgn("MYC", m)}, {m}->E2F1 {sgn(m, "E2F1")}')
for m in M200:
    P(f'ZEB2 / {m}: reciprocal ZEB2 -| {m} ({sgn("ZEB2", m)}) and {m} -| ZEB2 ({sgn(m, "ZEB2")}); common targets of ZEB2 and {m} '
      f'in the analysed network: {sorted(names[t] for t in out_[ix["ZEB2"]] & out_[ix[m]])}')
open(f'{HERE}/s5d_positive_controls.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out[-14:]))
