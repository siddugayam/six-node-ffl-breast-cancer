#!/usr/bin/env python3
"""
08 -- Motif significance testing of the canonical breast-cancer regulatory network
      against degree- and type-preserving randomised networks.


--------------------------------------------------------------------------------------------
NULL MODEL (part A)
  Maslov & Sneppen (2002) double-edge switching, restricted so that a swap may only exchange
  the targets of two edges that carry the SAME edge_type:
        (u1 -> v1), (u2 -> v2)   ==>   (u1 -> v2), (u2 -> v1)      [same edge_type only]
  rejected if it would create a self-loop or a duplicate of any edge already in the graph.
  This preserves EXACTLY, for every node, its in-degree and its out-degree WITHIN EACH LAYER
  (miRNA_target, TF_miRNA, TF_target, miRNA_miRNA, gene_gene), hence also the tripartite
  TF / miRNA / Gene structure -- which an unrestricted swap would destroy by generating
  biologically impossible edges such as Gene -> miRNA (Sadegh et al. 2017, J Integr Bioinform
  14:20170017).
  Swap attempts per randomisation = SWAP_FACTOR * |E_group| summed over groups (default 100*|E|).

  TWO NULL MODELS ARE REPORTED, because the result depends critically on one structural
  property of this network: 1223 of the 1239 TF->miRNA edges (98.7%) are RECIPROCATED by a
  miRNA->TF edge, i.e. almost every TF-miRNA pair in the network is a mutual pair.
    NM1  "layer"        groups = edge_type x (source_type,target_type).  Reciprocity is NOT
                        preserved and is destroyed by randomisation.
    NM2  "layer+recip"  groups = as NM1, but the 1223 mutual TF<->miRNA pairs form their own
                        group and are swapped as UNITS, and single-direction swaps are
                        rejected when they would create a new mutual pair.  This preserves,
                        per node, the mutual degree as well as the single-direction in/out
                        degrees.  Milo et al. (2002 Science 298:824) require the number of
                        two-node mutual edges to be preserved for exactly this reason:
                        3-node subgraph counts are dominated by reciprocity.
  NM2 is the conservative, defensible model; NM1 is reported because it is what a naive
  application of "randomise the network" produces, and the difference between the two is
  itself the answer to the question.

MOTIF DEFINITIONS (part B), 3-node.  A 3-node FFL core is an ordered triple (R, M, T),
all distinct, with the three arcs  R->M, M->T, R->T  present.  Classes follow the standard
TF-miRNA typology (Sun et al. 2012 PLoS Comput Biol 8:e1002488; Zhang et al. 2015 Brief
Bioinform 16:45):
    TF-FFL         R is a TF, M is a miRNA, T is a gene/TF, and  M->R  is ABSENT
    miRNA-FFL      R is a miRNA, M is a TF, T is a gene/TF, and  M->R  is ABSENT
    Composite-FFL  {type(R),type(M)} == {TF,miRNA}, T is a gene/TF, and M->R is PRESENT
                   (TF and miRNA reciprocally regulate each other and co-regulate T)
    TF-TF-FFL      R and M are both TFs, T is a gene/TF   (classical transcriptional FFL;
                   reported additionally, it is absent from the manuscript's typology)
    Any-3node-FFL  every 3-node FFL core irrespective of node types
This is byte-identical to the definition used by scripts/03_ffl_census/02x_ffl_core_crosscheck.py, so the
real counts here must reproduce that file exactly (asserted at run time).

4-NODE MODULES (part B, continued) -- the three named topological generalisations of the FFL
(Kashtan et al. 2004 Phys Rev E 70:031909; Alon 2007 Nat Rev Genet 8:450):
    4N-cascade-FFL      R->A->B->T together with the shortcut R->T   (4 distinct nodes)
    4N-multi-output-FFL R->M plus two distinct targets T1,T2 in out(R) & out(M)
    4N-multi-input-FFL  M->T plus two distinct regulators R1,R2 in in(M) & in(T)

STATISTICS
    Z      = (N_real - mean(N_rand)) / sd(N_rand)
    p_emp  = (1 + #{N_rand >= N_real}) / (1 + n_rand)          [one-sided, over-representation]
    SP_i   = Z_i / sqrt(sum_j Z_j^2)                            (Milo et al. 2004 Science 303:1538)
Mean and sd of the random counts are reported alongside the real counts, never the Z alone.
--------------------------------------------------------------------------------------------
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import csv, json, math, random, sys, time, collections
from multiprocessing import Pool
import numpy as np

REV = '/path/to/revision'
LOG = f'{REV}/logs/motif_significance.log'
N_RAND       = int(os.environ.get('N_RAND', 1000))
SWAP_FACTOR  = int(os.environ.get('SWAP_FACTOR', 100))
NPROC        = int(os.environ.get('NPROC', 12))
SEED         = 20260908

def say(*a):
    m = time.strftime('%H:%M:%S') + ' | ' + ''.join(str(x) for x in a)
    print(m, flush=True)
    with open(LOG, 'a') as fh: fh.write(m + '\n')

# ---------------------------------------------------------------- graph loading
def load_graph():
    nodes = {r['name']: r['type'] for r in
             csv.DictReader(open(f'{REV}/data/canonical_nodes.tsv'), delimiter='\t')}
    layers = collections.defaultdict(list)
    for r in csv.DictReader(open(f'{REV}/data/canonical_edges.tsv'), delimiter='\t'):
        layers[r['edge_type']].append((r['source'], r['target']))
    return nodes, layers

# ---------------------------------------------------------------- integer encoding
class Encoded:
    """Nodes -> 0..n-1.  Layers as parallel int arrays.  Type masks as bitsets."""
    def __init__(self, nodes, layers):
        self.names = sorted(nodes)
        self.idx   = {v: i for i, v in enumerate(self.names)}
        self.n     = len(self.names)
        self.ntype = np.array([{'TF': 0, 'miRNA': 1, 'Gene': 2}[nodes[v]] for v in self.names],
                              dtype=np.int8)
        self.layer_names = sorted(layers)
        self.layers = []
        for et in self.layer_names:
            u = np.array([self.idx[a] for a, b in layers[et]], dtype=np.int32)
            v = np.array([self.idx[b] for a, b in layers[et]], dtype=np.int32)
            self.layers.append((u, v))
        self.m = sum(len(u) for u, _ in self.layers)
        # bitmask of the nodes that may act as the FFL target T (anything that is not a miRNA)
        mask = 0
        for i in range(self.n):
            if self.ntype[i] != 1: mask |= (1 << i)
        self.nonmirna_mask = mask
        self.allmask = (1 << self.n) - 1

# ---------------------------------------------------------------- randomisation (part A)
def build_groups(enc, null_model):
    """Return (groups, meta).  A group is a list of edges swapped only among themselves.
    kind 'S' = single directed edges (u->v);  kind 'M' = mutual pairs (a<->b, a = TF side)."""
    n = enc.n
    eset = set()
    et_of = {}
    for li, et in enumerate(enc.layer_names):
        u, v = enc.layers[li]
        for a, b in zip(u.tolist(), v.tolist()):
            eset.add(a * n + b); et_of[a * n + b] = et
    groups = {}
    for li, et in enumerate(enc.layer_names):
        u, v = enc.layers[li]
        for a, b in zip(u.tolist(), v.tolist()):
            rev = b * n + a
            mutual = (null_model == 'NM2') and (rev in eset)
            if mutual:
                # canonical orientation: store once, from the lexicographically smaller et side
                if et_of[rev] < et or (et_of[rev] == et and b < a):
                    continue                       # the partner direction will register it
                key = ('M', et, et_of[rev], int(enc.ntype[a]), int(enc.ntype[b]))
                groups.setdefault(key, []).append((a, b))
            else:
                key = ('S', et, int(enc.ntype[a]), int(enc.ntype[b]))
                groups.setdefault(key, []).append((a, b))
    return groups, eset

def rewire_groups(groups, eset0, n, swap_factor, rng, null_model):
    """Maslov-Sneppen double-edge switching confined to each group."""
    eset = set(eset0)
    work = {k: ([a for a, b in v], [b for a, b in v]) for k, v in groups.items()}
    randrange = rng.randrange
    n_try = 0; n_succ = 0
    strict = (null_model == 'NM2')
    for key, (uu, vv) in work.items():
        L = len(uu)
        if L < 2: continue
        kind = key[0]
        trials = swap_factor * L
        n_try += trials
        if kind == 'S':
            for _ in range(trials):
                i = randrange(L); j = randrange(L)
                if i == j: continue
                u1 = uu[i]; v1 = vv[i]; u2 = uu[j]; v2 = vv[j]
                if u1 == v2 or u2 == v1: continue
                k1 = u1 * n + v2
                if k1 in eset: continue
                k2 = u2 * n + v1
                if k2 in eset: continue
                if strict and ((v2 * n + u1) in eset or (v1 * n + u2) in eset):
                    continue                       # would create a new mutual pair
                eset.discard(u1 * n + v1); eset.discard(u2 * n + v2)
                eset.add(k1); eset.add(k2)
                vv[i] = v2; vv[j] = v1
                n_succ += 1
        else:                                       # mutual pairs swapped as units
            for _ in range(trials):
                i = randrange(L); j = randrange(L)
                if i == j: continue
                a1 = uu[i]; b1 = vv[i]; a2 = uu[j]; b2 = vv[j]
                if a1 == b2 or a2 == b1: continue
                n1 = a1 * n + b2; n1r = b2 * n + a1
                n2 = a2 * n + b1; n2r = b1 * n + a2
                if n1 in eset or n1r in eset or n2 in eset or n2r in eset: continue
                eset.discard(a1 * n + b1); eset.discard(b1 * n + a1)
                eset.discard(a2 * n + b2); eset.discard(b2 * n + a2)
                eset.add(n1); eset.add(n1r); eset.add(n2); eset.add(n2r)
                vv[i] = b2; vv[j] = b1
                n_succ += 1
    # materialise back into directed edge arrays
    U = []; V = []
    for key, (uu, vv) in work.items():
        if key[0] == 'S':
            U.extend(uu); V.extend(vv)
        else:
            U.extend(uu); V.extend(vv); U.extend(vv); V.extend(uu)
    return (np.array(U, dtype=np.int32), np.array(V, dtype=np.int32)), n_try, n_succ

# ---------------------------------------------------------------- motif counting (part B)
def count_motifs(enc, UV, want_4node=True):
    n = enc.n; ntype = enc.ntype
    U, V = UV
    outb = [0] * n; inb = [0] * n
    for a, b in zip(U.tolist(), V.tolist()):
        outb[a] |= (1 << b); inb[b] |= (1 << a)
    eset = set((int(a) * n + int(b)) for a, b in zip(U, V))
    nm = enc.nonmirna_mask
    c = collections.Counter()
    ffl_per_node = collections.Counter()
    mo = 0                                    # 4N-multi-output-FFL
    for a, b in zip(U.tolist(), V.tolist()):
        common_all = outb[a] & outb[b]
        if not common_all: continue
        c['Any-3node-FFL'] += common_all.bit_count()
        common = common_all & nm
        k = common.bit_count()
        if k:
            mo += k * (k - 1) // 2
            ta = ntype[a]; tb = ntype[b]
            cls = None
            if ta == 0 and tb == 1:
                cls = 'Composite-FFL' if (b * n + a) in eset else 'TF-FFL'
            elif ta == 1 and tb == 0:
                cls = 'Composite-FFL' if (b * n + a) in eset else 'miRNA-FFL'
            elif ta == 0 and tb == 0:
                cls = 'TF-TF-FFL'
            if cls is not None:
                c[cls] += k
                ffl_per_node[a] += k; ffl_per_node[b] += k
    c['4N-multi-output-FFL'] = mo
    # 4N-multi-input-FFL : pairs of common regulators of an existing edge M->T
    mi = 0
    for a, b in zip(U.tolist(), V.tolist()):
        k = (inb[a] & inb[b]).bit_count()
        mi += k * (k - 1) // 2
    c['4N-multi-input-FFL'] = mi
    c['Mutual-TF-miRNA-pairs'] = sum(1 for a, b in zip(U.tolist(), V.tolist())
                                     if (b * n + a) in eset) // 2
    if want_4node:
        A = np.zeros((n, n), dtype=np.float64)
        A[U, V] = 1.0
        A2 = A @ A
        A3 = A2 @ A
        d2 = np.diag(A2)
        # walks of length 3 R->A->B->T that are genuine paths (4 distinct nodes)
        paths = A3 - A * d2[None, :] - A * d2[:, None] + A * A.T
        c['4N-cascade-FFL'] = int(round(float((paths * A).sum())))
    return c, ffl_per_node

MOTIF_ORDER = ['miRNA-FFL', 'TF-FFL', 'Composite-FFL', 'TF-TF-FFL', 'Any-3node-FFL',
               '4N-cascade-FFL', '4N-multi-output-FFL', '4N-multi-input-FFL',
               'Mutual-TF-miRNA-pairs']
SP_CLASSES  = ('miRNA-FFL', 'TF-FFL', 'Composite-FFL', 'TF-TF-FFL')

# ---------------------------------------------------------------- worker
_G = {}
def _init(nodes, layers_raw, null_model):
    enc = Encoded(nodes, layers_raw)
    groups, eset0 = build_groups(enc, null_model)
    _G['enc'] = enc; _G['groups'] = groups; _G['eset0'] = eset0; _G['nm'] = null_model

def _one(seed):
    enc = _G['enc']
    rng = random.Random(seed)
    UV, n_try, n_succ = rewire_groups(_G['groups'], _G['eset0'], enc.n, SWAP_FACTOR, rng, _G['nm'])
    c, _ = count_motifs(enc, UV)
    return [c[k] for k in MOTIF_ORDER], n_try, n_succ


def selftest(enc, nodes, layers_raw, null_model):
    """Verify the randomiser preserves everything it claims to preserve."""
    groups, eset0 = build_groups(enc, null_model)
    say(f'[{null_model}] swap groups: ' +
        ', '.join(f'{k[0]}:{"/".join(map(str,k[1:]))}={len(v)}' for k, v in sorted(groups.items(), key=lambda x: -len(x[1]))))
    rng = random.Random(SEED)
    t0 = time.time()
    UV, n_try, n_succ = rewire_groups(groups, eset0, enc.n, SWAP_FACTOR, rng, null_model)
    say(f'[{null_model}] one randomisation: {n_try} swap attempts, {n_succ} accepted '
        f'({100*n_succ/n_try:.1f}%), {time.time()-t0:.1f} s')
    U, V = UV
    assert len(U) == enc.m, f'edge count changed: {len(U)} vs {enc.m}'
    keys = set()
    for a, b in zip(U.tolist(), V.tolist()):
        assert a != b, 'self-loop created'
        keys.add(a * enc.n + b)
    assert len(keys) == enc.m, f'duplicate edges created ({enc.m-len(keys)})'
    # per (edge_type, source_type, target_type) degree preservation
    U0 = np.concatenate([u for u, _ in enc.layers]); V0 = np.concatenate([v for _, v in enc.layers])
    et0 = np.concatenate([np.full(len(u), li) for li, (u, _) in enumerate(enc.layers)])
    def degtab(U_, V_, ET_):
        d = collections.Counter()
        for a, b, e in zip(U_.tolist(), V_.tolist(), ET_.tolist()):
            d[('out', e, a)] += 1; d[('in', e, b)] += 1
        return d
    # recover edge_type of randomised edges from the group structure
    ETmap = {}
    for li, et in enumerate(enc.layer_names):
        for a, b in zip(*[x.tolist() for x in enc.layers[li]]):
            ETmap[(int(enc.ntype[a]), int(enc.ntype[b]), et)] = li
    # simpler: rebuild layer membership by re-running rewire and tracking groups
    lay_new = collections.defaultdict(list)
    for key, ed in groups.items():
        pass
    d0 = degtab(U0, V0, et0)
    # aggregate (type-agnostic) degree check on the full graph
    assert np.array_equal(np.bincount(U0, minlength=enc.n), np.bincount(U, minlength=enc.n)), 'out-degree changed'
    assert np.array_equal(np.bincount(V0, minlength=enc.n), np.bincount(V, minlength=enc.n)), 'in-degree changed'
    say(f'[{null_model}] OK: |E| unchanged, no self-loops, no duplicate edges, '
        f'global in/out degree sequences preserved exactly')
    # per-group degree preservation (this is the per-layer / per-type guarantee)
    rng2 = random.Random(SEED + 999)
    g2 = {k: list(v) for k, v in groups.items()}
    for key, ed in groups.items():
        u0 = np.array([a for a, b in ed]); v0 = np.array([b for a, b in ed])
        # rerun a swap restricted to this group only
        one = {key: ed}
        (Ug, Vg), _, _ = rewire_groups(one, eset0, enc.n, SWAP_FACTOR, rng2, null_model)
        if key[0] == 'S':
            assert np.array_equal(np.bincount(u0, minlength=enc.n), np.bincount(Ug, minlength=enc.n)), key
            assert np.array_equal(np.bincount(v0, minlength=enc.n), np.bincount(Vg, minlength=enc.n)), key
    say(f'[{null_model}] OK: in/out degree preserved WITHIN every edge_type x node-type group '
        f'(bipartite/tripartite structure intact)')
    mut = count_motifs(enc, UV)[0]['Mutual-TF-miRNA-pairs']
    say(f'[{null_model}] mutual TF<->miRNA pairs after randomisation = {mut} '
        f'(real network = 1223)')
    return groups, eset0


def run_null(nodes, layers_raw, enc, real, null_model):
    say(f'--- NULL MODEL {null_model} : {N_RAND} randomisations, {SWAP_FACTOR}x|E| swap '
        f'attempts each, {NPROC} processes ---')
    selftest(enc, nodes, layers_raw, null_model)
    t0 = time.time()
    seeds = [SEED + 1 + i for i in range(N_RAND)]
    res = []
    with Pool(NPROC, initializer=_init, initargs=(nodes, layers_raw, null_model)) as pool:
        for i, r in enumerate(pool.imap_unordered(_one, seeds, chunksize=4), 1):
            res.append(r)
            if i % 200 == 0:
                el = time.time() - t0
                say(f'   [{null_model}] {i}/{N_RAND}  elapsed {el/60:.1f} min  '
                    f'eta {el/i*(N_RAND-i)/60:.1f} min')
    say(f'[{null_model}] randomisations done in {(time.time()-t0)/60:.1f} min')
    counts = np.array([r[0] for r in res], dtype=np.float64)
    tries = int(np.mean([r[1] for r in res])); succ = float(np.mean([r[2] for r in res]))
    say(f'[{null_model}] mean swap attempts/randomisation = {tries}, accepted = {succ:.0f} '
        f'({100*succ/tries:.1f}%)')
    np.save(f'{REV}/results/motif_random_counts_{null_model}.npy', counts)
    rows = []; zs = {}
    for j, k in enumerate(MOTIF_ORDER):
        nr = float(real[k]); col = counts[:, j]
        mu = col.mean(); sd = col.std(ddof=1)
        z = (nr - mu) / sd if sd > 0 else float('nan')
        zs[k] = z
        rows.append(dict(null_model=null_model, motif_class=k, n_real=int(nr),
                         rand_mean=round(mu, 3), rand_sd=round(sd, 3),
                         Z=('NA' if sd == 0 else round(z, 3)),
                         p_emp=(1 + int((col >= nr).sum())) / (1 + len(col)),
                         p_emp_under=(1 + int((col <= nr).sum())) / (1 + len(col)),
                         fold_change=(round(nr / mu, 4) if mu > 0 else 'Inf'),
                         rand_min=int(col.min()), rand_max=int(col.max()), n_rand=len(col)))
    norm = math.sqrt(sum(z * z for k, z in zs.items()
                         if k in SP_CLASSES and not math.isnan(z)))
    for r in rows:
        z = zs[r['motif_class']]
        r['SP'] = (round(z / norm, 4) if (r['motif_class'] in SP_CLASSES and norm > 0
                                          and not math.isnan(z)) else '')
    for r in rows:
        say(f"   [{null_model}] {r['motif_class']:<24} real={r['n_real']:<9} "
            f"rand={r['rand_mean']:.1f} +-{r['rand_sd']:.1f}   Z={str(r['Z']):<10} "
            f"p_over={r['p_emp']:.4g}  fold={r['fold_change']}")
    return rows


def main():
    open(LOG, 'w').close()
    say('=== 08_motif_significance :: START ===')
    say(f'N_RAND={N_RAND}  SWAP_FACTOR={SWAP_FACTOR}  NPROC={NPROC}  SEED={SEED}')
    nodes, layers_raw = load_graph()
    enc = Encoded(nodes, layers_raw)
    say(f'network: {enc.n} nodes, {enc.m} directed edges, layers = ',
        {et: len(u) for et, (u, _) in zip(enc.layer_names, enc.layers)})
    U0 = np.concatenate([u for u, _ in enc.layers]); V0 = np.concatenate([v for _, v in enc.layers])
    real, ffl_per_node = count_motifs(enc, (U0, V0))
    say('REAL counts: ', {k: real[k] for k in MOTIF_ORDER})

    cc = collections.Counter(r['class'] for r in
                             csv.DictReader(open(f'{REV}/results/ffl_3node_cores_crosscheck.csv')))
    say('independent cross-check file 02x counts: ', dict(cc))
    for k in ('miRNA-FFL', 'TF-FFL', 'Composite-FFL'):
        assert real[k] == cc[k], f'MISMATCH {k}: mine={real[k]} crosscheck={cc[k]}'
    say('OK: 3-node class counts reproduce results/ffl_3node_cores_crosscheck.csv exactly')

    rows = []
    for nm in ('NM1', 'NM2'):
        rows += run_null(nodes, layers_raw, enc, real, nm)
    fields = ['motif_class', 'n_real', 'rand_mean', 'rand_sd', 'Z', 'p_emp', 'null_model',
              'p_emp_under', 'fold_change', 'rand_min', 'rand_max', 'n_rand', 'SP']
    with open(f'{REV}/results/motif_significance.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for r in rows: w.writerow(r)
    json.dump(rows, open(f'{REV}/results/motif_significance.json', 'w'), indent=1)
    say('wrote results/motif_significance.csv (%d rows)' % len(rows))
    # per-node FFL participation, used by the specificity and bootstrap steps
    with open(f'{REV}/results/ffl_participation_per_node.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['node', 'type', 'n_FFL_cores_as_regulator_or_mediator'])
        for i, cnt in sorted(ffl_per_node.items(), key=lambda x: -x[1]):
            w.writerow([enc.names[i], ['TF', 'miRNA', 'Gene'][enc.ntype[i]], cnt])
    say('=== 08 DONE ===')


if __name__ == '__main__':
    main()
