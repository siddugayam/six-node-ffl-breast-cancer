#!/usr/bin/env python3
"""A2 of SETTINGS.md, item 2: the twelve Bhat networks under the two cluster variants, on the observed graph and on every
replicate graph of NULL-B, NULL-C and NULL-L (section O's program and settings; each replicate's arc list is written by
bhat_null.c through BHAT_DUMP and counted here with make_ffl_networks.py's enumeration).

Counts per graph and network (instances and node sets):
  base      unmodified (check: equal to the C program's counts on every replicate)
  noMIR17   (a) instances with no miR-17~92 member (hsa-miR-17, -18a, -19a, -20a, -19b, -92a)
  clusets   (b) distinct cluster-level node sets: each miRNA replaced by its 10 kb cluster (the observed cluster map)
  merged    (b) instances and node sets on the merged graph (each observed cluster one miRNA node, union of its members'
            arcs, self-loops dropped); five- and six-node merged counts are reported but not tested (zero on the observed
            graph by construction)
Statistics and decision rule as section O (p_upper = (#{x >= o} + 1) / (R + 1); over-represented if p_upper < 0.05 under
both NULL-C and NULL-L; explained by the three-node cores if under NULL-C only; depleted if p_lower < 0.05 under both).

usage: python3 cluster_nulls.py <analysis root> <network deposit folder> <original SIF folder> <nulls inputs folder>
                                <output folder> <workers>
"""
import os, sys, subprocess, collections, multiprocessing as mp
import numpy as np, pandas as pd
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bhat_common as bc

R, SEED, SPE = 1000, 20250908, 100
NULLS = [('NULL-B', 'B'), ('NULL-C', 'C'), ('NULL-L', 'L')]
M17 = {'hsa-miR-17', 'hsa-miR-18a', 'hsa-miR-19a', 'hsa-miR-20a', 'hsa-miR-19b', 'hsa-miR-92a'}
LT = {0: 'TF_miRNA', 1: 'miRNA_target', 2: 'TF_target', 3: 'TF_target', 4: 'gene_gene', 5: 'gene_gene', 6: 'miRNA_miRNA'}
CK = {'miRNA_FFL': 'mirFFL', 'TF_FFL': 'TFFFL', 'composite_FFL': 'comp'}
G = {}                                                  # per-process globals


def clusters(arcs, nodes):
    adj = collections.defaultdict(set)
    for u, v, l in arcs:
        if l == 6 and u != v: adj[u].add(v); adj[v].add(u)
    cm, seen = {}, set()
    for v in sorted(adj):
        if v in seen: continue
        comp, st = set(), [v]
        while st:
            x = st.pop()
            if x not in comp: comp.add(x); st.extend(adj[x] - comp)
        seen |= comp
        name = 'cluster:' + sorted(comp)[0]
        for x in comp: cm[x] = name
    return cm


def count(arcs):
    mfn, nodes, cm = G['mfn'], G['nodes'], G['cm']
    census = {}
    for u, v, l in arcs: census.setdefault((u, v), (LT[l],))
    bh = mfn.Bhat(nodes, census)
    out = {}
    for n in bc.SIZES:
        for k in bc.CLASSES:
            I = bh.enumerate(n, k); tag = f'{n}node_{k}'
            out[f'base|{tag}|inst'] = len(I); out[f'base|{tag}|sets'] = len({frozenset(i) for i in I})
            J = [i for i in I if not (set(i) & M17)]
            out[f'noMIR17|{tag}|inst'] = len(J); out[f'noMIR17|{tag}|sets'] = len({frozenset(i) for i in J})
            out[f'clusets|{tag}|sets'] = len({frozenset(cm.get(x, x) for x in i) for i in I})
    mc = {}
    for (u, v), t in census.items():
        a, b = cm.get(u, u), cm.get(v, v)
        if a != b: mc.setdefault((a, b), t)
    nm = dict(nodes); nm.update({c: 'miRNA' for c in set(cm.values())})
    bm = mfn.Bhat(nm, mc)
    for n in bc.SIZES:
        for k in bc.CLASSES:
            I = bm.enumerate(n, k); tag = f'{n}node_{k}'
            out[f'merged|{tag}|inst'] = len(I); out[f'merged|{tag}|sets'] = len({frozenset(i) for i in I})
    return out


def init(root, netdir, sifdir, cm):
    G['mfn'] = bc.builder(netdir); G['nodes'] = G['mfn'].load(root)[0]; G['cm'] = cm


def work(task):
    r, u, v, l = task
    names = G['names']
    return r, count([(names[a], names[b], int(c)) for a, b, c in zip(u, v, l)])


def main(root, netdir, sifdir, inp, out, workers):
    workers = int(workers)
    os.makedirs(os.path.join(out, 'runs'), exist_ok=True)
    log = []; say = lambda *a: (print(*a, flush=True), log.append(' '.join(map(str, a))))
    N = bc.networks(root, netdir, sifdir)
    say('networks rebuilt byte-identical to the deposit')
    GR, LB, NM = (os.path.join(inp, f) for f in ('graph_nolegacy_null.txt', 'graph_nolegacy_labels.txt', 'node_names.txt'))
    names = [l.strip() for l in open(NM) if l.strip()]
    L = open(GR).read().split('\n'); nv, ne, _ = map(int, L[0].split())
    labs = [int(x) for x in open(LB).read().split()]
    obs = [(names[int(a)], names[int(b)], lab) for (a, b, c), lab in zip((l.split() for l in L[2:2 + ne]), labs)]
    cm = clusters(obs, N['nodes'])
    cl = collections.defaultdict(list)
    for x, c in cm.items(): cl[c].append(x)
    say(f'10 kb clusters: {len(cl)} clusters, {len(cm)} miRNAs, sizes {sorted(collections.Counter(len(v) for v in cl.values()).items())}')
    pd.DataFrame(sorted((c, x) for c, xs in cl.items() for x in xs), columns=['cluster', 'miRNA']).to_csv(os.path.join(out, 'cluster_map.csv'), index=False)
    G['names'] = names; init(root, netdir, sifdir, cm)
    O = count(obs)
    for net in bc.NETWORKS:
        assert O[f'base|{net}|inst'] == len(N['inst'][net]) and O[f'base|{net}|sets'] == len({frozenset(i) for i in N['inst'][net]}), net
    say('CHECK observed: base counts equal the network builder\'s for all twelve networks: True')
    say('observed merged-graph five- and six-node instances: ' + ', '.join(f'{net} {O[f"merged|{net}|inst"]}' for net in bc.NETWORKS if net[0] in '56'))
    pd.DataFrame([dict(statistic=k, observed=v) for k, v in O.items()]).to_csv(os.path.join(out, 'observed.csv'), index=False)
    exe = os.path.join(out, 'bhat_null'); subprocess.run(['gcc', '-O2', '-o', exe, os.path.join(HERE, 'bhat_null.c')], check=True)
    rec = np.dtype([('r', '<i4'), ('u', '<i2', ne), ('v', '<i2', ne), ('l', '<i2', ne)])
    rows = {}
    for nm, m in NULLS:
        dump = os.path.join(out, f'dump_{m}.bin')
        with open(os.path.join(out, 'runs', f'null_{m}.tsv'), 'w') as fo, open(os.path.join(out, 'runs', f'null_{m}.err'), 'w') as fe:
            subprocess.run([exe, GR, LB, m, str(R), str(SEED), str(SPE), '0', '1', '0'], stdout=fo, stderr=fe, check=True, env=dict(os.environ, BHAT_DUMP=dump))
        D = np.fromfile(dump, dtype=rec); assert len(D) == R
        C = pd.read_csv(os.path.join(out, 'runs', f'null_{m}.tsv'), sep='\t'); C['rep'] = C.rep.astype(int); C = C.set_index('rep')
        with mp.Pool(workers, initializer=_winit, initargs=(root, netdir, sifdir, cm, names)) as pool:
            res = dict(pool.imap_unordered(work, ((int(d['r']), d['u'], d['v'], d['l']) for d in D), chunksize=4))
        bad = 0
        for r, x in res.items():
            for net in bc.NETWORKS:
                for w in ('inst', 'sets'):
                    if x[f'base|{net}|{w}'] != int(C.loc[r, f"bhat{net[0]}_{CK[net.split('_', 1)[1]]}_{w}"]): bad += 1
        say(f'CHECK {nm}: {len(res)} replicates; base counts equal to the C program on every replicate and network: {bad == 0} ({bad} differences)')
        assert bad == 0
        rows[m] = pd.DataFrame([dict(rep=r, **x) for r, x in sorted(res.items())])
        rows[m].to_csv(os.path.join(out, 'runs', f'variant_counts_{m}.csv'), index=False)
        os.remove(dump)
    os.remove(exe)
    T = []
    for nm, m in NULLS:
        X = rows[m]
        for stat in O:
            if stat.startswith('base'): continue
            x = X[stat].to_numpy(float); o = float(O[stat]); mu = x.mean(); sd = x.std(ddof=1)
            T.append(dict(statistic=stat, null=nm, observed=o, null_mean=mu, null_sd=sd, obs_over_mean=o / mu if mu > 0 else np.nan,
                          Z=(o - mu) / sd if sd > 0 else np.nan, p_upper=(np.sum(x >= o) + 1) / (R + 1), p_lower=(np.sum(x <= o) + 1) / (R + 1),
                          replicates=len(x)))
    T = pd.DataFrame(T); T.to_csv(os.path.join(out, 'variant_null_summary.csv'), index=False)
    dec = []
    for stat in O:
        if stat.startswith('base'): continue
        variant, net, w = stat.split('|')
        c = T[(T.statistic == stat) & (T.null == 'NULL-C')].iloc[0]; l = T[(T.statistic == stat) & (T.null == 'NULL-L')].iloc[0]
        b = T[(T.statistic == stat) & (T.null == 'NULL-B')].iloc[0]
        if variant == 'merged' and net[0] in '56': d = 'not tested (zero on the observed graph by construction)'
        elif c.p_upper < 0.05 and l.p_upper < 0.05: d = 'over-represented'
        elif c.p_upper < 0.05: d = 'explained by the three-node cores'
        elif c.p_lower < 0.05 and l.p_lower < 0.05: d = 'depleted'
        else: d = 'none of the three (not over-represented, not depleted)'
        dec.append(dict(variant=variant, network=net, count='instances' if w == 'inst' else 'node sets', observed=c.observed,
                        NULL_B_ratio=b.obs_over_mean, NULL_B_p_upper=b.p_upper, NULL_C_ratio=c.obs_over_mean, NULL_C_p_upper=c.p_upper,
                        NULL_C_p_lower=c.p_lower, NULL_L_ratio=l.obs_over_mean, NULL_L_p_upper=l.p_upper, NULL_L_p_lower=l.p_lower, decision=d))
    D = pd.DataFrame(dec); D.to_csv(os.path.join(out, 'variant_decision.csv'), index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    say('\n' + D.to_string(index=False))
    open(os.path.join(out, 'cluster_nulls.log'), 'w').write('\n'.join(log) + '\n')


def _winit(root, netdir, sifdir, cm, names):
    init(root, netdir, sifdir, cm); G['names'] = names


if __name__ == '__main__':
    main(*sys.argv[1:7])
