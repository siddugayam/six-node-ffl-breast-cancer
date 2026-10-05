#!/usr/bin/env python3
"""Section O of SETTINGS.md: over-representation of the twelve Bhat networks under NULL-B, NULL-C and NULL-L.

Steps:
  1. Build bhat_null.c: the S1 null program with the Bhat three- to five-node counters added.
  2. Check 1: the input graph (S1's inputs) equals the census graph of make_ffl_networks.py, arc by arc and layer by
     layer.
  3. Observed graph: the counts of all twelve networks, and check 2: the instance lists equal the network builder's.
  4. Nulls B, C, L, LT, LG and LM: 1,000 replicates, seed 20250908, 100 swap attempts per rewired item, run in parallel.
  5. Check 3: in every replicate, S1's columns equal S1's stored runs. Check 4: the new three-node composite count
     equals S1's bhat3_comp.
  6. Statistics and the decision rule, as S1's summariser.

usage: python3 nulls.py <analysis root> <network deposit folder> <original SIF folder> <S1 stored runs folder> <output folder>
"""
import os, sys, csv, subprocess, collections
import numpy as np, pandas as pd
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bhat_common as bc

R, SEED, SPE = 1000, 20250908, 100
NULLS = [('NULL-B', 'B'), ('NULL-C', 'C'), ('NULL-L', 'L'), ('NULL-L TRRUST layer only', 'LT'),
         ('NULL-L STRING layer only', 'LG'), ('NULL-L 10 kb layer only', 'LM')]
S1COLS = ['total3', 'comp3', 'tf3', 'mir3', 'other3', 'recip', 'TTrecip', 'bhat3_comp', 'bhat6_mirFFL_inst', 'bhat6_mirFFL_sets',
          'bhat6_TFFFL_inst', 'bhat6_TFFFL_sets', 'bhat6_comp_inst', 'bhat6_comp_sets', 'model6_full_inst', 'model6_full_sets',
          'model6_def_inst', 'model6_def_sets']
CK = {'miRNA_FFL': 'mirFFL', 'TF_FFL': 'TFFFL', 'composite_FFL': 'comp'}
LABEL = {0: 'TF_miRNA', 1: 'miRNA_target', 2: 'TF_target', 3: 'TF_target', 4: 'gene_gene', 5: 'gene_gene', 6: 'miRNA_miRNA'}
LAYER = {0: 'canonical', 1: 'canonical', 2: 'canonical', 3: 'layer_TF_target', 4: 'canonical', 5: 'layer_gene_gene',
         6: 'layer_miRNA_miRNA'}


def col(net, what): return f"bhat{net[0]}_{CK[net.split('_', 1)[1]]}_{what}"


def main(root, netdir, sifdir, s1runs, out):
    os.makedirs(os.path.join(out, 'runs'), exist_ok=True); os.makedirs(os.path.join(out, 'obs'), exist_ok=True)
    inp = os.path.join(out, 'inputs'); G, L, NM = (os.path.join(inp, f) for f in
                                                  ('graph_nolegacy_null.txt', 'graph_nolegacy_labels.txt', 'node_names.txt'))
    log = []; say = lambda *a: (print(*a, flush=True), log.append(' '.join(map(str, a))))
    exe = os.path.join(out, 'bhat_null')
    subprocess.run(['gcc', '-O2', '-o', exe, os.path.join(HERE, 'bhat_null.c')], check=True)

    # check 1: the input graph is the census graph of the network builder
    N = bc.networks(root, netdir, sifdir)
    names = [l.strip() for l in open(NM) if l.strip()]
    lines = open(G).read().split('\n'); nv, ne, _ = map(int, lines[0].split())
    types = list(map(int, lines[1].split())); labs = [int(x) for x in open(L).read().split()]
    arcs = {}
    for l, lab in zip(lines[2:2 + ne], labs):
        u, v, _ = map(int, l.split()); arcs[(names[u], names[v])] = lab
    tmap = {'miRNA': 0, 'TF': 1, 'Gene': 2}
    ok1 = (nv == len(names) == 587 and ne == len(labs) == len(arcs) == 9226 and set(arcs) == set(N['census'])
           and all(tmap[N['nodes'][names[i]]] == types[i] for i in range(nv))
           and all(N['census'][e][0] == LABEL[b] and N['census'][e][1] == LAYER[b] for e, b in arcs.items()))
    say(f'CHECK 1 input graph = census graph of make_ffl_networks.py (9,226 arcs, types and layers): {ok1}')
    assert ok1

    # observed counts and check 2
    obs_tsv = os.path.join(out, 'obs', 'observed.tsv'); lp = os.path.join(out, 'obs', 'observed')
    with open(obs_tsv, 'w') as f:
        subprocess.run([exe, G, L, 'OBS', '0', str(SEED), str(SPE), '0', '1', '0', lp], stdout=f, check=True)
    O = pd.read_csv(obs_tsv, sep='\t').iloc[0]
    ok2 = True
    for net in bc.NETWORKS:
        n, k = int(net[0]), net.split('_', 1)[1]
        fn = f"{lp}_bhat{n}_{'composite' if (n == 6 and k == 'composite_FFL') else k}_instances.tsv"
        mine = [tuple(names[int(x)] for x in l.split()) for l in open(fn) if l.strip()]
        py = N['inst'][net]
        same = (sorted(mine) == sorted(py) and int(O[col(net, 'inst')]) == len(py)
                and int(O[col(net, 'sets')]) == len({frozenset(i) for i in py}))
        ok2 &= same
        say(f'CHECK 2 {net:20s} observed instances {int(O[col(net, "inst")]):6d} sets {int(O[col(net, "sets")]):6d}; '
            f'instance list equal to the builder\'s: {same}')
    assert ok2
    for fn in os.listdir(os.path.join(out, 'obs')):                  # keep only the counts and the Bhat lists
        if 'census6c' in fn or 'model6' in fn: os.remove(os.path.join(out, 'obs', fn))

    # the nulls, in parallel
    procs = []
    for nm, m in NULLS:
        fo = open(os.path.join(out, 'runs', f'null_{m}.tsv'), 'w'); fe = open(os.path.join(out, 'runs', f'null_{m}.err'), 'w')
        procs.append((m, subprocess.Popen([exe, G, L, m, str(R), str(SEED), str(SPE), '0', '1', '0'], stdout=fo, stderr=fe), fo, fe))
    for m, p, fo, fe in procs:
        assert p.wait() == 0, m
        fo.close(); fe.close()
    runs = {}
    for nm, m in NULLS:
        X = pd.read_csv(os.path.join(out, 'runs', f'null_{m}.tsv'), sep='\t'); X['rep'] = X.rep.astype(int)
        runs[m] = X.sort_values('rep').reset_index(drop=True)
        S = pd.read_csv(os.path.join(s1runs, f'null_{m}.tsv'), sep='\t')
        S = S[S.rep.astype(str) != 'obs'].copy(); S['rep'] = S.rep.astype(int)
        M = runs[m].merge(S, on='rep', suffixes=('', '_s1'))
        ok3 = len(runs[m]) == R and len(M) == R and all((M[c] == M[c + '_s1']).all() for c in S1COLS)
        ok4 = bool((runs[m].bhat3_comp_inst == runs[m].bhat3_comp).all())
        say(f'CHECK 3 {nm}: {len(runs[m])} replicates; S1 columns equal to the stored S1 run in every replicate: {ok3}')
        say(f'CHECK 4 {nm}: new three-node composite count equal to S1 bhat3_comp in every replicate: {ok4}')
        assert ok3 and ok4

    # statistics and decision rule (S1's summariser)
    rows = []
    for nm, m in NULLS:
        X = runs[m]
        for net in bc.NETWORKS:
            for what, lab in (('inst', 'instances'), ('sets', 'node sets')):
                c = col(net, what); x = X[c].to_numpy(float); o = float(O[c]); mu = x.mean(); sd = x.std(ddof=1)
                rows.append(dict(network=net, count=lab, column=c, null=nm, observed=o, null_mean=mu, null_sd=sd,
                                 obs_over_mean=o / mu if mu > 0 else np.nan, Z=(o - mu) / sd if sd > 0 else np.nan,
                                 p_upper=(np.sum(x >= o) + 1) / (R + 1), p_lower=(np.sum(x <= o) + 1) / (R + 1),
                                 replicates=len(x), null_min=x.min(), null_max=x.max()))
    T = pd.DataFrame(rows); T.to_csv(os.path.join(out, 'null_summary.csv'), index=False)
    dec = []
    for net in bc.NETWORKS:
        for what, lab in (('inst', 'instances'), ('sets', 'node sets')):
            c = T[(T.column == col(net, what)) & (T.null == 'NULL-C')].iloc[0]; l = T[(T.column == col(net, what)) & (T.null == 'NULL-L')].iloc[0]
            if c.p_upper < 0.05 and l.p_upper < 0.05: d = 'over-represented'
            elif c.p_upper < 0.05: d = 'explained by the three-node cores'
            elif c.p_lower < 0.05 and l.p_lower < 0.05: d = 'depleted'
            else: d = 'none of the three (not over-represented, not depleted)'
            b = T[(T.column == col(net, what)) & (T.null == 'NULL-B')].iloc[0]
            dec.append(dict(network=net, count=lab, observed=c.observed, NULL_B_ratio=b.obs_over_mean, NULL_B_p_upper=b.p_upper,
                            NULL_C_ratio=c.obs_over_mean, NULL_C_p_upper=c.p_upper, NULL_C_p_lower=c.p_lower,
                            NULL_L_ratio=l.obs_over_mean, NULL_L_p_upper=l.p_upper, NULL_L_p_lower=l.p_lower,
                            NULL_L_null_sd=l.null_sd, decision=d))
    D = pd.DataFrame(dec); D.to_csv(os.path.join(out, 'decision.csv'), index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    say('\n' + D.to_string(index=False))
    for m in ('LT', 'LG', 'LM'):
        say(f'\nNULL-{m} observed / mean:', ', '.join(f"{r.network} {r.obs_over_mean:.3f} (p_upper {r.p_upper:.4f})"
                                                  for r in T[(T.null == dict((b, a) for a, b in NULLS)[m]) & (T['count'] == 'instances')].itertuples()))
    open(os.path.join(out, 'nulls.log'), 'w').write('\n'.join(log) + '\n')
    os.remove(exe)


if __name__ == '__main__':
    main(*sys.argv[1:6])
