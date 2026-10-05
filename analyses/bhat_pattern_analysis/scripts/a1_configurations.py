#!/usr/bin/env python3
"""A1 of SETTINGS.md: the sign configurations of the Bhat patterns to model, and the module list.

Per network (class x size), the compulsory arcs whose sign the model takes from the instances (network edge files'
sign; '?' = unsigned, compatible with either sign):
  three to five nodes  TF-FFL and composite: TF1 -> miR1 (s_TM), TF1 -> target (s_TY: G1 at three nodes, G2 at four and
                       five);  miRNA-FFL: s_TY only (s_TM = 0: no TF1 -> miR1 arc)
  six nodes            composite: s_TM, TF1 -> TF2 (s_T12), TF2 -> G2 (s_T2Y);  TF-FFL: s_TM, s_T2Y (s_T12 = 0);
                       miRNA-FFL: s_T12, s_T2Y (s_TM = 0);  s_TY = +1 (the paper's; no TF1 -> gene arc is compulsory)
Greedy cover: at each step the fully signed configuration compatible with the most instances not yet covered (ties:
+1 before -1, arcs in the order listed); stop at >= 90 % of the network's instances.
TF2 -> TF1: six-node TF-FFL and composite configurations are run with s_T21 = +1, -1 and 0 (A4 rule); miRNA-FFL: 0.
Comparison modules added: for every six-node configuration, the same class's core, GG and GG+MM modules with the same
s_TM and s_TY = +1 (McNemar pairs). Modules identical to the paper's stored modules are not rerun (the stored per-set
files are used, as S2 did; the nested check verifies the extended code reproduces them).
Time estimate: CPU-hours per 16,384-set module from analyses/six_node_pattern S2's log (core 1.0, GG 2.4, GG+MM 2.7,
GG+MM+TT 3.2), divided by the workers.
usage: python3 a1_configurations.py <analysis root> <network deposit folder> <original SIF folder> <output folder> <workers>
"""
import os, sys, itertools, collections
import pandas as pd
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bhat_common as bc

LAY = {3: 'core', 4: 'GG', 5: 'GG+MM', 6: 'GG+MM+TT'}
CPUH = {'core': 1.0, 'GG': 2.4, 'GG+MM': 2.7, 'GG+MM+TT': 3.2}
STORED = {('composite_FFL', -1, 1): 'COMP_C2_toggle', ('composite_FFL', 1, 1): 'COMP_I1_negfeedback', ('TF_FFL', 1, 1): 'I1_miRNA_FFL'}


def arcs_for(n, k):
    if n < 6: return ['s_TY'] if k == 'miRNA_FFL' else ['s_TM', 's_TY']
    return {'composite_FFL': ['s_TM', 's_T12', 's_T2Y'], 'TF_FFL': ['s_TM', 's_T2Y'], 'miRNA_FFL': ['s_T12', 's_T2Y']}[k]


def main(root, netdir, sifdir, out, workers):
    workers = int(workers); os.makedirs(out, exist_ok=True)
    N = bc.networks(root, netdir, sifdir)
    mfn, census, table = N['mfn'], N['census'], N['table']
    sg = lambda a, b: {'1': 1, '-1': -1, '': 0}[mfn.attributes((a, b), census, table)[1]] if (a, b) in census else None
    rows, mods, log = [], {}, []
    for net in bc.NETWORKS:
        n, k = int(net[0]), net.split('_', 1)[1]
        A = arcs_for(n, k)
        pats = []
        for inst in N['inst'][net]:
            r = dict(zip(mfn.ROLES[n], inst))
            v = {'s_TM': sg(r['T1'], r['M1']), 's_TY': sg(r['T1'], r['G1'] if n == 3 else r['G2']) if n < 6 else None,
                 's_T12': sg(r['T1'], r['T2']) if n == 6 else None, 's_T2Y': sg(r['T2'], r['G2']) if n == 6 else None}
            assert all(v[a] is not None for a in A), (net, inst)          # compulsory arcs are present
            pats.append(tuple(v[a] for a in A))
        left = collections.Counter(pats); tot = len(pats); chosen = []
        cands = list(itertools.product((1, -1), repeat=len(A)))           # +1 before -1, arcs in order
        while sum(c for p, c in left.items()) > 0.1 * tot:
            comp = lambda cfg: sum(c for p, c in left.items() if all(x == 0 or x == y for x, y in zip(p, cfg)))
            best = max(cands, key=lambda cfg: (comp(cfg), -cands.index(cfg)))
            got = comp(best); chosen.append((best, got))
            left = collections.Counter({p: c for p, c in left.items() if not all(x == 0 or x == y for x, y in zip(p, best))})
        cov = tot - sum(left.values())
        log.append(f'{net}: {tot} instances; {len(set(pats))} sign patterns; {len(chosen)} configurations cover {cov} ({100 * cov / tot:.1f} %)')
        print(log[-1], flush=True)
        for i, (cfg, got) in enumerate(chosen, 1):
            d = dict(zip(A, cfg)); sTM = d.get('s_TM', 0); sTY = d.get('s_TY', 1)
            rows.append(dict(network=net, config=i, **{a: d.get(a, '') for a in ('s_TM', 's_TY', 's_T12', 's_T2Y')},
                             instances_newly_covered=got, cumulative_pct=None))
            t21s = (1, -1, 0) if (n == 6 and k != 'miRNA_FFL') else (0,)
            for t21 in t21s:
                spec = dict(cls=k, size=n, s_TM=sTM, s_TY=sTY if n < 6 else 1, s_T12=d.get('s_T12', 0 if (n == 6 and k == 'TF_FFL') else 1),
                            s_T21=t21, s_T2Y=d.get('s_T2Y', 1))
                mods[tuple(sorted(spec.items()))] = spec
            if n == 6:                                                      # comparison modules
                for m in (3, 4, 5):
                    spec = dict(cls=k, size=m, s_TM=sTM, s_TY=1, s_T12=1, s_T21=0, s_T2Y=1)
                    mods[tuple(sorted(spec.items()))] = spec
        cum = 0
        for r_ in rows:
            if r_['network'] == net: cum += r_['instances_newly_covered']; r_['cumulative_pct'] = round(100 * cum / tot, 1)
    pd.DataFrame(rows).to_csv(os.path.join(out, 'a1_configurations.csv'), index=False)
    ML = []
    for spec in mods.values():
        k, n = spec['cls'], spec['size']; lay = LAY[n]
        stored = STORED.get((k, spec['s_TM'], spec['s_TY'])) if (spec['s_T21'] == 0 and spec['s_T12'] == 1 and spec['s_T2Y'] == 1) else None
        tag = {'composite_FFL': 'Bhat_comp', 'TF_FFL': 'Bhat_TF', 'miRNA_FFL': 'Bhat_miR'}[k]
        sig = lambda x: {1: '+', -1: '-', 0: '0'}[x]
        fam = f"{tag}_TM{sig(spec['s_TM'])}_TY{sig(spec['s_TY'])}"
        mod = lay + ('' if n < 6 else f"_T12{sig(spec['s_T12'])}_T21{sig(spec['s_T21'])}_T2Y{sig(spec['s_T2Y'])}")
        ML.append(dict(family=fam, module=mod, registry='', size_key=lay, mir_to_TF=str(k != 'TF_FFL'),
                       input_node='MIR1' if k == 'miRNA_FFL' else 'TF1', s_TM=spec['s_TM'], s_TY=spec['s_TY'], s_T12=spec['s_T12'],
                       s_T21=spec['s_T21'], s_T2Y=spec['s_T2Y'], s_T2M=1, stored_as=f'{stored} {lay}' if stored else '',
                       cpu_hours=0.0 if stored else CPUH[lay]))
    ML = pd.DataFrame(ML).sort_values(['size_key', 'family', 'module']).reset_index(drop=True)
    ML.to_csv(os.path.join(out, 'a1_modules.csv'), index=False)
    run = ML[ML.stored_as == '']
    tot_h = run.cpu_hours.sum()
    log.append(f'modules: {len(ML)} ({len(ML) - len(run)} available as the paper\'s stored modules, {len(run)} to run); '
               f'CPU-hours {tot_h:.0f}; with {workers} workers about {tot_h / workers:.1f} h, plus limit-cycle certification')
    print(log[-1])
    open(os.path.join(out, 'a1_configurations.log'), 'w').write('\n'.join(log) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:6])
