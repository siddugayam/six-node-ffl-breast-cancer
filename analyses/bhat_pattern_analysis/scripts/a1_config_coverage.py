#!/usr/bin/env python3
"""A1 (added 2026-09-29; SETTINGS.md, Changes and decisions, 11:24): how many instances each chosen six-node sign
configuration covers. a1_configurations.csv gives the greedy cover's newly-covered count: instances compatible with
the configuration and with none chosen before it. That count is a lower bound on the instances compatible with the
configuration. This script reads the instances and signs exactly as a1_configurations.py does and gives, per chosen
configuration:
  compatible        instances compatible with it (an unsigned arc matches either sign)
  exact             instances whose compulsory arcs all carry exactly its signs (no unsigned arc)
  newly_covered     the greedy cover's count (a1_configurations.csv)
usage: python3 a1_config_coverage.py <analysis root> <network deposit folder> <original SIF folder> <dynamics folder>
Output: <dynamics folder>/a1_config_coverage.csv"""
import os, sys
import pandas as pd
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bhat_common as bc
from a1_configurations import arcs_for


def main(root, netdir, sifdir, dyn):
    N = bc.networks(root, netdir, sifdir)
    mfn, census, table = N['mfn'], N['census'], N['table']
    sg = lambda a, b: {'1': 1, '-1': -1, '': 0}[mfn.attributes((a, b), census, table)[1]] if (a, b) in census else None
    CF = pd.read_csv(os.path.join(dyn, 'a1_configurations.csv'), dtype=str).fillna('')
    rows = []
    for net in [n for n in bc.NETWORKS if n.startswith('6')]:
        k = net.split('_', 1)[1]; A = arcs_for(6, k)
        pats = []
        for inst in N['inst'][net]:
            r = dict(zip(mfn.ROLES[6], inst))
            v = {'s_TM': sg(r['T1'], r['M1']), 's_T12': sg(r['T1'], r['T2']), 's_T2Y': sg(r['T2'], r['G2'])}
            pats.append(tuple(v[a] for a in A))
        for c in CF[CF.network == net].to_dict('records'):
            cfg = tuple(int(c[a]) for a in A)
            comp = sum(all(x == 0 or x == y for x, y in zip(p, cfg)) for p in pats)
            exact = sum(p == cfg for p in pats)
            rows.append(dict(network=net, config=int(c['config']), **{a: c[a] for a in ('s_TM', 's_T12', 's_T2Y')}, instances=len(pats),
                             compatible=comp, exact=exact, newly_covered=int(c['instances_newly_covered'])))
    R = pd.DataFrame(rows); R.to_csv(os.path.join(dyn, 'a1_config_coverage.csv'), index=False); print(R.to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:5])
