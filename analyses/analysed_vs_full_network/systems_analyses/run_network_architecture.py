#!/usr/bin/env python3
"""Re-run the v3 architecture / controllability / information-flow scripts (Supplementary Note S4)
unchanged, with only two things swapped in scripts/05_network_architecture/netlib.py at run time:
  RES/FIG/LOG -> out_<mode>/ in this folder (the project is not written to)
  load_edges  -> mode 'full': unchanged (the 6,859-edge deposit, to check reproduction)
                 mode 'noexemplar': the deposit minus its 30 miRNA_miRNA edges (GeneMANIA exemplar)
usage: python3 run_network_architecture.py <full|noexemplar> <script.py> [...]"""
import sys, os, runpy
sys.dont_write_bytecode = True
V3 = '/path/to/six-node-ffl-breast-cancer/scripts/05_network_architecture'
HERE = os.path.dirname(os.path.abspath(__file__))
mode = sys.argv[1]; assert mode in ('full', 'noexemplar')
OUT = os.path.join(HERE, f"out_{mode}_hash{os.environ.get('PYTHONHASHSEED', 'random')}")
for d in (OUT, f'{OUT}/fig', f'{OUT}/log'): os.makedirs(d, exist_ok=True)
sys.path.insert(0, V3)
import netlib
netlib.RES, netlib.FIG, netlib.LOG = OUT, f'{OUT}/fig', f'{OUT}/log'
_load = netlib.load_edges
if mode == 'noexemplar':
    def load_edges():
        E = _load(); keep = [e for e in E if e['edge_type'] != 'miRNA_miRNA']
        assert len(E) - len(keep) == 30, len(E) - len(keep)
        return keep
    netlib.load_edges = load_edges
for s in sys.argv[2:]:
    print(f'\n######## {mode}: {s}', flush=True)
    sys.argv = [os.path.join(V3, s)]
    runpy.run_path(os.path.join(V3, s), run_name='__main__')
