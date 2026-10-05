#!/usr/bin/env python3
"""Runs one analysis step and appends a line to RUN_LOG.tsv (next to SETTINGS.md):
analysis, command, start, end, seconds, exit status, and the md5 of SETTINGS.md at the start.
In the logged command, a path inside this folder is written relative to it, and any other path by its last component.

usage: python3 run_logged.py <analysis> <program> [arguments ...]"""
import os, sys, time, hashlib, subprocess
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG = os.path.join(HERE, 'RUN_LOG.tsv')


def show(a):
    if not os.path.isabs(a) or not os.path.exists(a): return a
    r = os.path.relpath(a, HERE)
    return r if not r.startswith('..') else '<' + os.path.basename(a.rstrip('/')) + '>'


def main(analysis, cmd):
    s5 = hashlib.md5(open(os.path.join(HERE, 'SETTINGS.md'), 'rb').read()).hexdigest()
    t0 = time.time(); start = time.strftime('%Y-%m-%d %H:%M:%S %z')
    st = subprocess.call(cmd)
    end = time.strftime('%Y-%m-%d %H:%M:%S %z')
    new = not os.path.exists(LOG)
    with open(LOG, 'a') as f:
        if new: f.write('analysis\tcommand\tstart\tend\tseconds\texit_status\tSETTINGS_md5_at_start\n')
        f.write('\t'.join([analysis, ' '.join(show(a) for a in cmd), start, end, f'{time.time() - t0:.0f}', str(st), s5]) + '\n')
    return st


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2:]))
