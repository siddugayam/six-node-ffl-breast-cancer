#!/usr/bin/env python3
"""A1b of SETTINGS.md: certification of the sustained oscillation flagged in the A1b modules. It runs a1_certify.py
unchanged, with a1b_run.py's extended layers: first the validation against the stored six-node certification (83 of
108, set by set), then every flagged set (1,000 tau from two initial conditions).
usage: python3 a1b_certify.py <analysis root> <module list csv> <perset folder> <output csv> [workers]"""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a1b_run                       # the extended layers (reads the analysis root from argv[1])
import a1_certify

if __name__ == '__main__':
    a1_certify.main(*sys.argv[1:6])
