#!/usr/bin/env python3
"""Prepare this repository for running its scripts.

The scripts address the analysis root by the placeholder /path/to/revision. This script
1. replaces that placeholder with the absolute path of this repository in every script under scripts/ and analyses/,
2. links each network file of data/network/ into data/, where the scripts read it.
It changes nothing else, and `python3 set_root.py --undo` restores the placeholder (the links stay).

Other placeholders mark locations outside the repository that you set yourself (see the README):
/path/to/home/Desktop/DD/R_GPR/ (raw inputs of the first analysis round), /path/to/scratch (a temporary folder) and
/path/to/home/bin/bigBedToBed (the UCSC bigBedToBed program).

Usage: python3 set_root.py [--undo] [--dry-run]"""
import os, sys

PLACEHOLDER = "/path/to/revision"
ROOT = os.path.dirname(os.path.abspath(__file__))
CODE = (".py", ".R", ".sh", ".Rmd", ".c", ".cpp")
undo, dry = "--undo" in sys.argv, "--dry-run" in sys.argv
old, new = (ROOT, PLACEHOLDER) if undo else (PLACEHOLDER, ROOT)

changed = 0
for top in ("scripts", "analyses"):
    for d, _, files in os.walk(os.path.join(ROOT, top)):
        for f in files:
            if not f.endswith(CODE):
                continue
            p = os.path.join(d, f)
            with open(p, encoding="utf8", errors="surrogateescape", newline="") as fh:
                t = fh.read()
            if old in t:
                changed += 1
                if not dry:
                    with open(p, "w", encoding="utf8", errors="surrogateescape", newline="") as fh:
                        fh.write(t.replace(old, new))
print(f"{'would update' if dry else 'updated'} {changed} scripts: {old} -> {new}")

if not undo:
    net, linked = os.path.join(ROOT, "data", "network"), 0
    for f in sorted(os.listdir(net)):
        target, link = os.path.join(net, f), os.path.join(ROOT, "data", f)
        if f == "README.md" or os.path.lexists(link):
            continue
        linked += 1
        if not dry:
            try:
                os.symlink(os.path.join("network", f), link)
            except OSError:                      # file systems without symbolic links: copy instead
                import shutil
                shutil.copy2(target, link)
    print(f"{'would link' if dry else 'linked'} {linked} network files into data/")
