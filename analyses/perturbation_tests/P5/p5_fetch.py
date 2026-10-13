#!/usr/bin/env python3
"""P5: fetch the FANTOM5 miRNA atlas tables (de Rie et al. 2017; the atlas viewer's human data package)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _dl
B = 'https://fantom.gsc.riken.jp/5/suppl/De_Rie_et_al_2017/vis_viewer/data/human/'
LIC = 'CC BY 4.0 (licence link on https://fantom.gsc.riken.jp/5/suppl/De_Rie_et_al_2017/)'
for f in ('datapackage.json', 'human.srna.cpm.txt', 'human.srna.samples.tsv', 'human.srna.cellontology.tsv', 'human.mirna.cellontology.tsv'):
    print(_dl.fetch(B + f, 'P5/fantom5', f, 'P5', 'FANTOM5 miRNA atlas (de Rie et al. 2017)', f, LIC))
