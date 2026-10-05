#!/usr/bin/env python3
"""A1b of SETTINGS.md: a1_run.py with the four layers of S2's factorial that A1 did not need. They are defined as
analyses/six_node_pattern S2 (s2_run_modules.py, factorial()) defines them, by layer switches only:
  MM    = core + miRNA2          (gene_gene none,   mir_mir True,  tf_tf False; 4 nodes)
  TT    = core + TF2             (gene_gene none,   mir_mir False, tf_tf True;  4 nodes)
  GG+TT = core + gene-gene + TF2 (gene_gene mutual, mir_mir False, tf_tf True;  5 nodes)
  MM+TT = core + miRNA2 + TF2    (gene_gene none,   mir_mir True,  tf_tf True;  5 nodes)
Nothing else changes: the model (dyn_models_bhat.py), the 16,384-set sample (seed 20260908; Kx21 and nx21 from seed
20260909), run_module() and the chunked runner are a1_run.py's. Other A1b scripts import this module for the layers.
usage: python3 a1b_run.py <analysis root> <module list csv> <output folder> [workers] [chunk]"""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a1_run as A

A.LAYERS.update({"MM": dict(gene_gene="none", mir_mir=True, tf_tf=False),
                 "TT": dict(gene_gene="none", mir_mir=False, tf_tf=True),
                 "GG+TT": dict(gene_gene="mutual", mir_mir=False, tf_tf=True),
                 "MM+TT": dict(gene_gene="none", mir_mir=True, tf_tf=True)})
A.NN.update({"MM": 4, "TT": 4, "GG+TT": 5, "MM+TT": 5})

if __name__ == "__main__":
    A.main(*sys.argv[1:6])
