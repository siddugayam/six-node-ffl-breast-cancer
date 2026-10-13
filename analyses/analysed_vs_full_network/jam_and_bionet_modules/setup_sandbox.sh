#!/bin/bash
# Sandboxes for the BioNet (01, 02, 04) and jActiveModules (01, 02, 05) chains of scripts/14_module_detection.
# Scripts are copied with the project root replaced by the sandbox root; inputs are copied
# (summary tables only); nothing is written to the project.
set -e
REV=/path/to/six-node-ffl-breast-cancer
HERE=$(cd "$(dirname "$0")" && pwd)
IN="data/canonical_nodes.tsv data/ffl_module_sets.rds results/expression_validation/DE_genes.csv results/expression_validation/DE_mirnas.csv
    results/BRCA_DEX_mirnas.csv results/prioritisation_and_stroma/tables/prioritised_30.csv results/prioritisation_and_stroma/tables/node_prioritisation_table.csv
    cache/prioritisation_and_stroma/farmer_signature_sets.rds cache/prioritisation_and_stroma/farmer/farmer2009_stromal_signature_50.txt cache/celltype/caf_signature_scrna.txt"
for m in full noexemplar; do
  S=$HERE/sandbox_$m
  mkdir -p $S/scripts/14_module_detection $S/results/module_detection $S/figures/module_detection $S/cache/module_detection/jam
  for f in $IN; do mkdir -p $S/$(dirname $f); cp -p $REV/$f $S/$f; done
  if [ $m = full ]; then cp -p $REV/data/canonical_edges.tsv $S/data/
  else awk -F'\t' 'NR==1 || $3!="miRNA_miRNA"' $REV/data/canonical_edges.tsv > $S/data/canonical_edges.tsv; fi
  echo "$m: $(($(wc -l < $S/data/canonical_edges.tsv)-1)) edges"
  cp -p $REV/scripts/14_module_detection/jam_engine.cpp $S/scripts/14_module_detection/
  for s in jam_common.R jam_01_search.R jam_02_characterise.R jam_05_circularity_check.R \
           bionet_01_bum_and_scan.R bionet_02_main.R bionet_04_stability_and_nulls.R; do
    sed "s#$REV#$S#g" $REV/scripts/14_module_detection/$s > $S/scripts/14_module_detection/$s
    [ -z "$(grep -o "$REV[^\"]*" $S/scripts/14_module_detection/$s | grep -v "jam_and_bionet_modules/sandbox_")" ] || { echo "project path left in $s"; exit 1; }
  done
done
