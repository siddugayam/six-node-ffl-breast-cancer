#!/bin/bash
# Sandboxes for the BioNet (01, 02, 04) and jActiveModules (01, 02, 05) chains of scripts/v7.
# Scripts are copied with the project root replaced by the sandbox root; inputs are copied
# (summary tables only); nothing is written to the project.
set -e
REV=/path/to/revision
HERE=$(cd "$(dirname "$0")" && pwd)
IN="data/canonical_nodes.tsv data/ffl_module_sets.rds results/v2/v2_DE_genes.csv results/v2/v2_DE_mirnas.csv
    results/BRCA_DEX_mirnas.csv results/v5/tables/Table4_prioritised_30.csv results/v5/tables/TableS4_node_prioritisation_full.csv
    cache/v5/farmer_signature_sets.rds cache/v5/farmer/farmer2009_stromal_signature_50.txt cache/celltype/caf_signature_scrna.txt"
for m in full nolegacy; do
  S=$HERE/sandbox_$m
  mkdir -p $S/scripts/v7 $S/results/v7 $S/figures/v7 $S/cache/v7/jam
  for f in $IN; do mkdir -p $S/$(dirname $f); cp -p $REV/$f $S/$f; done
  if [ $m = full ]; then cp -p $REV/data/canonical_edges.tsv $S/data/
  else awk -F'\t' 'NR==1 || $3!="miRNA_miRNA"' $REV/data/canonical_edges.tsv > $S/data/canonical_edges.tsv; fi
  echo "$m: $(($(wc -l < $S/data/canonical_edges.tsv)-1)) edges"
  cp -p $REV/scripts/v7/jam_engine.cpp $S/scripts/v7/
  for s in jam_common.R jam_01_search.R jam_02_characterise.R jam_05_circularity_check.R \
           bionet_01_bum_and_scan.R bionet_02_main.R bionet_04_stability_and_nulls.R; do
    sed "s#$REV#$S#g" $REV/scripts/v7/$s > $S/scripts/v7/$s
    [ -z "$(grep -o "$REV[^\"]*" $S/scripts/v7/$s | grep -v "legacy_audit/v7b/sandbox_")" ] || { echo "project path left in $s"; exit 1; }
  done
done
