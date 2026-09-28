#!/usr/bin/env python3
"""Aggregate ffl_census_composition partition outputs into one summary (JSON + readable text).

usage: aggregate_census.py out_prefix part1.txt [part2.txt ...]
Exhaustive runs: the partitions are disjoint sets of source vertices, so counts are summed.
RAND-ESU runs: pass ONE file per call (one seed); estimates are counts / prod_q.
"""
import sys, json, collections

CLS = ['TF->Gene', 'TF->TF', 'TF->miRNA', 'miRNA->Gene', 'miRNA->TF', 'Gene->Gene', 'miRNA-miRNA']
FLAGS = ['a_2TF_2miR_2Gene', 'b_nodetype_TFTF_GeneGene_TFmiRGene', 'b_paperclass_TFTF_GeneGene_TFmiRGene',
         'c_literal_TF1toTF2_TF2toT_TF2toM_TF1recipM', 'c_full_compositecore_plus_TF2']


def mask_names(m):
    return '+'.join(CLS[i] for i in range(7) if (m >> i) & 1)


def main():
    prefix, files = sys.argv[1], sys.argv[2:]
    K = None; prod = None; found = 0; visited = 0; parts = []
    comp = collections.Counter(); mask = collections.Counter(); flags = collections.Counter()
    exflag = collections.defaultdict(list); exmask = collections.defaultdict(list)
    for fn in files:
        for line in open(fn):
            w = line.rstrip('\n').split(' ')
            if w[0] == 'K': K = int(w[1])
            elif w[0] == 'prod_q': prod = float(w[1])
            elif w[0] == 'found': found += int(w[1])
            elif w[0] == 'visited': visited += int(w[1])
            elif w[0] == 'part': parts.append(int(w[1]))
            elif w[0] == 'comp': comp[(int(w[1]), int(w[2]), int(w[3]))] += int(w[4])
            elif w[0] == 'mask': mask[int(w[1])] += int(w[2])
            elif w[0] == 'flags': flags[int(w[1])] += int(w[2])
            elif w[0] == 'exflag': exflag[int(w[1])].append(' '.join(w[2:]))
            elif w[0] == 'exmask': exmask[int(w[1])].append(' '.join(w[2:]))
    scale = 1.0 / prod
    assert sum(comp.values()) == found == sum(mask.values()) == sum(flags.values())
    flag_tot = {FLAGS[f]: sum(v for k, v in flags.items() if (k >> f) & 1) for f in range(5)}
    ncls = collections.Counter()
    for m, v in mask.items():
        ncls[bin(m).count('1')] += v
    maxc = max(ncls)
    out = dict(
        K=K, n_partition_files=len(files), prod_q=prod, found=found, estimated_modules=found * scale,
        visited=visited,
        composition={f"TF{a}_miR{b}_G{c}": v for (a, b, c), v in sorted(comp.items())},
        flags_found={k: v for k, v in flag_tot.items()},
        flags_estimated={k: v * scale for k, v in flag_tot.items()},
        flags_pct={k: 100.0 * v / found for k, v in flag_tot.items()},
        flag_pattern_counts={'+'.join(FLAGS[f] for f in range(5) if (k >> f) & 1) or 'none': v
                             for k, v in sorted(flags.items())},
        n_classes_hist={str(k): v for k, v in sorted(ncls.items())},
        max_classes=maxc,
        class_sets_at_max={mask_names(m): v for m, v in sorted(mask.items(), key=lambda x: -x[1])
                           if bin(m).count('1') == maxc},
        class_sets_with_6_or_more={mask_names(m): v for m, v in sorted(mask.items(), key=lambda x: -x[1])
                                   if bin(m).count('1') >= 6},
        examples_by_flag={FLAGS[f]: exflag[f][:3] for f in range(5)},
        examples_by_class_set={mask_names(m): exmask[m][:3] for m in sorted(exmask)},
    )
    json.dump(out, open(prefix + '.json', 'w'), indent=1)
    with open(prefix + '.txt', 'w') as fh:
        fh.write(f"K={K}  files={len(files)}  prod_q={prod}  found={found:,}  estimate={found*scale:,.1f}  visited={visited:,}\n")
        fh.write("composition (TF, miRNA, Gene): found\n")
        for (a, b, c), v in sorted(comp.items(), key=lambda x: -x[1]):
            fh.write(f"  TF={a} miRNA={b} Gene={c}: {v:,}  ({100*v/found:.4f} %)\n")
        fh.write("flags: found, estimate, %\n")
        for k, v in flag_tot.items():
            fh.write(f"  {k}: {v:,}  est {v*scale:,.1f}  ({100*v/found:.4f} %)\n")
        fh.write("distinct paper edge classes per module: " +
                 ', '.join(f"{k}: {v:,}" for k, v in sorted(ncls.items())) + "\n")
        fh.write("class sets with >= 6 classes:\n")
        for m, v in sorted(mask.items(), key=lambda x: -x[1]):
            if bin(m).count('1') >= 6:
                fh.write(f"  [{bin(m).count('1')}] {mask_names(m)}: {v:,}\n")
        for f in range(5):
            fh.write(f"examples {FLAGS[f]}:\n")
            for e in exflag[f][:3]:
                fh.write(f"  {e}\n")
        for m in sorted(exmask):
            fh.write(f"examples class set {mask_names(m)}:\n")
            for e in exmask[m][:3]:
                fh.write(f"  {e}\n")
    print(open(prefix + '.txt').read())


if __name__ == '__main__':
    main()
