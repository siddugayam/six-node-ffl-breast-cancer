#!/usr/bin/env python3
"""Variant networks for the reruns of A2 and A3 (SETTINGS.md), written in the network deposit's file format so that
enrichment.R, concordance.R, survival.R and prioritisation.R can read them:
  <network>_nodes.tsv  node, type, n_ffl
  <network>_edges.tsv  source, target, edge_type, sign, sign_source, evidence_tier, layer, n_ffl
                       (compulsory edges of the kept instances; attributes as make_ffl_networks.py writes them)
  observed_bhat<n>_<class>_instances.tsv  the kept instances as node indices, in the order of the instance listing
                       they come from (for survival.R's first-listed choice)
  MANIFEST.tsv         path, md5, size_bytes
Variants:
  noMIR17    (A2a) section O's observed listings without the instances that contain a miR-17~92 member
  clusters   (A2b) section O's observed listings reduced to one instance per distinct cluster-level node set (the first
             listed); n_ffl = the number of distinct cluster-level instances containing the node (a miRNA counts those
             containing its cluster); no edge files (items 3 and 4 are identical to sections E and X by construction)
  physical   (A3a) the observed listings of the physical-STRING graph (evidence nulls)
  validated  (A3b) the observed listings of the validated-target graph (evidence nulls)
usage: python3 variant_networks.py <analysis root> <network deposit folder> <original SIF folder> <variant>
                                   <listing folder> <node names file> <output folder> [cluster map csv]
"""
import os, sys, csv, hashlib, collections
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bhat_common as bc

M17 = {'hsa-miR-17', 'hsa-miR-18a', 'hsa-miR-19a', 'hsa-miR-20a', 'hsa-miR-19b', 'hsa-miR-92a'}


def listing_name(net):
    n, k = net[0], net.split('_', 1)[1]
    return f"observed_bhat{n}_{'composite' if (n == '6' and k == 'composite_FFL') else k}_instances.tsv"


def main(root, netdir, sifdir, variant, lst, names_file, out, cmap=None):
    os.makedirs(out, exist_ok=True)
    N = bc.networks(root, netdir, sifdir)
    mfn, nodes, census, table = N['mfn'], N['nodes'], N['census'], N['table']
    names = [l.strip() for l in open(names_file) if l.strip()]; idx = {x: i for i, x in enumerate(names)}
    cm = {r['miRNA']: r['cluster'] for r in csv.DictReader(open(cmap))} if cmap else {}
    log = []
    for net in bc.NETWORKS:
        n, k = int(net[0]), net.split('_', 1)[1]
        rows = [tuple(names[int(x)] for x in l.split()) for l in open(os.path.join(lst, listing_name(net))) if l.strip()]
        if variant in ('noMIR17', 'clusters'):
            assert sorted(rows) == sorted(N['inst'][net]), net          # section O's listing = the deposit's instances
        if variant == 'noMIR17': keep = [i for i in rows if not (set(i) & M17)]
        elif variant == 'clusters':
            seen, keep = set(), []
            for i in rows:
                c = frozenset(cm.get(x, x) for x in i)
                if c not in seen: seen.add(c); keep.append(i)
        else: keep = rows
        with open(os.path.join(out, listing_name(net)), 'w') as f:
            for i in keep: f.write('\t'.join(str(idx[x]) for x in i) + '\n')
        V, E = collections.Counter(), collections.Counter()
        if variant == 'clusters':
            ci = {tuple(cm.get(x, x) for x in i) for i in rows}                # distinct cluster-level instances
            for t in ci:
                for x in set(t): V[x] += 1
            members = collections.defaultdict(set)
            for x, c in cm.items(): members[c].add(x)
            V = collections.Counter({m: v for x, v in V.items() for m in (members[x] if x in members else {x})})
        else:
            for i in keep:
                V.update(set(i)); E.update(set(mfn.Bhat.edges_of(n, k, i)))
        with open(os.path.join(out, net + '_nodes.tsv'), 'w') as f:
            f.write('node\ttype\tn_ffl\n')
            for v in sorted(V): f.write(f'{v}\t{nodes[v]}\t{V[v]}\n')
        if variant != 'clusters':
            er = []
            for (a, b, und), c in E.items():
                if und:
                    arcs = sorted((e for e in ((a, b), (b, a)) if e in census), key=lambda e: census[e][1] != 'canonical')
                    at = mfn.attributes(arcs[0], census, table)
                else: at = mfn.attributes((a, b), census, table)
                er.append((a, b) + at + (c,))
            er.sort(key=lambda r: (r[0], r[1]))
            with open(os.path.join(out, net + '_edges.tsv'), 'w') as f:
                f.write('source\ttarget\tedge_type\tsign\tsign_source\tevidence_tier\tlayer\tn_ffl\n')
                for r in er: f.write('\t'.join(map(str, r)) + '\n')
        log.append(f'{variant} {net}: instances {len(rows)} -> kept {len(keep)}; nodes {len(V)}; edges {len(E) if variant != "clusters" else "not written"}')
        print(log[-1], flush=True)
    md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
    with open(os.path.join(out, 'MANIFEST.tsv'), 'w') as f:
        f.write('path\tmd5\tsize_bytes\n')
        for fn in sorted(os.listdir(out)):
            if fn != 'MANIFEST.tsv' and os.path.isfile(os.path.join(out, fn)): f.write(f'{fn}\t{md5(os.path.join(out, fn))}\t{os.path.getsize(os.path.join(out, fn))}\n')
    open(os.path.join(out, 'variant.log'), 'w').write('\n'.join(log) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:9])
