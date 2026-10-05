#!/usr/bin/env python3
"""Shared loader for the analyses of the twelve Bhat FFL networks (see SETTINGS.md).

networks(root, netdir, sifdir) rebuilds the twelve networks with the network deposit's make_ffl_networks.py, imported
read-only from netdir, into a temporary folder. It stops unless every file listed in netdir/MANIFEST.tsv is rebuilt
byte-identical (md5). It returns the builder, the node types, the analysed network, the census graph, the interaction
table, the instances of every network and the deposited edge and node tables.

  root    the analysis root (data/, results/)
  netdir  the unzipped network deposit (the twelve networks, README.md, MANIFEST.tsv, make_ffl_networks.py)
  sifdir  the original submission's SIF folder (read by the builder for its README)
"""
import os, sys, io, csv, hashlib, tempfile, importlib.util, contextlib
sys.dont_write_bytecode = True
CLASSES = ('miRNA_FFL', 'TF_FFL', 'composite_FFL')
SIZES = (3, 4, 5, 6)
NETWORKS = [f'{n}node_{k}' for n in SIZES for k in CLASSES]
UNDIRECTED = ('gene_gene', 'miRNA_miRNA')           # links; every other edge type is a directed arc


def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()


def builder(netdir):
    spec = importlib.util.spec_from_file_location('make_ffl_networks', os.path.join(netdir, 'make_ffl_networks.py'))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def edge_key(source, target, edge_type):
    """an undirected link is one unordered pair; a directed arc is an ordered pair"""
    return frozenset((source, target)) if edge_type in UNDIRECTED else (source, target)


def tsv(path): return list(csv.DictReader(open(path, newline=''), delimiter='\t'))


def networks(root, netdir, sifdir):
    mfn = builder(netdir)
    man = {r['path']: r['md5'] for r in tsv(os.path.join(netdir, 'MANIFEST.tsv'))}
    assert all(f'{n}{s}' in man for n in NETWORKS for s in ('.sif', '_edges.tsv', '_nodes.tsv'))
    with tempfile.TemporaryDirectory() as t:
        with contextlib.redirect_stdout(io.StringIO()):
            mfn.main(root, sifdir, t)
        bad = [p for p, h in man.items() if md5(os.path.join(t, p)) != h]
    assert not bad, ('rebuilt network files differ from the deposit', bad)
    nodes, analysed, census, table = mfn.load(root)
    bh = mfn.Bhat(nodes, census)
    inst = {f'{n}node_{k}': bh.enumerate(n, k) for n in SIZES for k in CLASSES}
    edges = {net: tsv(os.path.join(netdir, net + '_edges.tsv')) for net in NETWORKS}
    vtab = {net: tsv(os.path.join(netdir, net + '_nodes.tsv')) for net in NETWORKS}
    for net in NETWORKS:                              # the instances give the deposited node and edge sets
        n = int(net[0]); k = net.split('_', 1)[1]
        V = {v for i in inst[net] for v in i}
        E = {(a, b) for i in inst[net] for (a, b, u) in bh.edges_of(n, k, i)}
        assert V == {r['node'] for r in vtab[net]} and E == {(r['source'], r['target']) for r in edges[net]}, net
    return dict(mfn=mfn, nodes=nodes, analysed=analysed, census=census, table=table, bhat=bh, inst=inst,
                edges=edges, vtab=vtab, manifest=man)


if __name__ == '__main__':                            # python3 bhat_common.py check <root> <netdir> <sifdir>
    assert sys.argv[1] == 'check'
    networks(*sys.argv[2:5])
    print('networks rebuilt byte-identical to the deposit')
