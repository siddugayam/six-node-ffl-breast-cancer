"""
ffl_graph.py -- construction of the canonical DIRECTED, SIGNED, TYPED regulatory
graph G used for the n-node feed-forward-loop (FFL) census.

Sources (all under revision/data/):
  canonical_edges.tsv     canonical network (miRNA_target arcs are taken from here)
  layer_TF_target.tsv     TRRUST directed TF/regulator -> target, signed
  layer_TF_miRNA.tsv      TransmiR directed TF -> miRNA, signed
  layer_gene_gene.tsv     TRRUST directed tier + STRING undirected tier
  layer_miRNA_miRNA.tsv   miRBase polycistronic co-transcription (UNDIRECTED)

Edge representation
-------------------
Every edge is stored as one or two arcs in `arcs`:  (u,v) -> dict
    layer      : which source layer it came from
    iclass     : interaction class, one of
                 TF-TF, TF-gene, TF-miRNA, miRNA-gene, miRNA-TF,
                 gene-TF, gene-gene, miRNA-miRNA, assoc-STRING
    sign       : +1 / -1 / 0 (0 = curated but unsignable, or unknown)
    sign_src   : 'TRRUST' | 'TransmiR' | 'miRNA_repression' |
                 'assumed_activation' | 'assumed_cotranscription' | 'unsigned'
    undirected : True if the underlying biological evidence is undirected
                 (miRNA-miRNA co-transcription, STRING association).  Undirected
                 edges are stored as a reciprocal ARC PAIR.
"""
import csv, os, sys
from collections import defaultdict

DATA = "/path/to/revision/data"


def _rd(path, header=True):
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t", quotechar='"')
        rows = list(r)
    return rows[1:] if header else rows


def iclass_of(ntype_u, ntype_v):
    m = {"TF": "TF", "Gene": "gene", "miRNA": "miRNA"}
    return "%s-%s" % (m[ntype_u], m[ntype_v])


def build(include_string=False, mirna_mirna_threshold="10kb",
          include_author_mirna_mirna=False, verbose=True):
    """Return (nodes, ntype, arcs, meta)."""
    log = []

    # ---- nodes -----------------------------------------------------------
    ntype = {}
    for row in _rd(os.path.join(DATA, "canonical_nodes.tsv")):
        ntype[row[0]] = row[1]
    nodes = sorted(ntype)
    log.append("nodes: %d (TF %d, Gene %d, miRNA %d)" % (
        len(nodes),
        sum(1 for n in nodes if ntype[n] == "TF"),
        sum(1 for n in nodes if ntype[n] == "Gene"),
        sum(1 for n in nodes if ntype[n] == "miRNA")))

    arcs = {}

    def add(u, v, layer, sign, sign_src, undirected=False):
        if u not in ntype or v not in ntype or u == v:
            return False
        key = (u, v)
        if key in arcs:
            # keep the more informative sign / the directed evidence
            old = arcs[key]
            if old["sign"] == 0 and sign != 0:
                old["sign"], old["sign_src"] = sign, sign_src
            if old["undirected"] and not undirected:
                old["undirected"] = False
                old["layer"] = layer
            return False
        arcs[key] = dict(layer=layer, iclass=iclass_of(ntype[u], ntype[v]),
                         sign=sign, sign_src=sign_src, undirected=undirected)
        return True

    # ---- 1. TRRUST directed regulator -> target --------------------------
    n1 = 0
    for row in _rd(os.path.join(DATA, "layer_TF_target.tsv")):
        s, t, mode, sign = row[0], row[1], row[2], int(row[3])
        src = "TRRUST" if sign != 0 else "unsigned"
        n1 += add(s, t, "TF_target", sign, src)
    log.append("layer_TF_target.tsv -> %d arcs" % n1)

    # ---- 2. TransmiR directed TF -> miRNA --------------------------------
    n2 = 0
    for row in _rd(os.path.join(DATA, "layer_TF_miRNA.tsv")):
        s, t, mode, sign = row[0], row[1], row[2], int(row[3])
        src = "TransmiR" if sign != 0 else "unsigned"
        n2 += add(s, t, "TF_miRNA", sign, src)
    log.append("layer_TF_miRNA.tsv -> %d arcs" % n2)

    # ---- 3. gene-gene TRRUST directed tier -------------------------------
    n3 = 0
    n_string = 0
    for row in _rd(os.path.join(DATA, "layer_gene_gene.tsv")):
        s, t, ev, directed, score = row[0], row[1], row[2], row[3], row[4]
        if directed == "TRUE":
            mode = score
            sign = {"Activation": 1, "Repression": -1}.get(mode, 0)
            src = "TRRUST" if sign != 0 else "unsigned"
            n3 += add(s, t, "gene_gene_TRRUST", sign, src)
        else:
            if include_string:
                a = add(s, t, "STRING_assoc", 0, "unsigned", undirected=True)
                b = add(t, s, "STRING_assoc", 0, "unsigned", undirected=True)
                n_string += (a + b)
    log.append("layer_gene_gene.tsv TRRUST-directed -> %d new arcs" % n3)
    log.append("layer_gene_gene.tsv STRING-undirected -> %d arcs (include_string=%s)"
               % (n_string, include_string))

    # ---- 4. miRNA -> target (canonical, always repression) ---------------
    n4 = 0
    n_can_tf_t = n_can_tf_m = 0
    author_mm = []
    for row in _rd(os.path.join(DATA, "canonical_edges.tsv")):
        s, t, et, sign = row[0], row[1], row[2], int(row[3])
        if et == "miRNA_target":
            n4 += add(s, t, "miRNA_target", -1, "miRNA_repression")
        elif et == "TF_target":
            n_can_tf_t += add(s, t, "TF_target", 1, "assumed_activation")
        elif et == "TF_miRNA":
            n_can_tf_m += add(s, t, "TF_miRNA", 1, "assumed_activation")
        elif et == "miRNA_miRNA":
            author_mm.append((s, t))
        elif et == "gene_gene":
            pass  # COL1A1->COL3A1: no directed evidence; handled by layer file
    log.append("canonical miRNA_target -> %d arcs" % n4)
    log.append("canonical TF_target arcs NOT already in TRRUST layer -> %d" % n_can_tf_t)
    log.append("canonical TF_miRNA arcs NOT already in TransmiR layer -> %d" % n_can_tf_m)
    log.append("author miRNA_miRNA edges in canonical file: %d (used=%s)"
               % (len(author_mm), include_author_mirna_mirna))

    # ---- 5. miRNA-miRNA polycistronic co-transcription (UNDIRECTED) ------
    col = {"3kb": 6, "10kb": 5, "50kb": 7}[mirna_mirna_threshold]
    n5 = 0
    n_mm_pairs = 0
    for row in _rd(os.path.join(DATA, "layer_miRNA_miRNA.tsv")):
        if row[col] != "TRUE":
            continue
        a, b = row[0], row[1]
        n_mm_pairs += 1
        n5 += add(a, b, "miRNA_miRNA", 1, "assumed_cotranscription", undirected=True)
        n5 += add(b, a, "miRNA_miRNA", 1, "assumed_cotranscription", undirected=True)
    log.append("layer_miRNA_miRNA.tsv @%s -> %d pairs, %d arcs"
               % (mirna_mirna_threshold, n_mm_pairs, n5))

    if include_author_mirna_mirna:
        n5b = 0
        for a, b in author_mm:
            n5b += add(a, b, "miRNA_miRNA_author", 1, "assumed_cotranscription", True)
            n5b += add(b, a, "miRNA_miRNA_author", 1, "assumed_cotranscription", True)
        log.append("author miRNA_miRNA -> %d extra arcs" % n5b)

    meta = dict(log=log, n_arcs=len(arcs))
    if verbose:
        for l in log:
            print("[graph] " + l)
        print("[graph] TOTAL ARCS: %d" % len(arcs))
    return nodes, ntype, arcs, meta


def adjacency(nodes, arcs):
    out = defaultdict(set)
    inn = defaultdict(set)
    for (u, v) in arcs:
        out[u].add(v)
        inn[v].add(u)
    return out, inn


if __name__ == "__main__":
    nodes, ntype, arcs, meta = build()
    from collections import Counter
    print(Counter(a["iclass"] for a in arcs.values()))
    print(Counter(a["layer"] for a in arcs.values()))
    print(Counter(a["sign"] for a in arcs.values()))
