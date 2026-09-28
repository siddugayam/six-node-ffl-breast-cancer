#!/usr/bin/env python3
"""Erratum check for B5 (six-node composite count sent on 2026-09-26, INBOX_2026-09-26b/b5/).
For each of the four census graphs: the B5 value as sent, the same rule with TF1/TF2 restricted to TF-typed
nodes, and the literal Fig. 1d rule (composite needs TF1<->TF2).  Uses build_bhat_networks.G, which reads
any graph exported to INBOX_2026-09-26b/graphs/."""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_bhat_networks import G

b5 = {o['graph']: o['classes']['composite_FFL (TF1<->M1)'] for o in
      json.load(open(os.path.join(HERE, '..', '..', 'INBOX_2026-09-26b', 'b5', 'B5_bhat2024_six_node_pattern.json')))}
rows = []
for tag in ('table3', 'nolegacy', 'table3_nostring', 'nolegacy_nostring'):
    g = G(tag)
    relaxed = g.enumerate(6, 'Composite_FFL', relaxed=True); literal = g.enumerate(6, 'Composite_FFL')
    r = dict(graph=tag, B5_as_sent_instances=b5[tag]['instances'], B5_as_sent_sets=b5[tag]['distinct_six_node_sets'],
             either_direction_TF_typed_instances=len(relaxed), either_direction_TF_typed_sets=len({frozenset(i) for i in relaxed}),
             fig1d_literal_instances=len(literal), fig1d_literal_sets=len({frozenset(i) for i in literal}))
    rows.append(r); print(r)
with open(os.path.join(HERE, 'b5_erratum_six_node_composite.csv'), 'w') as f:
    f.write(','.join(rows[0]) + '\n')
    for r in rows: f.write(','.join(str(v) for v in r.values()) + '\n')
