#!/usr/bin/env python3
"""
C2: feedback-loop structure of the modelled topologies, derived from dyn_models.py.

03f_n6_limit_cycle_certification.py states loop structure in its docstring but prints none, so the
loops are derived here from the model itself: for each topology the Jacobian of the right-hand side
is evaluated at the operating point (input S = 2.0, state integrated for 50 tau from zero, as in
03_higher_order_sweep.py) for the first 1,024 sets of the same Sobol design (seed 20260908).  An
interaction j -> i exists when dF_i/dx_j is non-zero in any set; its sign is recorded as +, - or
mixed.  Self-terms (degradation) are excluded.  Every elementary cycle of the interaction graph is
listed with its sign (product of edge signs; 'mixed' if any edge is mixed).  Reads the project
read-only.
"""
import sys, os, itertools
sys.dont_write_bytecode = True
import numpy as np, pandas as pd

REV = "/path/to/revision"
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{REV}/scripts/v3")
import dyn_models as D, dyn_sample as SMP
import dataclasses

NS = 1024
P = {k: np.asarray(v, float)[:NS] for k, v in SMP.sample(16384, seed=20260908).items()}
NAMES = ["TF1", "miR1", "G1_mRNA", "G1_protein", "G2_mRNA", "G2_protein", "miR2", "TF2"]
TOL = 1e-9


def interactions(topo):
    rhs = D.make_rhs(topo, P); act = topo.active_states()
    Z = D.integrate(rhs, np.zeros((D.NSTATE, NS)), lambda t: 2.0, 50.0, 0.01, active=act)
    f0 = rhs(0.0, Z, 2.0)
    E = {}
    for j in act:
        Zp = Z.copy(); h = 1e-6 * np.maximum(np.abs(Z[j]), 1.0); Zp[j] = Zp[j] + h
        d = (rhs(0.0, Zp, 2.0) - f0) / h[None, :]
        for i in act:
            if i == j:
                continue
            col = d[i]
            pos = (col > TOL).mean(); neg = (col < -TOL).mean()
            if pos == 0 and neg == 0:
                continue
            sign = "+" if neg == 0 else "-" if pos == 0 else "mixed"
            E[(j, i)] = dict(sign=sign, frac_nonzero=float(pos + neg), frac_pos=float(pos), frac_neg=float(neg))
    return act, E


def cycles(act, E):
    out = []
    for L in range(2, len(act) + 1):
        for combo in itertools.permutations(act, L):
            if combo[0] != min(combo):
                continue
            ok = all((combo[k], combo[(k + 1) % L]) in E for k in range(L))
            if not ok:
                continue
            signs = [E[(combo[k], combo[(k + 1) % L])]["sign"] for k in range(L)]
            if "mixed" in signs:
                s = "mixed"
            else:
                s = "positive" if signs.count("-") % 2 == 0 else "negative"
            frac = min(E[(combo[k], combo[(k + 1) % L])]["frac_nonzero"] for k in range(L))
            out.append(dict(length=L, loop=" -> ".join(NAMES[c] for c in combo) + " -> " + NAMES[combo[0]],
                            edge_signs=" ".join(signs), loop_sign=s, min_frac_sets_edge_active=frac))
    return out


if __name__ == "__main__":
    fams = {"COMP_C2_toggle": D.higher_order_registry("AND", "C2_MIR", composite=True),
            "I1_miRNA_FFL": D.higher_order_registry("AND", "I1_MIR", composite=False)}
    reg_c = fams["COMP_C2_toggle"]
    extra = {"n4tf": dataclasses.replace(reg_c["n6"], name="COMP_C2_MIR_n4tf", n_nodes=4,
                                         gene_gene="none", mir_mir=False)}
    rows_e, rows_c = [], []
    for fam, reg in fams.items():
        mods = dict((m, reg[m]) for m in ("n3", "n4", "n5", "n6"))
        if fam == "COMP_C2_toggle":
            mods.update(extra)
        for m, topo in mods.items():
            act, E = interactions(topo)
            for (j, i), v in sorted(E.items()):
                rows_e.append(dict(family=fam, module=m, source=NAMES[j], target=NAMES[i], **v))
            cyc = cycles(act, E)
            for c in cyc:
                rows_c.append(dict(family=fam, module=m, **c))
            npos = sum(c["loop_sign"] == "positive" for c in cyc); nneg = sum(c["loop_sign"] == "negative" for c in cyc)
            nmix = sum(c["loop_sign"] == "mixed" for c in cyc)
            print(f"\n{fam} {m}: {len(cyc)} feedback loops ({npos} positive, {nneg} negative, {nmix} mixed-sign)")
            for c in cyc:
                print(f"   [{c['loop_sign']:8s}] len {c['length']}  {c['loop']}   ({c['edge_signs']}; edge active in >= {100*c['min_frac_sets_edge_active']:.0f}% of sets)")
    pd.DataFrame(rows_e).to_csv(f"{OUT}/c2_interaction_signs.csv", index=False)
    pd.DataFrame(rows_c).to_csv(f"{OUT}/c2_feedback_loops.csv", index=False)
