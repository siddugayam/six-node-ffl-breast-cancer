#!/usr/bin/env python3
"""
S2 / S3 / S6 analysis (analyses/six_node_pattern).  All per-set results through s2_common.load_all(); certified
sustained oscillation through s2_common.certified_flags() (03f rule; stored 03f result for
COMP_C2_toggle GG+MM+TT, s2_limit_cycle_certification.csv for every other module).

Per family, module and behaviour: % of the 16,384 sets, sets gained / lost relative to the family's
three-node core (same parameter set), net change, exact McNemar p against the core.
Behaviours: memory, ultrasensitivity, bistability, pulse generation, noise rejection, damped oscillation,
certified sustained oscillation (flagged-but-uncertified sets are reported separately).
Joint behaviours (same parameter set): memory AND pulse; memory AND damped oscillation; memory AND
certified sustained oscillation; memory AND either oscillation (damped OR certified sustained; the two are
mutually exclusive by construction: damped needs Re(lead) < 0, sustained needs Re(lead) > 1e-4).
For each: % of all sets; and, among the sets in which the core LACKS that joint behaviour, % that have it.
Six-node-specific: sets in which GG+MM+TT shows the behaviour and none of the seven smaller factorial
modules does (only for families with all eight modules; a family with modules not run, e.g. S6 beyond its time
limit, gets the McNemar tests and decision columns against the modules available, listed in smaller_modules_not_run).  McNemar (exact binomial, two-sided): GG+MM+TT against each smaller module, per joint
behaviour; b = six-node only, c = smaller only.
Decision rule (fixed in the analysis plan): the sentence "of the modelled modules, only the six-node composite module
combines [pulse generation and oscillation] with the memory gain" survives only if, in COMP_C2_toggle, the
six-node module's joint prevalence (memory AND pulse; memory AND either oscillation) exceeds that of every
smaller module with McNemar p < 0.05.  Otherwise the smallest module (fewest nodes) whose joint prevalence
equals or exceeds the six-node module's is named.  Applied to COMP_I1 as well (and reported for the others).
S3: the COMP_C2_toggle six-node sign variants (and S3c on TT) against the unmodified modules and the core.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from scipy.stats import binomtest
import s2_common as C

SIZE = {"core": 3, "GG": 4, "MM": 4, "TT": 4, "GG+MM": 5, "GG+TT": 5, "MM+TT": 5, "GG+MM+TT": 6}
D = C.load_all()
cert = C.certified_flags(D)
D["certified_sustained"] = 0
for (f, m), s in cert.items():
    msk = (D.family == f) & (D.module == m) & D.param_set.isin(s); D.loc[msk, "certified_sustained"] = 1
flag_uncert = {(f, m) for (f, m), g in D.groupby(["family", "module"]) if g.is_sustained_osc.sum() > 0 and (f, m) not in cert}
if flag_uncert: print("WARNING: flagged sustained oscillation without certification for", sorted(flag_uncert))
BEH = [("memory", "is_memory"), ("ultrasensitivity", "is_ultrasensitive"), ("bistability", "is_bistable"),
       ("pulse generation", "is_pulse"), ("noise rejection", "is_noise_rejecting"), ("damped oscillation", "is_damped_osc"),
       ("certified sustained oscillation", "certified_sustained"), ("flagged sustained oscillation (sweep)", "is_sustained_osc")]
M = {}
for (f, m), g in D.groupby(["family", "module"]):
    g = g.set_index("param_set").sort_index(); assert len(g) == 16384 and (g.index == np.arange(16384)).all(), (f, m, len(g))
    b = {nm: (g[c].astype(float).values > 0.5) for nm, c in BEH}
    b["memory AND pulse"] = b["memory"] & b["pulse generation"]
    b["memory AND damped oscillation"] = b["memory"] & b["damped oscillation"]
    b["memory AND certified sustained oscillation"] = b["memory"] & b["certified sustained oscillation"]
    b["memory AND either oscillation"] = b["memory"] & (b["damped oscillation"] | b["certified sustained oscillation"])
    assert not (b["damped oscillation"] & b["certified sustained oscillation"]).any()
    M[(f, m)] = b
JOINT = ["memory AND pulse", "memory AND damped oscillation", "memory AND certified sustained oscillation", "memory AND either oscillation"]
ALLB = [nm for nm, _ in BEH] + JOINT


def mcn(x, y):
    b = int((x & ~y).sum()); c = int((~x & y).sum())
    return b, c, (binomtest(b, b + c, 0.5).pvalue if b + c else 1.0)


rows = []; jrows = []
for (f, m), b in sorted(M.items()):
    core = M[(f, "core")]
    for nm in ALLB:
        x, c0 = b[nm], core[nm]
        gain = int((x & ~c0).sum()); loss = int((~x & c0).sum())
        rows.append(dict(family=f, module=m, nodes=SIZE.get(m.split("_S3")[0], ""), behaviour=nm, n_sets=int(x.sum()), pct=100 * x.mean(),
                         pct_core=100 * c0.mean(), n_gained_vs_core=gain, n_lost_vs_core=loss, pct_gained=100 * gain / 16384,
                         pct_lost=100 * loss / 16384, net_pct=100 * (gain - loss) / 16384,
                         mcnemar_p_vs_core=(binomtest(gain, gain + loss, 0.5).pvalue if gain + loss else 1.0)))
        if nm in JOINT:
            lack = ~c0
            jrows.append(dict(family=f, module=m, joint_behaviour=nm, n_sets=int(x.sum()), pct_all_sets=100 * x.mean(),
                              n_sets_core_lacks=int(lack.sum()), n_with_joint_where_core_lacks=int((x & lack).sum()),
                              pct_where_core_lacks=100 * (x & lack).sum() / max(1, lack.sum())))
R = pd.DataFrame(rows); J = pd.DataFrame(jrows)
R.to_csv(f"{HERE}/s2_prevalence_gain_loss.csv", index=False); J.to_csv(f"{HERE}/s2_joint_behaviours.csv", index=False)

# six-node-specific and McNemar six vs smaller
srows, mrows, drows = [], [], []
for f in sorted({f for f, _ in M}):
    if (f, "GG+MM+TT") not in M: print("no six-node module, skipped for six-node tests:", f); continue
    missing = [m for m in SIZE if (f, m) not in M]
    if missing: print("family incomplete: six-node tests against the available smaller modules only;", f, "missing", missing)
    six = M[(f, "GG+MM+TT")]; smaller = [m for m in SIZE if m != "GG+MM+TT" and (f, m) in M]
    for nm in ALLB if not missing else []:          # six-node-specific needs all seven smaller modules
        anys = np.zeros(16384, bool)
        for m in smaller: anys |= M[(f, m)][nm]
        srows.append(dict(family=f, behaviour=nm, n_six_node=int(six[nm].sum()), n_six_node_only=int((six[nm] & ~anys).sum()),
                          pct_six_node_only_of_all_sets=100 * (six[nm] & ~anys).mean()))
    for nm in JOINT:
        for m in smaller:
            b_, c_, p = mcn(six[nm], M[(f, m)][nm])
            mrows.append(dict(family=f, joint_behaviour=nm, smaller_module=m, nodes=SIZE[m], pct_six=100 * six[nm].mean(),
                              pct_smaller=100 * M[(f, m)][nm].mean(), b_six_only=b_, c_smaller_only=c_, mcnemar_p=p))
    for nm in ("memory AND pulse", "memory AND either oscillation", "memory AND damped oscillation", "memory AND certified sustained oscillation"):
        ps = 100 * six[nm].mean()
        beats = all(100 * M[(f, m)][nm].mean() < ps and mcn(six[nm], M[(f, m)][nm])[2] < 0.05 for m in smaller)
        ge = [m for m in smaller if 100 * M[(f, m)][nm].mean() >= ps]
        smin = min([SIZE[m] for m in ge]) if ge else None
        drows.append(dict(family=f, joint_behaviour=nm, pct_six_node=ps, six_node_exceeds_every_smaller_module_p05=beats,
                          smallest_modules_with_equal_or_higher_prevalence=';'.join(m for m in ge if SIZE[m] == smin) if ge else '',
                          their_pct=';'.join(f"{100 * M[(f, m)][nm].mean():.3f}" for m in ge if SIZE[m] == smin) if ge else '',
                          smaller_modules_not_significantly_below=';'.join(
                              f"{m} ({100 * M[(f, m)][nm].mean():.3f} %, p = {mcn(six[nm], M[(f, m)][nm])[2]:.3g})" for m in smaller
                              if not (100 * M[(f, m)][nm].mean() < ps and mcn(six[nm], M[(f, m)][nm])[2] < 0.05)),
                          smaller_modules_not_run=';'.join(missing)))
pd.DataFrame(srows).to_csv(f"{HERE}/s2_six_node_specific.csv", index=False)
pd.DataFrame(mrows).to_csv(f"{HERE}/s2_mcnemar_six_vs_smaller.csv", index=False)
DR = pd.DataFrame(drows); DR.to_csv(f"{HERE}/s2_decision.csv", index=False)
# S3: each sign variant against its unmodified module (McNemar, every behaviour and joint behaviour)
s3rows = []
for var, ref in (("GG+MM+TT_S3a_TF2-|miR", "GG+MM+TT"), ("GG+MM+TT_S3b_TF2-|G1", "GG+MM+TT"), ("GG+MM+TT_S3c_TF1-|TF2", "GG+MM+TT"),
                 ("TT_S3c_TF1-|TF2", "TT")):
    if ("COMP_C2_toggle", var) not in M: continue
    for nm in ALLB:
        x, y = M[("COMP_C2_toggle", var)][nm], M[("COMP_C2_toggle", ref)][nm]
        b_, c_, pv = mcn(x, y)
        s3rows.append(dict(variant=var, reference=ref, behaviour=nm, pct_variant=100 * x.mean(), pct_reference=100 * y.mean(),
                           n_variant_only=b_, n_reference_only=c_, mcnemar_p=pv))
pd.DataFrame(s3rows).to_csv(f"{HERE}/s3_variant_vs_unmodified.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20); pd.set_option("display.max_rows", 500)
print(DR.to_string(index=False))
for f in ("COMP_C2_toggle", "COMP_I1_negfeedback"):
    d = DR[(DR.family == f) & DR.joint_behaviour.isin(["memory AND pulse", "memory AND either oscillation"])]
    if len(d) == 2:
        ok = bool(d.six_node_exceeds_every_smaller_module_p05.all())
        print(f"DECISION {f}: sentence survives = {ok}; " + "; ".join(
            f"{r.joint_behaviour}: six-node {r.pct_six_node:.3f} %, smallest module(s) at or above: {r.smallest_modules_with_equal_or_higher_prevalence or 'none'} ({r.their_pct})"
            for r in d.itertuples()))
