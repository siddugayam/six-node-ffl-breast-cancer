#!/usr/bin/env python3
"""
Shared loader for S2/S3/S4 (analyses/six_node_pattern): every per-set result of the 16,384-set Sobol design
(seed 20260908), stored and new, under the factorial module labels.  Read-only on the project.

  source                                                              family                module labels
  results/v3/dynamics_higher_order_persetset.csv.gz (2026-09-09)      I1_miRNA_FFL,         n3 core, n4 GG, n5 GG+MM,
                                                                      COMP_C2_toggle        n6 GG+MM+TT (n4d not used)
  results/v3/dynamics_higher_order_persetset_compI1.csv.gz (09-09)    COMP_I1_negfeedback   the same
  analyses/dynamics_checks/c3_persetset.csv.gz               COMP_C2_toggle        n4tf = TT (C3 control)
  analyses/six_node_pattern/S2/perset/<family>__<module>.csv.gz               all                   new modules (this run)
  analyses/six_node_pattern/S6/perset/COMP_I3_allneg_TT__<module>.csv.gz       COMP_I3_allneg_TT     S6 modules
"""
import os, glob
import numpy as np, pandas as pd

REV = "/path/to/revision"
HERE = os.path.dirname(os.path.abspath(__file__))
STORED_MAP = {"n3": "core", "n4": "GG", "n5": "GG+MM", "n6": "GG+MM+TT"}
FACTORIAL = ["core", "GG", "MM", "TT", "GG+MM", "GG+TT", "MM+TT", "GG+MM+TT"]
FAMILIES = ["COMP_C2_toggle", "COMP_I1_negfeedback", "I1_miRNA_FFL"]
BEH = {"memory": "is_memory", "ultrasensitivity": "is_ultrasensitive", "bistability": "is_bistable",
       "pulse generation": "is_pulse", "noise rejection": "is_noise_rejecting",
       "damped oscillation": "is_damped_osc", "sustained oscillation (flagged)": "is_sustained_osc"}


def load_all(include_verify=False):
    parts = []
    A = pd.read_csv(f"{REV}/results/v3/dynamics_higher_order_persetset.csv.gz")
    B = pd.read_csv(f"{REV}/results/v3/dynamics_higher_order_persetset_compI1.csv.gz")
    for X, src in ((A, "results/v3/dynamics_higher_order_persetset.csv.gz"),
                   (B, "results/v3/dynamics_higher_order_persetset_compI1.csv.gz")):
        X = X[X.module.isin(STORED_MAP)].copy()
        X["module"] = X.module.map(STORED_MAP); X["source"] = src; parts.append(X)
    C = pd.read_csv(f"{REV}/analyses/dynamics_checks/c3_persetset.csv.gz")
    C = C[C.module == "n4tf"].copy(); C["module"] = "TT"; C["family"] = "COMP_C2_toggle"
    C["source"] = "analyses/dynamics_checks/c3_persetset.csv.gz"; parts.append(C)
    for fn in sorted(glob.glob(f"{HERE}/perset/*.csv.gz")) + sorted(glob.glob(f"{HERE}/../S6/perset/*.csv.gz")):
        X = pd.read_csv(fn)
        if X.module.iloc[0].endswith("@verify") and not include_verify:
            continue
        X["source"] = "analyses/six_node_pattern/" + ("S6" if "/S6/" in fn else "S2") + "/perset/" + os.path.basename(fn); parts.append(X)
    D = pd.concat(parts, ignore_index=True)
    dup = D.duplicated(["family", "module", "param_set"])
    assert not dup.any(), D[dup][["family", "module"]].drop_duplicates()
    return D


def certified_flags(D):
    """per-set certified sustained oscillation: stored 03f result for COMP_C2_toggle GG+MM+TT, else the
    S2 certification file (s2_certify_limit_cycles.py) when present.  Returns dict (family, module) ->
    set of certified param_set; modules with flagged sets but no certification are absent."""
    out = {}
    c = pd.read_csv(f"{REV}/results/v3/dynamics_n6_limit_cycle_certified.csv")
    out[("COMP_C2_toggle", "GG+MM+TT")] = set(c.loc[c.certified_limit_cycle.astype(bool), "param_set"].astype(int))
    fn = f"{HERE}/s2_limit_cycle_certification.csv"
    if os.path.exists(fn):
        s = pd.read_csv(fn)
        for (f, m), g in s.groupby(["family", "module"]):
            if (f, m) == ("COMP_C2_toggle", "GG+MM+TT"):
                continue
            out[(f, m)] = set(g.loc[g.certified_limit_cycle.astype(bool), "param_set"].astype(int))
    return out
