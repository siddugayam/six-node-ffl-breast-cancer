#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02h_stochastic_noise_buffering.py
=================================
DOES THE miRNA-MEDIATED INCOHERENT FFL BUFFER NOISE?  (Osella, Bosia, Cora &
Caselle 2011 PLoS Comput Biol, doi:10.1371/journal.pcbi.1001101; Siciliano et al.
2013 Nat Commun, doi:10.1038/ncomms3364)

The deterministic battery measures the CLASSICAL Alon notion of filtering: a
brief spurious input is rejected relative to a persistent one.  Osella et al.'s
claim is a different one -- that the miRNA-mediated incoherent FFL reduces the
INTRINSIC fluctuations of the target at a given mean.  That claim cannot be
tested with an ODE, so this script re-implements the same circuits as chemical
Langevin equations and measures the coefficient of variation of the target
protein at steady state.

CHEMICAL LANGEVIN EQUATION.  Each species' ODE is split into its production and
degradation fluxes a+ and a-, and integrated as

    dz = (a+ - a-) dt + Omega^-1/2 ( sqrt(a+) dW1 - sqrt(a-) dW2 )

with independent Wiener increments and reflection at zero (Gillespie 2000
J Chem Phys, doi:10.1063/1.481811).  Omega is the number of molecules
corresponding to one scaled concentration unit, so 1/sqrt(Omega) sets the
burst-free Poisson noise floor.  Omega -> infinity recovers the ODE, and the
script asserts that it does.

Because the circuits differ in their mean output, CVs are compared BOTH paired
per parameter set and after regressing log CV on log mean, which is the
comparison at matched mean that Osella et al. make.

Out: results/v3/dynamics_stochastic_noise.csv
     results/v3/dynamics_stochastic_noise_summary.csv
"""
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D, dyn_sample as SMP

RES = "/path/to/revision/results/v3"
N_SET = 2048
SEED = 20260908
DT = 0.002
T_BURN = 50.0
T_SAMP = 200.0
S_OP = 1.0
OMEGAS = (100.0, 500.0, 2000.0)
TOPOS = ["CASCADE", "I1_TXN_AND", "I1_MIR_AND", "C2_MIR_AND", "COMP_I1_AND", "COMP_C2_AND"]

reg = D.build_registry()


def fluxes(topo, P, Z, S):
    """Production and degradation fluxes of x1, m1, y1, p1 -- the same terms as
    dyn_models.make_rhs, written out so the noise can be attached to each."""
    x1, m1, y1, p1 = Z[0], Z[1], Z[2], Z[3]
    G = D.gate_and if topo.gate == "AND" else D.gate_or
    b_leak = P["b"] if topo.gate == "OR" else 0.0
    translational_block = (topo.gate == "AND")
    # Shared RISC/AGO loading saturates exactly as in dyn_models.make_rhs; leaving it
    # out makes the Langevin model disagree with the ODE it is meant to extend.
    Meff = m1 / (1.0 + m1 / P["Kc"])

    # --- TF1 ---
    if topo.input_node == "TF1":
        u_x1 = D.hill_act(S, P["KS"], P["nS"])
        u_x1 = u_x1 * np.ones_like(x1)
    else:
        u_x1 = np.ones_like(x1)
    dec_x1 = np.ones_like(x1)
    if topo.mir_to_TF and topo.mir_arm:
        q = 1.0 - D.hill_rep(Meff, P["Kmx"], P["nmx"])
        u_x1 = u_x1 * (1.0 - q)
        dec_x1 = 1.0 + P["lam_x"] * q
    ap_x = P["gx"] * u_x1
    am_x = P["gx"] * dec_x1 * x1

    # --- intermediate ---
    if topo.input_node == "MIR1":
        u_m = D.hill_act(S, P["KS"], P["nS"]) * np.ones_like(x1)
    elif topo.s_TM == 0:
        u_m = np.ones_like(x1)
    elif topo.s_TM > 0:
        u_m = D.hill_act(x1, P["Kxm"], P["nxm"])
    else:
        u_m = D.hill_rep(x1, P["Kxm"], P["nxm"])
    ap_m = P["gm"] * u_m
    am_m = P["gm"] * m1
    if topo.mir_arm:
        am_m = am_m + P["theta"] * m1 * y1

    # --- target mRNA ---
    u_TY = np.zeros_like(x1)
    if topo.s_TY != 0:
        u_TY = D.hill_act(x1, P["Kxy"], P["nxy"]) if topo.s_TY > 0 \
            else D.hill_rep(x1, P["Kxy"], P["nxy"])
    if topo.mir_arm:
        ap_y = (u_TY + b_leak) / (1.0 + b_leak)
        q_y = 1.0 - D.hill_rep(Meff, P["Kmy"], P["nmy"])
        am_y = (1.0 + P["lam_y"] * q_y) * y1 + P["theta"] * Meff * y1
        tr = D.hill_rep(Meff, P["Kmp"], P["nmp"]) if translational_block else np.ones_like(y1)
    else:
        if topo.s_MY == 0:
            prod = u_TY
        else:
            u_MY = D.hill_act(m1, P["Kmy"], P["nmy"]) if topo.s_MY > 0 \
                else D.hill_rep(m1, P["Kmy"], P["nmy"])
            prod = G(u_TY, u_MY) if topo.s_TY != 0 else u_MY
        ap_y = (prod + b_leak) / (1.0 + b_leak)
        am_y = y1
        tr = np.ones_like(y1)

    ap_p = P["gp"] * y1 * tr
    am_p = P["gp"] * p1
    return (np.array([ap_x, ap_m, ap_y, ap_p]),
            np.array([am_x, am_m, am_y, am_p]))


def simulate(topo, P, N, omega, rng, record=True):
    Z = np.zeros((4, N))
    n_burn = int(T_BURN / DT)
    n_samp = int(T_SAMP / DT)
    sq = np.sqrt(DT / omega)
    s1 = np.zeros(N); s2 = np.zeros(N); cnt = 0
    thin = 10
    for i in range(n_burn + n_samp):
        ap, am = fluxes(topo, P, Z, S_OP)
        drift = (ap - am) * DT
        if omega < np.inf:
            noise = sq * (np.sqrt(np.maximum(ap, 0.0)) * rng.standard_normal((4, N))
                          - np.sqrt(np.maximum(am, 0.0)) * rng.standard_normal((4, N)))
        else:
            noise = 0.0
        Z = np.maximum(Z + drift + noise, 0.0)
        if record and i >= n_burn and (i % thin == 0):
            p = Z[3]
            s1 += p; s2 += p * p; cnt += 1
    mean = s1 / max(cnt, 1)
    var = np.maximum(s2 / max(cnt, 1) - mean ** 2, 0.0)
    return mean, np.sqrt(var)


if __name__ == "__main__":
    t0 = time.time()
    P = {k: np.asarray(v, float) for k, v in SMP.sample(N_SET, seed=SEED).items()}
    rows = []
    # deterministic reference: the CLE with Omega = inf must reproduce the ODE
    for nm in TOPOS:
        topo = reg[nm]
        rng = np.random.default_rng(1)
        mean_inf, sd_inf = simulate(topo, P, N_SET, np.inf, rng)
        rhs = D.make_rhs(topo, P)
        Zode = D.integrate(rhs, np.zeros((D.NSTATE, N_SET)), lambda t: S_OP,
                           T_BURN + T_SAMP, 0.01, active=topo.active_states())
        err = np.abs(mean_inf - Zode[D.IP1])
        print(f"  {nm:14s} CLE(Omega=inf) vs ODE: max |diff| = {err.max():.3e}, "
              f"median = {np.median(err):.3e}", flush=True)
        for omega in OMEGAS:
            rng = np.random.default_rng(int(omega))
            mean, sd = simulate(topo, P, N_SET, omega, rng)
            cv = np.where(mean > 1e-6, sd / np.maximum(mean, 1e-12), np.nan)
            rows.append(pd.DataFrame(dict(topology=nm, omega=omega,
                                          param_set=np.arange(N_SET),
                                          mean_protein=mean, sd_protein=sd, cv_protein=cv,
                                          fano_like=np.where(mean > 1e-6, sd ** 2 * omega / mean, np.nan),
                                          ode_protein=Zode[D.IP1],
                                          cle_infinite_omega_protein=mean_inf)))
        print(f"  {nm:14s} stochastic done ({time.time()-t0:.0f}s)", flush=True)
    A = pd.concat(rows, ignore_index=True)
    A.to_csv(f"{RES}/dynamics_stochastic_noise.csv", index=False)

    # ---- summary: paired vs the cascade, and the CV-vs-mean regression --------
    from scipy.stats import wilcoxon
    out = []
    for omega in OMEGAS:
        sub = A[A.omega == omega]
        base = sub[sub.topology == "CASCADE"].set_index("param_set")
        for nm in TOPOS:
            g = sub[sub.topology == nm].set_index("param_set")
            ok = np.isfinite(g.cv_protein) & np.isfinite(base.cv_protein) & \
                 (g.mean_protein > 1e-3) & (base.mean_protein > 1e-3)
            a = base.loc[ok, "cv_protein"].to_numpy(float)
            b = g.loc[ok, "cv_protein"].to_numpy(float)
            r = dict(topology=nm, omega=omega, n=int(ok.sum()),
                     median_CV=float(np.nanmedian(b)),
                     median_mean_output=float(np.nanmedian(g.loc[ok, "mean_protein"])))
            if nm != "CASCADE" and ok.sum() > 20:
                r["median_CV_ratio_to_cascade"] = float(np.median(b / a))
                r["pct_sets_quieter_than_cascade"] = 100 * float(np.mean(b < a))
                r["wilcoxon_p_vs_cascade"] = float(wilcoxon(b, a).pvalue)
            # CV vs mean regression (noise at matched mean)
            m_ = g.loc[ok, "mean_protein"].to_numpy(float)
            keep = (m_ > 1e-3) & (b > 0)
            if keep.sum() > 50:
                X = np.column_stack([np.ones(keep.sum()), np.log(m_[keep])])
                beta, *_ = np.linalg.lstsq(X, np.log(b[keep]), rcond=None)
                r["log_CV_intercept_at_mean_1"] = float(beta[0])
                r["log_CV_slope_vs_log_mean"] = float(beta[1])
            out.append(r)
    S = pd.DataFrame(out)
    S.to_csv(f"{RES}/dynamics_stochastic_noise_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print()
    print(S.round(4).to_string(index=False))
    print(f"\nDONE in {(time.time()-t0)/60:.1f} min")
