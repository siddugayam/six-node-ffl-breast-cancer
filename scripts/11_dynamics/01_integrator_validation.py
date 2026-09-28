#!/usr/bin/env python3
"""
01_integrator_validation.py
Cross-check the vectorised fixed-step RK4 used everywhere else against
scipy.integrate.solve_ivp (LSODA, rtol=1e-9, atol=1e-12) on every topology,
and against itself at half the step size.  Out: results/v3/dynamics_integrator_validation.csv
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyn_models as D
from scipy.integrate import solve_ivp

OUT = "/path/to/revision/results/v3/dynamics_integrator_validation.csv"

def run_scipy(topo, P, S, t_end):
    rhs = D.make_rhs(topo, P)
    act = topo.active_states()
    def f(t, z):
        Z = z.reshape(D.NSTATE, 1)
        dZ = rhs(t, Z, S)
        out = np.zeros(D.NSTATE)
        out[act] = dZ[act, 0]
        return out
    sol = solve_ivp(f, (0, t_end), np.zeros(D.NSTATE), method="LSODA",
                    rtol=1e-9, atol=1e-12, dense_output=True)
    return sol.y[:, -1]

rows = []
reg = D.build_registry()
ho  = D.higher_order_registry("AND", "I1_MIR", composite=True)
reg.update({t.name: t for t in ho.values()})
P = D.default_params(1)
P["theta"] = np.array([2.0]); P["kappa"] = np.array([1.5])   # exercise the nonlinear terms
for S in (0.5, 2.0):
    for nm, t in reg.items():
        ref = run_scipy(t, P, S, 80.0)
        for dt in (0.02, 0.01, 0.005):
            Z = D.steady_state(D.make_rhs(t, P), np.zeros((D.NSTATE, 1)), S, 80.0, dt,
                               active=t.active_states())
            err = np.abs(Z[:, 0] - ref)
            rows.append(dict(topology=nm, S=S, dt=dt,
                             max_abs_err=float(err.max()),
                             max_rel_err=float((err / np.maximum(np.abs(ref), 1e-6)).max()),
                             p1_rk4=float(Z[D.IP1, 0]), p1_lsoda=float(ref[D.IP1])))
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(df.groupby("dt")[["max_abs_err", "max_rel_err"]].max())
print("\nworst rows at dt=0.01:")
print(df[df.dt == 0.01].nlargest(5, "max_abs_err").to_string(index=False))
print("\nwritten:", OUT)
