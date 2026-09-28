#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02j_noise_at_matched_mean.py
============================
NOISE COMPARED AT MATCHED MEAN OUTPUT.

The raw paired comparison in dynamics_stochastic_noise_summary.csv says the
incoherent circuits are noisier than a cascade, but they are also much dimmer
(median target protein 0.21 vs 0.91), and the coefficient of variation falls as
roughly mean^-1/2 for any Poisson-limited process.  Osella et al. (2011,
doi:10.1371/journal.pcbi.1001101) make their claim AT MATCHED MEAN, so this
script makes the same comparison two ways:

  (1) EXCESS NOISE.  Fit log CV = a + b log(mean) on the cascade alone; for every
      circuit parameter set report log CV(circuit) minus the cascade's predicted
      log CV at the circuit's own mean.  Positive = noisier than a cascade
      producing the same amount of protein.
  (2) BINNED.  Within deciles of mean output, compare the median CV of each
      circuit with the median CV of the cascade sets falling in the same bin.

Out: results/v3/dynamics_stochastic_noise_matched_mean.csv
"""
import numpy as np, pandas as pd

RES = "/path/to/revision/results/v3"
A = pd.read_csv(f"{RES}/dynamics_stochastic_noise.csv")
A = A[(A.mean_protein > 1e-3) & (A.cv_protein > 0) & np.isfinite(A.cv_protein)]

rows = []
for omega, sub in A.groupby("omega"):
    cas = sub[sub.topology == "CASCADE"]
    x = np.log(cas.mean_protein.to_numpy()); y = np.log(cas.cv_protein.to_numpy())
    b, a = np.polyfit(x, y, 1)                      # slope, intercept
    edges = np.percentile(np.log(cas.mean_protein), np.linspace(0, 100, 11))
    for nm, g in sub.groupby("topology"):
        lm = np.log(g.mean_protein.to_numpy())
        lc = np.log(g.cv_protein.to_numpy())
        excess = lc - (a + b * lm)
        # binned comparison
        binned = []
        for i in range(len(edges) - 1):
            m1 = (lm >= edges[i]) & (lm < edges[i + 1])
            m2 = (np.log(cas.mean_protein) >= edges[i]) & (np.log(cas.mean_protein) < edges[i + 1])
            if m1.sum() >= 20 and m2.sum() >= 20:
                binned.append(np.median(np.exp(lc[m1])) / np.median(cas.cv_protein[m2]))
        rows.append(dict(
            topology=nm, omega=omega, n=len(g),
            cascade_logCV_intercept=a, cascade_logCV_slope=b,
            median_mean_output=float(np.median(g.mean_protein)),
            median_CV=float(np.median(g.cv_protein)),
            median_excess_noise_log=float(np.median(excess)),
            median_CV_ratio_at_matched_mean=float(np.exp(np.median(excess))),
            pct_sets_quieter_at_matched_mean=100 * float(np.mean(excess < 0)),
            n_mean_bins_compared=len(binned),
            median_binned_CV_ratio=float(np.median(binned)) if binned else np.nan,
        ))
R = pd.DataFrame(rows).sort_values(["omega", "topology"])
R.to_csv(f"{RES}/dynamics_stochastic_noise_matched_mean.csv", index=False)
pd.set_option("display.width", 250)
print(R.round(4).to_string(index=False))
