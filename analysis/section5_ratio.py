"""Section 5: ratio of the triplet energy to the closure energy U of the nearest closed-gap geometry, the
one-dimensional model of Eq. (9) and the Spearman correlations of rho and U (Hueckel-SSH, N = 4 to 8).
Input: ../beyond_huckel/l2nearest.npy, columns N, E_T1, U (nearest closed gap, multi-start),
corner bound, rho, g0, and the gap at the triplet geometry (all in |beta|), written by beyond_huckel/l2nearest.py."""
import os, numpy as np
from scipy.stats import spearmanr, pearsonr
here = os.path.dirname(os.path.abspath(__file__))
N, E, U, Uc, rho, g0, _ = np.load(os.path.join(here, '..', 'beyond_huckel', 'l2nearest.npy')).T
m = (N >= 4) & (N <= 8)
r = E[m] / U[m]
print(f'classes with 4 <= N <= 8: {m.sum()}')
print(f'E_T1/U: min {r.min():.3f}, max {r.max():.3f}, mean {r.mean():.3f}')
for n in range(4, 9):
    k = N[m] == n; print(f'  N = {n}: mean {r[k].mean():.3f}')
x = g0[m] / (2 * U[m]); print(f'r = g0/2U: min {x.min():.3f}, max {x.max():.3f}, mean {x.mean():.3f}')
model = g0[m] * (1 - g0[m] / (4 * U[m]))
print(f'Eq. (9): largest relative error {np.max(np.abs(model - E[m]) / E[m]):.3f}, correlation {pearsonr(model, E[m])[0]:.4f}')
m5 = (N >= 5) & (N <= 8)
print(f'Spearman with E_T1, N = 5 to 8: rho {spearmanr(rho[m5], E[m5]).correlation:.4f}, U {spearmanr(U[m5], E[m5]).correlation:.4f}')
