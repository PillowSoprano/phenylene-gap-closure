"""Section 3: the threshold of biphenylene at its relaxed Hueckel-SSH ground state (0.310) expressed as
bond-length changes, compared with the spread of the relaxed bond lengths."""
import numpy as np
R0, X, Y, RREF = 1.523, 0.189, 0.2756, 1.397
RHO = 0.3100                                  # rho at the relaxed ground state (geometry/mech.py)
# biphenylene: atoms 0-5 and 6-11, bridging bonds 0-6 and 1-7
bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6), (0, 6), (1, 7)]
n = 12; nh = 6; ia = np.array([a for a, _ in bonds]); ib = np.array([b for _, b in bonds])
R = np.full(len(bonds), RREF)
for _ in range(20000):
    H = np.zeros((n, n)); t = 1 + (RREF - R) / Y; H[ia, ib] = -t; H[ib, ia] = -t
    e, C = np.linalg.eigh(H); w = np.zeros(n); w[:nh] = 2
    Rn = R0 - X * ((C[ia, :] * C[ib, :]) @ w)
    if np.max(abs(Rn - R)) < 1e-12: break
    R = 0.6 * R + 0.4 * Rn
t = 1 + (RREF - R) / Y
# a relative coupling change rho*t_e changes the length by rho*t_e*y (coupling linear in length)
print(f'relaxed bond lengths {R.min():.3f} to {R.max():.3f} A; bridging bonds {R[-1]:.3f} A')
print(f'largest bond-length change at rho = {RHO}: {np.max(RHO * t * Y):.3f} A; at equal couplings (rho = 0.1888): {0.1888 * Y:.3f} A')
