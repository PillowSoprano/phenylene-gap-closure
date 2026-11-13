"""Closure thresholds of zigzag and helical [4]phenylene at their relaxed Hueckel-SSH ground states (Section 7).
Exact by complete sign enumeration (Poljak-Rohn), about three minutes per molecule."""
import os, sys, itertools, numpy as np, networkx as nx
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'gapclosure', 'src'))
from gapclosure.skeleton import phenylene_tree
R0, X, Y, RREF = 1.523, 0.189, 0.2756, 1.397

def relaxed_couplings(bonds, n):
    nh = n // 2; ia = np.array([a for a, _ in bonds]); ib = np.array([b for _, b in bonds])
    R = np.full(len(bonds), RREF)
    for _ in range(20000):
        H = np.zeros((n, n)); t = 1 + (RREF - R) / Y; H[ia, ib] = -t; H[ib, ia] = -t
        e, C = np.linalg.eigh(H); w = np.zeros(n); w[:nh] = 2
        Rn = R0 - X * ((C[ia, :] * C[ib, :]) @ w)
        if np.max(abs(Rn - R)) < 1e-12: break
        R = 0.6 * R + 0.4 * Rn
    return 1 + (RREF - R) / Y

def rho(T):
    n = T.shape[0]; Ti = np.linalg.inv(T); best = 0
    signs = np.array(list(itertools.product([-1, 1], repeat=n)))
    for z in signs[signs[:, 0] == 1]:
        ev = np.linalg.eigvals(Ti[None] @ (T[None] * signs[:, :, None] * z[None, None, :]))
        best = max(best, np.where(np.abs(ev.imag) < 1e-9, np.abs(ev.real), 0).max())
    return 1 / best

for name, links in (('zigzag', [[0, 0, 1], [0, 2, 2], [1, 4, 3]]), ('helical', [[0, 0, 1], [0, 2, 2], [1, 2, 3]])):
    sk = phenylene_tree(links); bonds = [tuple(b) for b in sk.bonds]
    t = relaxed_couplings(bonds, sk.n)
    G = nx.Graph(bonds); col = nx.bipartite.color(G)
    S = [v for v in G if col[v] == 0]; U = [v for v in G if col[v] == 1]
    si = {v: i for i, v in enumerate(S)}; ui = {v: i for i, v in enumerate(U)}
    T = np.zeros((len(S), len(U)))
    for (a, b), w in zip(bonds, t):
        if col[a] == 1: a, b = b, a
        T[si[a], ui[b]] = w
    print(f'{name} [4]phenylene, relaxed ground state: rho = {rho(T):.5f}', flush=True)
