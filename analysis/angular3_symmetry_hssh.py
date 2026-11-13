"""Hueckel-SSH triplet of angular [3]phenylene: symmetric versus broken-symmetry relaxation (Section 8).
Prints the lowering of the triplet in |beta| and the sign agreement of its bond-length changes with the
closing patterns (fraction of bonds and fraction of the total change)."""
import json, os, numpy as np, networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
here = os.path.dirname(os.path.abspath(__file__))
R0, X, Y, RREF = 1.523, 0.189, 0.2756, 1.397          # Kwapisz-Stolarczyk parameters
K = 2 / (X * Y)                                        # energies in units of t0 = |beta|
d = {r['id']: r for r in json.load(open(os.path.join(here, '..', 'beyond_huckel', 'patterns_N2-9.json')))}['PT03-00001']
E = [(a, b) for a, b, _ in d['cp']]; s = np.array([c for *_, c in d['cp']], float)
n = 6 * d['N']; nh = n // 2; ia = np.array([a for a, _ in E]); ib = np.array([b for _, b in E])

def solve(R, occ):
    H = np.zeros((n, n)); t = 1 + (RREF - R) / Y; H[ia, ib] = -t; H[ib, ia] = -t
    e, C = np.linalg.eigh(H); w = np.zeros(n); w[:nh] = 2
    if occ == 'T1': w[nh - 1] = 1; w[nh] = 1
    return (C[ia, :] * C[ib, :]) @ w, float(w @ e)

def relax(occ, R):
    for _ in range(20000):
        p, _ = solve(R, occ); Rn = R0 - X * p
        if np.max(abs(Rn - R)) < 1e-12: break
        R = 0.6 * R + 0.4 * Rn
    p, eel = solve(R, occ)
    return R, eel + K / 2 * np.sum((R - R0) ** 2)

G = nx.Graph(E); idx = {frozenset(e): k for k, e in enumerate(E)}
perms = [np.array([idx[frozenset((a[u], a[v]))] for u, v in E]) for a in GraphMatcher(G, G).isomorphisms_iter()]

def symmetrize(R):                     # average over the automorphisms of the skeleton
    out = np.zeros_like(R)
    for pm in perms:
        Ri = np.empty_like(R); Ri[pm] = R; out += Ri
    return out / len(perms)

def relax_sym(occ):
    R = np.full(len(E), RREF)
    for _ in range(20000):
        p, _ = solve(R, occ); Rn = symmetrize(R0 - X * p)
        if np.max(abs(Rn - R)) < 1e-12: break
        R = 0.6 * R + 0.4 * Rn
    p, eel = solve(R, occ)
    return R, eel + K / 2 * np.sum((R - R0) ** 2)

RS, ES = relax('S0', np.full(len(E), RREF))
RTs, ETs = relax_sym('T1')                                           # best symmetric triplet
rng = np.random.default_rng(0); best = None
for k in range(5):
    R, Et = relax('T1', RS + rng.normal(0, 0.02, len(E)))
    if best is None or Et < best[1]: best = (R, Et)
RTb, ETb = best
print('E_T1 symmetric %.4f, broken %.4f |beta|, lowering %.4f |beta|' % (ETs - ES, ETb - ES, ETs - ETb))
# closing patterns: the stored one and its images under the automorphisms of the skeleton
dR = RTb - RS; tot = np.sum(np.abs(dR))
for perm in perms:
    sp = np.empty_like(s); sp[perm] = s
    ok = np.sign(-dR) == sp
    print('pattern image: agree on %d/%d bonds (%.0f%%), carrying %.0f%% of the change'
          % (ok.sum(), len(E), 100 * ok.mean(), 100 * np.sum(np.abs(dR[ok])) / tot))
