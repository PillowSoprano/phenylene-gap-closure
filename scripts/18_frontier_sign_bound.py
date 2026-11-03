"""Frontier-sign upper bound rho_F (Section 5.1): the Poljak-Rohn family for the single sign pair
y = sgn(u), z = sgn(v) of the frontier singular vectors. Compared with the exact rho (Table 1) and,
by L/A sequence for N <= 6, with the Kekule upper bound of Table 2."""
import sys, json, itertools
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import phenylene, triangular4, benzocyclobutadiene, bipartite

def rho_F(T):
    U, s, Vt = np.linalg.svd(T); u, v = U[:, -1], Vt[-1]
    y = np.where(u >= 0, 1., -1.); z = np.where(v >= 0, 1., -1.)
    lam = np.linalg.eigvals(np.linalg.solve(T, T * np.outer(y, z)))
    r = np.abs(lam.real[np.abs(lam.imag) < 1e-9]); r = r[r > 1]
    return 1 / r.max() if len(r) else 1.0

exact = json.load(open(ROOT / 'results' / 'thresholds.json'))
mols = {'benzocyclobutadiene': benzocyclobutadiene(), 'biphenylene': phenylene([]), 'linear[3]': phenylene(['p']),
        'angular[3]': phenylene(['L']), 'linear[4]': phenylene(['p', 'p']), 'bent[4]': phenylene(['p', 'L']),
        'helical[4]': phenylene(['L', 'R']), 'zigzag[4]': phenylene(['L', 'L']), 'triangular[4]': triangular4()}
out = {'skeletons': {}, 'classes': {}}
print('skeleton               rho        rho_F')
for k, g in mols.items():
    T, *_ = bipartite(g); r = rho_F(T); sv = np.linalg.svd(T, compute_uv=False)
    degenerate = bool(sv[-2] - sv[-1] < 1e-9)   # triangular [4]: u, v not unique, rho_F depends on the basis returned
    out['skeletons'][k] = dict(rho=exact[k]['rho'], rho_F=r, sigma_min_degenerate=degenerate)
    assert r >= exact[k]['rho'] - 1e-12
    print(f'{k:20s} {exact[k]["rho"]:.6f}   {r:.6f}' + ('   equal' if abs(r - exact[k]['rho']) < 1e-9 else '') + ('   (sigma_min degenerate)' if degenerate else ''))
old = json.load(open(ROOT / 'results' / 'bounds_by_LA_sequence.json'))
print('\nclass      lower     rho_Kek   rho_F(min over isomers)  new upper  width ratio')
for N in (2, 3, 4, 5, 6):
    cls = {}
    for ch in itertools.product('pLR', repeat=N - 2):
        seq = ''.join('L' if c == 'p' else 'A' for c in ch); cls.setdefault(min(seq, seq[::-1]), []).append(ch)
    for key, chs in sorted(cls.items()):
        rf = min(rho_F(bipartite(phenylene(list(c)))[0]) for c in chs); o = old[f'{N}:{key or "-"}']
        lo = max(o['L_gap_max'], o['L_BS_max']); up = min(o['U_Kek'], rf)
        out['classes'][f'{N}:{key or "-"}'] = dict(lower=lo, U_Kek=o['U_Kek'], rho_F=rf, upper=up, width_ratio=(up - lo) / (o['U_Kek'] - lo))
        print(f'{N}:{key or "-":6s} {lo:.4f}   {o["U_Kek"]:.4f}   {rf:.4f}                  {up:.4f}     {(up - lo) / (o["U_Kek"] - lo):.3f}')
json.dump(out, open(ROOT / 'results' / 'frontier_sign_bound.json', 'w'), indent=1)

# the frontier sign pattern y_i z_j is symmetric under the automorphism that exchanges the two ends
# for angular [3]phenylene (sequence A) and for both [5]phenylenes with sequence LAL
import networkx as nx
print('\nfrontier sign pattern symmetric under end exchange:')
for name, ch in (('angular[3] (A)', ['L']), ('[5] LAL, p-L-p', ['p', 'L', 'p']), ('[5] LAL, p-R-p', ['p', 'R', 'p'])):
    g = phenylene(ch); T, L, R, bix = bipartite(g); U, s, Vt = np.linalg.svd(T)
    pat = {frozenset(e): np.sign(U[i, -1] * Vt[-1, j]) for e, (i, j) in zip(g['edges'], bix)}
    hexes = [r for r in g['rings'] if len(r) == 6]; G = nx.Graph(g['edges'])
    m = next(m for m in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter() if set(m[a] for a in hexes[0]) == set(hexes[-1]))
    sym = all(pat[frozenset(e)] == pat[frozenset((m[e[0]], m[e[1]]))] for e in g['edges'])
    print(f'  {name}: {sym}, rho_F = {rho_F(T):.4f}'); out.setdefault('symmetric_pattern', {})[name] = bool(sym)
json.dump(out, open(ROOT / 'results' / 'frontier_sign_bound.json', 'w'), indent=1)
