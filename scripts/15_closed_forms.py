"""Closed forms of Theorem 2 and the spectrum of linear [3]phenylene (Sections 3.1, 5.1 and 5.3).

1. For a skeleton with one minority structure M-, Theorem 2 gives rho from sum_a q^(k_a) = 1, q = (1-rho)/(1+rho),
   where k_a is the number of double bonds of the a-th majority structure that M- does not contain.
   The k_a are read off the Kekule structures, the equation is solved, and the result is compared with the closed
   forms q^2 = 1/phi (benzocyclobutadiene) and q^3 + q = 1 (biphenylene) and with rho in results/thresholds.json.
2. Every Hueckel eigenvalue of benzocyclobutadiene is an eigenvalue of linear [3]phenylene, and the two skeletons
   have the same gap and the same threshold.
"""
import sys, json
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from scipy.optimize import brentq
from gapclosure import phenylene, benzocyclobutadiene, bipartite, kekule_terms

rhos = json.load(open(ROOT / 'results' / 'thresholds.json'))
phi = (1 + 5**0.5) / 2
q_bcb = phi**-0.5                                       # q^2 = 1/phi
q_biph = brentq(lambda q: q**3 + q - 1, 0, 1, xtol=1e-16)  # q^3 + q = 1
closed = {'benzocyclobutadiene': ('q^2 = 1/phi', q_bcb), 'biphenylene': ('q^3 + q = 1', q_biph)}

print('1. One minority structure (Theorem 2)')
for name, g in [('benzocyclobutadiene', benzocyclobutadiene()), ('biphenylene', phenylene([]))]:
    T, *_ = bipartite(g)
    terms = kekule_terms(T)
    minus = [set(enumerate(p)) for s, p in terms if s < 0]
    assert len(minus) == 1
    k = sorted(len(set(enumerate(p)) - minus[0]) for s, p in terms if s > 0)
    q = brentq(lambda q: sum(q**ka for ka in k) - 1, 1e-9, 1, xtol=1e-16)
    rho = (1 - q) / (1 + q)
    form, qc = closed[name]
    rho_c = (1 - qc) / (1 + qc)
    print(f'  {name}: K = {len(terms)}, k_a = {k}')
    print(f'    root of sum q^k_a = 1: q = {q:.15f}, rho = {rho:.15f}')
    print(f'    closed form {form}: q = {qc:.15f}, rho = {rho_c:.15f}')
    print(f'    rho in thresholds.json: {rhos[name]["rho"]:.15f}   |difference| = {abs(rho_c - rhos[name]["rho"]):.1e}')

print('\n2. Linear [3]phenylene and benzocyclobutadiene')
def spectrum(g):
    A = np.zeros((g['n'], g['n']))
    for u, v in g['edges']: A[u, v] = A[v, u] = 1
    return np.linalg.eigvalsh(A)
s_bcb, s_lin = spectrum(benzocyclobutadiene()), spectrum(phenylene(['p']))
r = lambda x: round(float(x), 9) + 0.0
missing = Counter(map(r, s_bcb)) - Counter(map(r, s_lin))
print(f'  spectrum of benzocyclobutadiene ({len(s_bcb)} eigenvalues) contained in that of linear [3]phenylene '
      f'({len(s_lin)} eigenvalues), with multiplicity: {not missing}')
gap = lambda s: s[s > 1e-12].min() - s[s < -1e-12].max()
print(f'  gap: benzocyclobutadiene {gap(s_bcb):.12f}, linear [3]phenylene {gap(s_lin):.12f}  (units of |beta|)')
print(f'  rho: benzocyclobutadiene {rhos["benzocyclobutadiene"]["rho"]:.12f}, '
      f'linear [3]phenylene {rhos["linear[3]"]["rho"]:.12f}')
