"""Linear phenylene polymer at equal couplings: band edges, and gaps of long linear, zigzag and helical chains."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import phenylene, bipartite
def H_k(k):                     # unit cell: one hexagon 0..5; para bonds (0,1) and (3,4); square 3-4-0'-1'
    H = np.zeros((6, 6), complex)
    for i in range(6): j = (i + 1) % 6; H[i, j] = H[j, i] = -1
    for a, b in ((3, 1), (4, 0)): H[a, b] += -np.exp(1j*k); H[b, a] += -np.exp(-1j*k)
    return H
ks = np.linspace(-np.pi, np.pi, 4001); E = np.array([np.linalg.eigvalsh(H_k(k)) for k in ks])
print(f'linear polymer: top of band 3 = {E[:,2].max():+.2e} at k = {ks[E[:,2].argmax()]:+.3f}; '
      f'bottom of band 4 = {E[:,3].min():+.2e} at k = {ks[E[:,3].argmin()]:+.3f}  (indirect gap {E[:,3].min()-E[:,2].max():.2e})')
def gap(seq):
    T, *_ = bipartite(phenylene(seq)); return 2*np.linalg.svd(T, compute_uv=False)[-1]
for N in (10, 20, 40, 80, 160):
    g = gap(['p']*(N-2)); print(f'linear  N={N:4d}  gap={g:.4e}  N^2*gap={N*N*g:.3f}')
for N in (6, 10, 20, 40):
    gz = gap(['L']*(N-2)); gh = gap((['L', 'R']*N)[:N-2])
    print(f'N={N:3d}  zigzag gap={gz:.4f}  helical gap={gh:.4f}')
