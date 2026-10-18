import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np, time
from itertools import product
from gapclosure import phenylene, triangular4, bipartite, kekule_terms
def all_patterns(g, tol=1e-9):
    T, L, R, bix = bipartite(g); n = len(T); Ti = np.linalg.inv(T)
    ys = np.array(list(product([-1, 1], repeat=n)), float)
    zs = np.hstack([np.ones((2**(n-1), 1)), np.array(list(product([-1, 1], repeat=n-1)), float)])
    A = np.einsum('ij,yj,jk->yik', Ti, ys, T); best, found = -np.inf, []
    zper = max(1, (1 << 15)//len(ys))
    for zb in range(0, len(zs), zper):
        Z = zs[zb:zb+zper]
        lam = np.linalg.eigvals((A[None]*Z[:, None, None, :]).reshape(-1, n, n))
        lm = np.where(np.abs(lam.imag) < 1e-9, np.abs(lam.real), -np.inf).max(1)
        m = lm.max()
        if m > best + tol: best, found = m, []
        for k in np.nonzero(lm > best - tol)[0]:
            iz, iy = divmod(k, len(ys)); found.append((ys[iy].copy(), Z[iz].copy()))
    rho = 1/best; out = {}
    for y, z in found:
        S = T*np.outer(y, z)
        # sign of the closing direction: det(T - eps S) = 0 at eps = 1/lambda; pick lambda sign
        lam = np.linalg.eigvals(Ti @ S); lam = lam[np.abs(lam.imag) < 1e-9].real
        sgn = 1 if np.isclose(np.abs(lam).max(), lam.max()) else -1
        W = T - sgn*rho*S
        U, s, Vt = np.linalg.svd(W); w = np.zeros(g['n']); w[L] = U[:, -1]**2; w[R] = Vt[-1]**2; w /= 2
        sh = np.array([sum(w[a] for a in rg) for rg in g['rings']]); sh /= sh.sum()
        out[tuple(np.round(sh, 2))] = s[-1]
    return rho, len(found), out
for name, g in (('bent[4]', phenylene(['p', 'L'])), ('triangular[4]', triangular4()), ('zigzag[4]', phenylene(['L','L'])), ('helical[4]', phenylene(['L','R']))):
    t = time.time(); rho, nf, out = all_patterns(g)
    print(f'{name}: rho={rho:.12f}, optimal sign-vector pairs={nf} (z1 fixed), distinct closing-orbital ring shares:', flush=True)
    for k, v in out.items(): print('    ', k, f'sigma_min={v:.1e}')
    print(f'   ({time.time()-t:.0f}s)', flush=True)
