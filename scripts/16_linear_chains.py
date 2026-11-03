"""Linear [N]phenylenes (Theorem 3, Appendix B): the ladder form of T, the block structure of its inverse,
the explicit closing witness, the analytic bounds on rho, and the gap asymptotics."""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import phenylene, bipartite

def ladder(N):
    """T of linear [N]phenylene with both atom sets ordered along the chain: cell i holds
    positions a_i, b_i, c_i; diagonal entries 1, 0, 1 are the rungs (the two bonds of hexagon i
    perpendicular to the chain); the off-diagonal entries are the rail bonds."""
    d = np.tile([1.0, 0.0, 1.0], N); n = 3 * N
    return np.diag(d) + np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)

def G_matrix(N):
    i = np.arange(1, N + 1)
    return np.minimum.outer(i, i) * (N + 1 - np.maximum.outer(i, i)) / (N + 1)

def frontier_sign_bound(T):
    U, s, Vt = np.linalg.svd(T); y = np.where(U[:, -1] >= 0, 1., -1.); z = np.where(Vt[-1] >= 0, 1., -1.)
    lam = np.linalg.eigvals(np.linalg.solve(T, T * np.outer(y, z)))
    r = np.abs(lam.real[np.abs(lam.imag) < 1e-9]); r = r[r > 1]
    return 1 / r.max() if len(r) else 1.0

exact = json.load(open(ROOT / 'results' / 'thresholds.json'))
exact_by_N = {2: exact['biphenylene']['rho'], 3: exact['linear[3]']['rho'], 4: exact['linear[4]']['rho']}
out = {}
print('ladder form equals T of linear [N]phenylene (singular values, |det|, bond count):')
for N in range(2, 9):
    T, *_ = bipartite(phenylene(['p'] * (N - 2))); J = ladder(N)
    ok = np.allclose(np.sort(np.linalg.svd(T, compute_uv=False)), np.sort(np.linalg.svd(J, compute_uv=False)))
    print(f'  N={N}: {ok}, |det T|={abs(np.linalg.det(T)):.6f} = N+1, bonds {np.count_nonzero(T)} = 8N-2: {np.count_nonzero(J) == 8*N-2}')
    assert ok and abs(abs(np.linalg.det(J)) - (N + 1)) < 1e-9 and np.count_nonzero(J) == 8 * N - 2
print('inverse in the order (a, c; b): endpoint block [[G,-G],[-G,G]], middle block -11^T/(N+1), mixed entries <= 1')
for N in (2, 3, 5, 10, 40, 160):
    Ji = np.linalg.inv(ladder(N)); a = list(range(0, 3*N, 3)); c = list(range(2, 3*N, 3)); b = list(range(1, 3*N, 3)); G = G_matrix(N)
    ok1 = np.allclose(Ji[np.ix_(a + c, a + c)], np.block([[G, -G], [-G, G]]))
    ok2 = np.allclose(Ji[np.ix_(b, b)], -np.ones((N, N)) / (N + 1)); mx = np.abs(Ji[np.ix_(a + c, b)]).max()
    ok3 = abs(np.linalg.norm(G, 2) - 1 / (4 * np.sin(np.pi / (2 * (N + 1)))**2)) < 1e-9 * N * N
    print(f'  N={N}: endpoint block {ok1}, middle block {ok2}, max mixed entry {mx:.4f}, ||G|| = l_N {ok3}')
    assert ok1 and ok2 and mx <= 1 and ok3
print('closing witness: rungs x (1-eps), connecting bonds x (1+eps), eps = tan^2(theta_N)')
for N in (2, 3, 5, 10, 40):
    th = np.pi / (2 * (N + 1)); e = np.tan(th)**2; M = ladder(N)
    for i in range(N): M[3*i, 3*i] = 1 - e; M[3*i+2, 3*i+2] = 1 - e
    for i in range(N - 1): M[3*i+2, 3*i+3] = 1 + e; M[3*i+3, 3*i+2] = 1 + e
    smin = np.linalg.svd(M, compute_uv=False)[-1]; print(f'  N={N}: eps={e:.6f}, sigma_min={smin:.1e}'); assert smin < 1e-12
print(f"{'N':>4} {'lower(analytic)':>15} {'BS lower':>10} {'rho exact':>10} {'rho_F':>10} {'tan^2':>10} {'(N+1)^2 rho_F':>14} {'(N+1)^2 g':>10}")
for N in (2, 3, 4, 5, 6, 10, 20, 40, 80, 160):
    J = ladder(N); th = np.pi / (2 * (N + 1)); l = 1 / (4 * np.sin(th)**2)
    lo = 1 / (4*l + 2*np.sqrt(l) + 4*N + 1); up = np.tan(th)**2
    bs = 1 / max(abs(np.linalg.eigvals(np.abs(np.linalg.inv(J)) @ J))); rf = frontier_sign_bound(J)
    g = 2 * np.linalg.svd(J, compute_uv=False)[-1]; ex = exact_by_N.get(N)
    assert lo <= bs + 1e-12 and rf <= up + 1e-12 and (ex is None or bs <= ex + 1e-12 <= rf + 2e-12)
    out[N] = dict(lower_analytic=lo, lower_BS=bs, rho_exact=ex, rho_F=rf, upper_analytic=up, gap=g,
                  Np1sq_rho_F=(N+1)**2*rf, Np1sq_gap=(N+1)**2*g, Np1sq_lower_BS=(N+1)**2*bs)
    print(f"{N:4d} {lo:15.6f} {bs:10.6f} {('%.6f' % ex) if ex else '':>10} {rf:10.6f} {up:10.6f} {(N+1)**2*rf:14.5f} {(N+1)**2*g:10.5f}")
print(f'limits: pi^2/4 = {np.pi**2/4:.5f}, pi^2 = {np.pi**2:.5f}')
json.dump({str(k): v for k, v in out.items()}, open(ROOT / 'results' / 'linear_chains.json', 'w'), indent=1)
