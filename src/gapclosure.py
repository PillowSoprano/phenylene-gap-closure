import numpy as np
from itertools import product
from scipy.optimize import brentq

# ---------- molecular graphs with drawing coordinates ----------
def hexagon():
    ang = np.pi/2 - np.arange(6)*np.pi/3
    return dict(n=6, edges=[(i,(i+1)%6) for i in range(6)],
                pos=[(np.cos(a), np.sin(a)) for a in ang],
                rings=[list(range(6))])

def fuse_polygon(g, u, v, k):
    """Fuse a regular k-gon onto existing bond (u,v), on the side away from the ring
    that already contains it. Returns (new graph, ring path [u, new..., v])."""
    pos = np.array(g['pos'], float)
    host = next(r for r in g['rings'] if u in r and v in r)
    c_host = pos[host].mean(0)
    pu, pv = pos[u], pos[v]
    mid, d = (pu+pv)/2, pv-pu
    nrm = np.array([-d[1], d[0]]); nrm /= np.linalg.norm(nrm)
    if np.dot(nrm, mid-c_host) < 0: nrm = -nrm
    L = np.linalg.norm(d); apo = L/(2*np.tan(np.pi/k)); c = mid + apo*nrm
    R = L/(2*np.sin(np.pi/k))
    a_u = np.arctan2(*(pu-c)[::-1]); a_v = np.arctan2(*(pv-c)[::-1])
    step = 2*np.pi/k
    def wrap(x): return (x+np.pi) % (2*np.pi) - np.pi
    sgn = 1 if abs(wrap(a_u + step*(k-1) - a_v)) < 1e-6 else -1
    new = list(range(g['n'], g['n']+k-2))
    newpos = [tuple(c + R*np.array([np.cos(a_u+sgn*step*(i+1)), np.sin(a_u+sgn*step*(i+1))]))
              for i in range(k-2)]
    path = [u]+new+[v]
    h = dict(n=g['n']+k-2, edges=g['edges']+list(zip(path[:-1], path[1:])),
             pos=list(g['pos'])+newpos, rings=g['rings']+[path])
    return h, path

def phenylene(choices):
    """[N]phenylene. choices: list of length N-2 with entries 'p' (para / linear),
    'L' or 'R' (angular to either side) for each internal benzene ring."""
    g = hexagon(); a, b = 2, 3          # bond of the first ring used for the first 4-ring
    for step in range(len(choices)+1):
        g, sq = fuse_polygon(g, a, b, 4)
        x, y = sq[2], sq[1]              # far bond of the square, oriented
        g, hx = fuse_polygon(g, x, y, 6)         # hx = [x, p1, p2, p3, p4, y]
        if step < len(choices):
            c = choices[step]
            a, b = {'p': (hx[2], hx[3]), 'L': (hx[1], hx[2]), 'R': (hx[3], hx[4])}[c]
    return g

def triangular4():
    g = hexagon()
    for a, b in [(0,1),(2,3),(4,5)]:
        g, sq = fuse_polygon(g, a, b, 4)
        g, _ = fuse_polygon(g, sq[2], sq[1], 6)
    return g

def benzocyclobutadiene():
    g = hexagon(); g, sq = fuse_polygon(g, 2, 3, 4); return g

def naphthalene():
    g = hexagon(); g, _ = fuse_polygon(g, 2, 3, 6); return g

def check_geometry(g, tol=1e-6):
    """Bonds all of unit length, no two atoms closer than 0.5 (rejects overlapping rings)."""
    P = np.array(g['pos'])
    for u, v in g['edges']:
        assert abs(np.linalg.norm(P[u]-P[v])-1) < tol
    D = np.linalg.norm(P[:,None]-P[None], axis=-1) + 10*np.eye(g['n'])
    assert D.min() > 0.5, 'overlapping atoms'
    deg = np.bincount(np.array(g['edges']).ravel(), minlength=g['n'])
    assert deg.max() <= 3

# ---------- bipartite block ----------
def bipartite(g):
    n, E = g['n'], g['edges']
    adj = [[] for _ in range(n)]
    for u, v in E: adj[u].append(v); adj[v].append(u)
    col, q = {0: 0}, [0]
    for u in q:
        for v in adj[u]:
            if v not in col: col[v] = 1-col[u]; q.append(v)
            assert col[v] != col[u], 'not bipartite'
    L = [i for i in range(n) if col[i] == 0]; R = [i for i in range(n) if col[i] == 1]
    assert len(L) == len(R), 'unbalanced'
    li = {a: i for i, a in enumerate(L)}; ri = {a: j for j, a in enumerate(R)}
    T = np.zeros((len(L), len(R))); bond_ix = []
    for u, v in E:
        if u in li: i, j = li[u], ri[v]
        else: i, j = li[v], ri[u]
        T[i, j] = 1; bond_ix.append((i, j))
    return T, L, R, bond_ix

def huckel_gap(T):
    return 2*np.linalg.svd(T, compute_uv=False)[-1]

# ---------- Kekulé structures and algebraic structure count ----------
def kekule_terms(T):
    m = len(T); out = []
    cols = [np.nonzero(T[i])[0] for i in range(m)]
    def rec(i, used, perm):
        if i == m:
            p = perm; inv = sum(p[a] > p[b] for a in range(m) for b in range(a+1, m))
            out.append(((-1)**inv, tuple(p))); return
        for j in cols[i]:
            if j not in used: rec(i+1, used | {j}, perm+[j])
    rec(0, frozenset(), [])
    sgn = int(np.sign(np.linalg.det(T)))
    return [(sgn*s, p) for s, p in out]

def structure_counts(T):
    terms = kekule_terms(T)
    K = len(terms); ASC = abs(sum(s for s, _ in terms))
    return K, ASC, sum(1 for s, _ in terms if s < 0)

# ---------- exact closure threshold rho (Rohn y,z-vertex criterion) ----------
def rho_rohn(T, chunk=1 << 15):
    """Exact regularity radius of the interval matrix [T - eps*T, T + eps*T] (zero pattern fixed).
    Rohn: the interval matrix is singular iff some A_yz = T - eps * D_y T D_z (y,z in {+-1}^n)
    is singular. For each (y,z), det(T - eps S) = 0  <=>  eps = 1/lambda, lambda a real
    eigenvalue of T^{-1} S. rho = min over (y,z) of 1/lambda_max (only lambda > 1 matter).
    Returns rho (1.0 if no closure for eps < 1) and the maximizing (y, z)."""
    n = len(T); Ti = np.linalg.inv(T)
    best, arg = 1.0, None
    ys = np.array(list(product([-1, 1], repeat=n)), float)          # 2^n
    zs = np.array(list(product([-1, 1], repeat=n-1)), float)        # z_0 fixed to +1
    zs = np.hstack([np.ones((len(zs), 1)), zs])
    A = np.einsum('ij,yj,jk->yik', Ti, ys, T)                        # Ti D_y T for every y
    zper = max(1, chunk // len(ys))
    for zb in range(0, len(zs), zper):
        Z = zs[zb:zb+zper]
        M = (A[None] * Z[:, None, None, :]).reshape(-1, n, n)        # (Ti D_y T) D_z
        lam = np.linalg.eigvals(M)
        lam = np.where(np.abs(lam.imag) < 1e-9, lam.real, -np.inf)
        lmax = lam.max(axis=1)
        k = int(np.argmax(lmax))
        if lmax[k] > max(1/best, 1 + 1e-4) + 1e-12:
            best = 1/lmax[k]
            iz, iy = divmod(k, len(ys)); arg = (ys[iy].copy(), Z[iz].copy())
    return best, arg

def witness_from_yz(T, rho, yz):
    y, z = yz
    W = T - rho*(T*np.outer(y, z))
    return W

# ---------- brute-force reference over all box vertices (small graphs only) ----------
def rho_vertices(T):
    terms = kekule_terms(T)
    edges = list(zip(*np.nonzero(T)))
    eix = {e: k for k, e in enumerate(edges)}
    tix = [(s, [eix[(i, p[i])] for i in range(len(T))]) for s, p in terms]
    if all(s > 0 for s, _ in tix): return 1.0
    m = len(edges)
    def mindet(eps):
        mn = np.inf
        for start in range(0, 1 << m, 1 << 16):
            idx = np.arange(start, min(start + (1 << 16), 1 << m))
            S = ((idx[:, None] >> np.arange(m)) & 1)*2-1
            d = sum(s*np.prod(1+eps*S[:, ix], axis=1) for s, ix in tix)
            mn = min(mn, d.min())
        return mn
    return brentq(mindet, 1e-6, 0.9999, xtol=1e-13)

def rho_local(T, restarts=200, seed=0):
    """Upper bound on rho by hill-climbing over Rohn sign vectors (y,z). Not a certificate."""
    n = len(T); Ti = np.linalg.inv(T); rng = np.random.default_rng(seed)
    def lam(y, z):
        l = np.linalg.eigvals(Ti @ (T*np.outer(y, z)))
        l = l.real[np.abs(l.imag) < 1e-9]
        return l.max() if len(l) else -np.inf
    best = -np.inf
    for _ in range(restarts):
        y = rng.choice([-1., 1.], n); z = rng.choice([-1., 1.], n); v = lam(y, z)
        improved = True
        while improved:
            improved = False
            for vec in (y, z):
                for k in range(n):
                    vec[k] *= -1; w = lam(y, z)
                    if w > v + 1e-12: v = w; improved = True
                    else: vec[k] *= -1
        best = max(best, v)
    return 1/best if best > 1 + 1e-6 else 1.0

# ---------- Hückel + sigma-spring (SSH-type) bond-length model ----------
from scipy.optimize import minimize

class BondModel:
    """t_b = t0 - alpha*delta_b ;  E = E_pi + sum_b (K/2 delta_b^2 + c delta_b),
    c = -2 alpha p_ref so that benzene (p = 2/3) is stationary at delta = 0.
    Units: eV, angstrom. delta_b is the change of bond b from the common reference length."""
    def __init__(self, g, t0=2.5, alpha=4.1, K=43.0, p_ref=2/3):
        self.T0, self.L, self.R, self.bix = bipartite(g)
        self.t0, self.a, self.K = t0, alpha, K
        self.c = -2*alpha*p_ref
        self.ii = np.array([i for i, j in self.bix]); self.jj = np.array([j for i, j in self.bix])
        self.m = len(self.bix)
    def Tw(self, d):
        T = np.zeros_like(self.T0); T[self.ii, self.jj] = self.t0 - self.a*d; return T
    def energy(self, d):
        U, s, Vt = np.linalg.svd(self.Tw(d))
        p = (U @ Vt)[self.ii, self.jj]                       # pi bond orders
        E = -2*s.sum() + np.sum(0.5*self.K*d**2 + self.c*d)
        grad = 2*self.a*p + self.K*d + self.c
        return E, grad
    def gap(self, d): return 2*np.linalg.svd(self.Tw(d), compute_uv=False)[-1]
    def bond_orders(self, d):
        U, s, Vt = np.linalg.svd(self.Tw(d)); return (U @ Vt)[self.ii, self.jj]
    def det(self, d):
        """det of the weighted block, normalized by t0^n, with its gradient (adjugate via SVD)."""
        A = self.Tw(d)/self.t0
        U, s, Vt = np.linalg.svd(A)
        sgnU, sgnV = np.linalg.det(U), np.linalg.det(Vt)
        D = sgnU*sgnV*np.prod(s)
        cof = np.array([np.prod(np.delete(s, k)) for k in range(len(s))])
        adjT = sgnU*sgnV*(U*cof) @ Vt                       # d det / dA  (= adj(A)^T)
        dDdd = adjT[self.ii, self.jj]*(-self.a/self.t0)
        return D, dDdd
    def relax(self, d0=None):
        d0 = np.zeros(self.m) if d0 is None else d0
        r = minimize(self.energy, d0, jac=True, method='BFGS', options=dict(gtol=1e-9, maxiter=5000))
        return r.x, r.fun
    def closure(self, starts, d_ref):
        """Lowest-energy geometry with det = 0 (Hückel gap exactly closed)."""
        best = None
        bnd = [(-0.8*self.t0/self.a, 0.8*self.t0/self.a)]*self.m
        for d0 in starts:
            r = minimize(self.energy, d0, jac=True, method='SLSQP', bounds=bnd,
                         constraints=[dict(type='eq', fun=lambda d: self.det(d)[0],
                                           jac=lambda d: self.det(d)[1])],
                         options=dict(ftol=1e-12, maxiter=2000))
            if r.success and abs(self.det(r.x)[0]) < 1e-8 and (best is None or r.fun < best.fun):
                best = r
        return best

# ---------- drawing ----------
def draw_bonds(ax, g, values=None, cmap=None, norm=None, lw=3.2, color='#334155', labels=False):
    P = np.array(g['pos'])
    for k, (u, v) in enumerate(g['edges']):
        c = color if values is None else cmap(norm(values[k]))
        ax.plot(*P[[u, v]].T, color=c, lw=lw, solid_capstyle='round', zorder=2)
    ax.scatter(P[:, 0], P[:, 1], s=14, color='#1f2937', zorder=3)
    if labels:
        for i, (x, y) in enumerate(P): ax.text(x, y, str(i), fontsize=6, ha='center', va='center')
    ax.set_aspect('equal'); ax.axis('off')
