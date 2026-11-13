"""Closed-gap geometry nearest to the relaxed Hueckel-SSH ground state (Section 5), N = 4 to 8.

The nearest point t = tS + g with det T(t) = 0 satisfies g proportional to u_i v_j, the product of the
left and right null vectors of T on the two atoms of each bond. A damped fixed-point iteration on this
condition is started from the rescaled closing pattern of rho, from the steepest-descent direction of
the gap at tS, and from NSTART randomly perturbed copies of both; the closest result is kept. Roots of
det T along a direction are found exactly as eigenvalues. Writes l2nearest.npy with columns
N, E_T1, U, U_corner, rho, g0, |g|  (energies in |beta|). Usage: NSTART=30 python l2nearest.py
(With IDS='PT08-00001 ...' only those classes are done and their points are written to l2nearest_points.json.)"""
import json, os, numpy as np
R0, x, y, Rref = 1.523, 0.189, 0.2756, 1.397; K = 2/(x*y)
NSTART = int(os.environ.get('NSTART', '30')); rng = np.random.default_rng(1108)
t_of_R = lambda R: 1.0 + (Rref - R)/y
R_of_t = lambda t: Rref - y*(t - 1.0)
here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, 'patterns_N2-9.json')))
H = {r['id']: r for r in json.load(open(os.path.join(here, 't1_scan.json')))}

def first_root(M0, M1):
    """smallest |c| with det(M0 + c M1) = 0, signed"""
    lam = np.linalg.eigvals(np.linalg.solve(M0, M1))
    c = [-1/l.real for l in lam if abs(l.imag) < 1e-9*max(1, abs(l)) and abs(l.real) > 1e-14]
    return min(c, key=abs) if c else None

IDS = set(os.environ.get('IDS', '').split()) - {''}
out = []; points = {}
for d in D:
    if d['N'] > 8 or (IDS and d['id'] not in IDS): continue
    cp = d['cp']; E = [(a, b) for a, b, _ in cp]; s = np.array([c for *_, c in cp], float); n = 6*d['N']; nh = n//2
    ia = np.array([a for a, _ in E]); ib = np.array([b for _, b in E])
    col = -np.ones(n, int); col[0] = 0; adj = {i: [] for i in range(n)}
    for a, b in E: adj[a].append(b); adj[b].append(a)
    st = [0]
    while st:
        u = st.pop()
        for v in adj[u]:
            if col[v] < 0: col[v] = 1 - col[u]; st.append(v)
    A = [v for v in range(n) if col[v] == 0]; B = [v for v in range(n) if col[v] == 1]
    ai = {v: i for i, v in enumerate(A)}; bi = {v: i for i, v in enumerate(B)}
    rr = np.array([ai[a] if col[a] == 0 else ai[b] for a, b in E]); cc = np.array([bi[b] if col[a] == 0 else bi[a] for a, b in E])
    def T(t): M = np.zeros((nh, nh)); M[rr, cc] = t; return M
    def ener(R):
        t = t_of_R(R); Hm = np.zeros((n, n)); Hm[ia, ib] = -t; Hm[ib, ia] = -t
        e, C = np.linalg.eigh(Hm); w = np.zeros(n); w[:nh] = 2
        return float(w @ e) + K/2*np.sum((R - R0)**2), C, w
    R = np.full(len(E), Rref)
    for it in range(3000):
        _, C, w = ener(R); p = (C[ia]*C[ib]) @ w; Rn = R0 - x*p
        if abs(Rn - R).max() < 1e-12: break
        R = 0.5*R + 0.5*Rn
    RS = R; ES = ener(RS)[0]; tS = t_of_R(RS); TS = T(tS)
    def project(dirn):
        dirn = dirn/np.linalg.norm(dirn); c = first_root(TS, T(dirn))
        return None if c is None else c*dirn
    def nearest(g):
        for it in range(500):
            Ut, S_, Vt = np.linalg.svd(T(tS + g)); dirn = Ut[rr, -1]*Vt[-1, cc]
            gn = project(dirn)
            if gn is None: return None
            if gn @ g < 0: gn = -gn if project(-dirn) is None else project(-dirn)
            if np.linalg.norm(gn - g) < 1e-11: return gn
            g = 0.5*g + 0.5*gn
        return g
    tc = 1 + d['rho']*s; tc = tc*(tc @ tS)/(tc @ tc)
    Ut, S_, Vt = np.linalg.svd(TS); gd = Ut[rr, -1]*Vt[-1, cc]
    seeds = [tc - tS, project(gd), project(-gd)]
    for _ in range(NSTART):
        seeds.append(project(tc - tS + rng.normal(0, .03, len(tc))))
        seeds.append(project(gd/np.linalg.norm(gd) + rng.normal(0, .3, len(tc))))
    gs = [nearest(g) for g in seeds if g is not None]
    g = min((g for g in gs if g is not None), key=np.linalg.norm)
    assert np.linalg.svd(T(tS + g), compute_uv=False)[-1] < 1e-8
    points[d['id']] = dict(edges=E, tS=tS.tolist(), t_near=(tS + g).tolist())
    Ec = ener(R_of_t(tS + g))[0]; Ecorner = ener(R_of_t(tc))[0]
    out.append((d['N'], H[d['id']]['ET1_adiabatic'], Ec - ES, Ecorner - ES, d['rho'], H[d['id']]['gapS0'], float(np.linalg.norm(g))))
A = np.array(out)
if IDS:
    json.dump(points, open('l2nearest_points.json', 'w')); print(A[:, 2])
else:
    np.save(os.path.join(here, 'l2nearest.npy'), A)
print('classes', len(A), ' bound E_T1 <= U holds:', bool(np.all(A[:, 1] <= A[:, 2] + 1e-9)))
