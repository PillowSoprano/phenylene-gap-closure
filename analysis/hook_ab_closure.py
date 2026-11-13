"""Closing pattern of the two hook molecules (Fig. 1): where the gap closes and where the triplet relaxes.

A = PT09-00009 ((()(A()))(A())(A())),  B = PT09-00036 ((A(()()))(L())(L())).
Needs the gapclosure package (../gapclosure/src) and beyond_huckel/patterns_N2-9.json.
Usage: python hook_ab_closure.py <path to gapclosure/src> <path to patterns_N2-9.json>
"""
import sys, json, itertools
import numpy as np
sys.path.insert(0, sys.argv[1])
from gapclosure.skeleton import _Builder
from gapclosure.threshold import closure_threshold

D = {d['id']: d for d in json.load(open(sys.argv[2]))}
MOLS = {'A': 'PT09-00009', 'B': 'PT09-00036'}


def build(links):
    """phenylene_tree that also returns the atoms of each hexagon and the parent of each hexagon"""
    b = _Builder(); ring = {0: list(range(6))}; parent = {0: None}; port = {}
    for i, a, j in links:
        c = ring[i]
        sq = b.fuse(c[a % 6], c[(a + 1) % 6], 4)
        hx = b.fuse(sq[2], sq[1], 6)
        ring[j] = [hx[5], hx[0], hx[1], hx[2], hx[3], hx[4]]
        parent[j] = i; port.setdefault(i, []).append(a % 6); port.setdefault(j, []).append(0)
    return b.skeleton(''), ring, parent, port


def ring_type(ports):
    """L/A label of an inner hexagon with two four-membered rings, B for three"""
    if len(ports) == 2:
        d = (ports[1] - ports[0]) % 6
        return 'L' if d == 3 else 'A'
    return {1: 'end', 3: 'branch'}.get(len(ports), '?')


# ---------------- Hueckel-SSH relaxation (Kwapisz-Stolarczyk) ----------------
R0, X, Y, RREF = 1.523, 0.189, 0.2756, 1.397


def relax(bonds, n, occ, R):
    for _ in range(20000):
        H = np.zeros((n, n))
        for (i, j) in bonds:
            H[i, j] = H[j, i] = -(1 + (RREF - R[(i, j)]) / Y)
        w, v = np.linalg.eigh(H); oc = np.zeros(n); oc[:n // 2] = 2
        if occ == 'T1':
            oc[n // 2 - 1] = 1; oc[n // 2] = 1
        p = {e: float(np.sum(oc * v[e[0]] * v[e[1]])) for e in bonds}
        Rn = {e: R0 - X * p[e] for e in bonds}
        if max(abs(Rn[e] - R[e]) for e in bonds) < 1e-11:
            break
        R = {e: 0.6 * R[e] + 0.4 * Rn[e] for e in bonds}
    E = float(np.sum(oc * w)) + sum((R[e] - R0) ** 2 for e in bonds) / (2 * X * Y)
    return R, E


out = {}
for name, cid in MOLS.items():
    d = D[cid]
    sk, ring, parent, port = build(d['links'])
    res = closure_threshold(sk)
    owner = {a: h for h, atoms in ring.items() for a in atoms}
    w_ring = {h: float(sum(res.orbital[a] for a in atoms)) for h, atoms in ring.items()}
    types = {h: ring_type(port.get(h, [])) for h in ring}
    # pattern: share of raised/lowered bonds per hexagon is not needed; record bonds inside each hexagon
    bonds = [tuple(sorted(b)) for b in sk.bonds]
    # relaxed ground state and lowest triplet (symmetric start + 4 perturbed starts)
    rng = np.random.default_rng(1)
    RS, ES = relax(bonds, sk.n, 'S0', {e: RREF for e in bonds})
    best = None
    for k in range(5):
        start = {e: RS[e] + (0 if k == 0 else rng.normal(0, 0.02)) for e in bonds}
        RT, ET = relax(bonds, sk.n, 'T1', start)
        if best is None or ET < best[1] - 1e-9:
            best = (RT, ET)
    RT, ET = best
    dR = {e: RT[e] - RS[e] for e in bonds}
    tot = sum(abs(v) for v in dR.values())
    # assign each bond change to hexagons: a hexagon bond to its hexagon, a four-ring bond half to each end
    t_ring = {h: 0.0 for h in ring}
    for (u, v), x in dR.items():
        hu, hv = owner[u], owner[v]
        if hu == hv:
            t_ring[hu] += abs(x)
        else:
            t_ring[hu] += abs(x) / 2; t_ring[hv] += abs(x) / 2
    t_ring = {h: x / tot for h, x in t_ring.items()}
    # sign agreement of triplet bond changes with the closing pattern (coupling up <-> bond shorter)
    pat = {tuple(sorted(b)): s for b, s in res.pattern.items()}
    agree = [e for e in bonds if pat.get(e, 0) != 0 and np.sign(-dR[e]) == pat[e]]
    out[name] = dict(id=cid, la=d['la_tree'], rho=res.rho, gap=res.gap, K=res.K, minority=res.minority,
                     ring_type=types, parent=parent, closing_orbital_by_ring=w_ring,
                     triplet_change_by_ring=t_ring, ET1_hssh=ET - ES,
                     sign_agree=[len(agree), len(bonds), float(sum(abs(dR[e]) for e in agree) / tot)])
    print(f"\n{name} = {cid} {d['la_tree']}: rho {res.rho:.4f}, gap {res.gap:.4f}, K {res.K}, minority {res.minority}, "
          f"HSSH ET1 {ET - ES:.4f} |beta|")
    print(f"  triplet bond changes with the sign of the closing pattern: {len(agree)}/{len(bonds)} bonds, "
          f"{100 * out[name]['sign_agree'][2]:.0f}% of the total change")
    print('  hexagon  parent  type    closing orbital  triplet change')
    for h in sorted(ring, key=lambda h: -w_ring[h]):
        print(f"  {h:6d}  {str(parent[h]):>6}  {types[h]:6s}  {w_ring[h]:15.3f}  {t_ring[h]:14.3f}")
json.dump(out, open('hook_ab_closure.json', 'w'), indent=1, default=str)
