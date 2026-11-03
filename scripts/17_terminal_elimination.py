"""Terminal elimination (Corollary 1, a consequence of the gluing formula) and the factorization of det T for
linear and angular [3]phenylene (Proposition 3): exact symbolic checks, the effective couplings
at the closing patterns, and the mirror symmetry that exchanges the two terminal rings."""
import sys, json, itertools
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np, sympy as sp, networkx as nx
from gapclosure import phenylene, benzocyclobutadiene, bipartite, kekule_terms

def poly(g, atoms, x):
    """Signed Kekule polynomial of the subgraph on `atoms`, in the bond variables x of g
    (overall sign fixed so that the value at equal couplings is positive)."""
    atoms = sorted(atoms); loc = {a: k for k, a in enumerate(atoms)}
    E = [(k, e) for k, e in enumerate(g['edges']) if e[0] in loc and e[1] in loc]
    if not E: return sp.Integer(1) if not atoms else sp.Integer(0)
    sub = dict(n=len(atoms), edges=[(loc[u], loc[v]) for _, (u, v) in E])
    T, L, R, bix = bipartite(sub); var = {bix[i]: x[k] for i, (k, _) in enumerate(E)}
    return sp.expand(sum(s * sp.prod([var[(i, j)] for i, j in enumerate(p)]) for s, p in kekule_terms(T)))

def eliminate(g, A, x):
    """Split g into A and B = rest, joined by two connecting bonds ac, bd. Returns P_A, the bond h = cd
    of B, and kappa_h = x_h - x_ac x_bd (dP_A/dx_ab) / P_A."""
    E = {frozenset(e): k for k, e in enumerate(g['edges'])}
    conn = [e for e in g['edges'] if (e[0] in A) != (e[1] in A)]; assert len(conn) == 2
    (a, c), (b, d) = [(u, v) if u in A else (v, u) for u, v in conn]
    PA = poly(g, A, x); xab = x[E[frozenset((a, b))]]; h = E[frozenset((c, d))]
    kappa = x[h] - x[E[frozenset((a, c))]] * x[E[frozenset((b, d))]] * sp.diff(PA, xab) / PA
    return PA, h, kappa

log = {}
# 1. Corollary: P_G = P_A * P_B(x_h -> kappa_h), checked on every phenylene with N <= 5 (first hexagon as A)
print('Corollary 1 (terminal elimination), exact symbolic check:')
cases = [('benzocyclobutadiene', benzocyclobutadiene())] + \
        [(f'[{N}]phenylene {"".join(ch) or "-"}', phenylene(list(ch))) for N in (2, 3, 4, 5) for ch in itertools.product('pLR', repeat=N - 2)]
for name, g in cases:
    x = sp.symbols(f'x0:{len(g["edges"])}'); A = set(g['rings'][0]); B = set(range(g['n'])) - A
    PG = poly(g, range(g['n']), x); PA, h, kappa = eliminate(g, A, x); PB = poly(g, B, x)
    rhs = sp.expand(sp.cancel(PA * PB.subs(x[h], kappa)))
    ok = sp.expand(PG - rhs) == 0 or sp.expand(PG + rhs) == 0
    assert ok, name
print(f'  holds for benzocyclobutadiene and all {len(cases) - 1} phenylene constructions with N <= 5')
log['corollary_cases'] = len(cases)

# 2. [3]phenylenes: eliminate both terminal hexagons; det T = A_L A_R P_H(kappa_L, kappa_R)
store = json.load(open(ROOT / 'results' / 'bond_model.json'))
for name, ch in (('linear[3]', ['p']), ('angular[3]', ['L'])):
    g = phenylene(ch); x = sp.symbols(f'x0:{len(g["edges"])}')
    hexes = [r for r in g['rings'] if len(r) == 6]; HL, HM, HR = hexes
    PL, hL, kL = eliminate(g, set(HL), x); PR, hR, kR = eliminate(g, set(HR), x)
    kL_, kR_ = sp.symbols('kappa_L kappa_R'); PH = poly(g, HM, x)
    reduced = sp.expand(PH.subs({x[hL]: kL_, x[hR]: kR_}))
    full = sp.expand(sp.cancel(PL * PR * reduced.subs({kL_: kL, kR_: kR})))
    PG = poly(g, range(g['n']), x); assert sp.expand(PG - full) == 0 or sp.expand(PG + full) == 0
    same = sp.Poly(reduced, kL_, kR_).coeff_monomial(kL_ * kR_) != 0
    print(f'{name}: det T = A_L A_R * ({reduced})   [attachment bonds in one Kekule structure of the middle ring: {same}]')
    T, Lb, Rb, bix = bipartite(g); s = store[name]; y, z = np.array(s['y']), np.array(s['z'])
    w = {x[k]: 1 - s['rho'] * y[i] * z[j] for k, (i, j) in enumerate(bix)}
    vals = dict(rho=s['rho'], kappa_L=float(kL.subs(w)), kappa_R=float(kR.subs(w)), kappa_equal=float(kL.subs({v: 1 for v in x})))
    if name == 'angular[3]':
        rest = sp.expand(reduced - sp.Poly(reduced, kL_, kR_).coeff_monomial(kL_ * kR_) * kL_ * kR_)
        tau = sp.Poly(reduced, kL_, kR_).coeff_monomial(kL_ * kR_)
        vals['Pi_over_tau'] = float((rest / tau).subs(w)); vals['kappa_L_kappa_R'] = vals['kappa_L'] * vals['kappa_R']
        # mirror symmetry exchanging the terminal rings maps A_L<->A_R, kappa_L<->kappa_R and fixes P and q
        G = nx.Graph(g['edges']); auts = [m for m in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter()
                                          if set(m[a] for a in HL) == set(HR)]
        Eix = {frozenset(e): k for k, e in enumerate(g['edges'])}
        perm = {x[k]: x[Eix[frozenset((auts[0][u], auts[0][v]))]] for k, (u, v) in enumerate(g['edges'])}
        mirror_ok = (sp.expand(PL.subs(perm, simultaneous=True) - PR) == 0 and sp.simplify(kL.subs(perm, simultaneous=True) - kR) == 0
                     and sp.expand(rest.subs(perm, simultaneous=True) - rest) == 0 and sp.expand(tau.subs(perm, simultaneous=True) - tau) == 0)
        vals['mirror_maps_L_to_R_and_fixes_Pi_tau'] = bool(mirror_ok); assert mirror_ok
    else:
        # all 16 singular corners of linear [3]: kappa_L = kappa_R = 0
        cc = json.load(open(ROOT / 'results' / 'closing_corners.json'))[name]
        free = {frozenset(e) for e in cc['free_bonds']}; ks = []
        for bits in itertools.product((-1, 1), repeat=4):
            ww = dict(w); it = iter(bits)
            for k, e in enumerate(g['edges']):
                if frozenset(e) in free: ww[x[k]] = 1 + next(it) * s['rho']
            ks.append((float(kL.subs(ww)), float(kR.subs(ww)), float(PG.subs(ww))))
        vals['max_abs_kappa_on_face'] = max(max(abs(a), abs(b)) for a, b, _ in ks); vals['max_abs_det_on_face'] = max(abs(c) for *_, c in ks)
    print('   ', {k: (round(v, 10) if isinstance(v, float) else v) for k, v in vals.items()})
    log[name] = dict(reduced=str(reduced), **vals)
json.dump(log, open(ROOT / 'results' / 'terminal_elimination.json', 'w'), indent=1)
