"""Tables 1 and 2 of the manuscript (written to results/table1.json and results/table2.json)."""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from scipy.stats import spearmanr
from gapclosure import (hexagon, naphthalene, benzocyclobutadiene, phenylene, triangular4, bipartite,
                        structure_counts, BondModel)

MOLS = {'benzene': hexagon(), 'naphthalene': naphthalene(), 'benzocyclobutadiene': benzocyclobutadiene(),
        'biphenylene': phenylene([]), 'linear [3]': phenylene(['p']), 'angular [3]': phenylene(['L']),
        'linear [4]': phenylene(['p', 'p']), 'bent [4]': phenylene(['p', 'L']),
        'helical [4]': phenylene(['L', 'R']), 'zigzag [4]': phenylene(['L', 'L']),
        'triangular [4]': triangular4()}

# exact thresholds from scripts/01_thresholds.py (Rohn criterion); values within 1e-4 of 1 mean "never closes"
_thr = json.load(open(ROOT/'results'/'thresholds.json'))
# rounded to 12 digits so that exact ties (helical = zigzag [4]phenylene, Section 3 of the paper) stay ties;
# unrounded they differ in the 16th digit, which would break the tie in rank correlations
RHO = {name: (1.0 if _thr[name.replace(' ', '')]['rho'] > 1 - 1e-4 else round(_thr[name.replace(' ', '')]['rho'], 12)) for name in MOLS}

def single_bond_increase(T, bix):
    """Smallest relative increase of one coupling that closes the gap (rank-one formula); inf if none."""
    Ti = np.linalg.inv(T)
    ts = [-1/Ti[j, i] for i, j in bix if Ti[j, i] < -1e-12]
    return min(ts) if ts else np.inf

def single_bond_any(T, bix):
    """Smallest |relative change| of one coupling (increase, or decrease that keeps it >= 0) that closes the gap."""
    Ti = np.linalg.inv(T)
    ts = [abs(-1/Ti[j, i]) for i, j in bix if abs(Ti[j, i]) > 1e-12 and -1/Ti[j, i] >= -1]
    return min(ts) if ts else np.inf

# ---------------- Table 1 ----------------
rows = []
for name, g in MOLS.items():
    T, L, R, bix = bipartite(g)
    K, ASC, neg = structure_counts(T)
    sv = np.linalg.svd(T, compute_uv=False)
    rows.append(dict(name=name, atoms=g['n'], bonds=len(g['edges']), K=K, ASC=int(ASC), neg=neg, gap=2*sv[-1],
                     rho=RHO[name], single_up=single_bond_increase(T, bix), single_any=single_bond_any(T, bix)))
json.dump(rows, open(ROOT/'results'/'table1.json', 'w'), indent=1, default=float)
for r in rows:
    print(f"{r['name']:20s} {r['atoms']:3d} {r['bonds']:3d} K={r['K']:3d} ASC={r['ASC']:2d} neg={r['neg']:2d} "
          f"gap={r['gap']:.3f} rho={r['rho']:.4f} single_up={r['single_up']:.3f} single_any={r['single_any']:.3f}")

# ---------------- Table 2: closure energies in the bond-length model ----------------
ssh = json.load(open(ROOT/'results'/'bond_model.json'))
mapn = {'benzocyclobutadiene': 'benzocyclobutadiene', 'biphenylene': 'biphenylene', 'linear [3]': 'linear[3]',
        'angular [3]': 'angular[3]', 'linear [4]': 'linear[4]', 'bent [4]': 'bent[4]', 'helical [4]': 'helical[4]',
        'zigzag [4]': 'zigzag[4]', 'triangular [4]': 'triangular[4]'}
t2 = []
for name, key in mapn.items():
    bm = BondModel(MOLS[name], K=43.0); d, E = bm.relax(); U, s, Vt = np.linalg.svd(bm.Tw(d)); g0 = 2*s[-1]
    eye = np.eye(bm.m); h = 1e-5
    F = np.column_stack([(bm.energy(d+h*eye[k])[1]-bm.energy(d-h*eye[k])[1])/(2*h) for k in range(bm.m)]); F = (F+F.T)/2
    Fi = np.linalg.inv(F)
    deg = int(np.sum(abs(s-s[-1]) < 1e-7))
    ths = [0.0] if deg == 1 else np.linspace(0, np.pi, 721)
    best = None
    for th in ths:
        vv = np.ones(1) if deg == 1 else np.array([np.cos(th), np.sin(th)])
        l = U[:, -deg:] @ vv; rr = Vt[-deg:].T @ vv
        a = -2*bm.a*l[bm.ii]*rr[bm.jj]; c = a @ Fi @ a; Pv = np.sum((l[bm.ii]*rr[bm.jj])**2)
        if best is None or c > best[0]: best = (c, Pv)
    t2.append(dict(name=name, g0=g0, P=best[1], E_hess=g0**2/(2*best[0]), g2P=g0**2/best[1],
                   dE=ssh[key]['K43']['dE_close'], lam=(float(np.linalg.eigvalsh(F).min()), float(np.linalg.eigvalsh(F).max()))))
json.dump(t2, open(ROOT/'results'/'table2.json', 'w'), indent=1)
dE = [x['dE'] for x in t2]
print('spearman E_hess', spearmanr(dE, [x['E_hess'] for x in t2]).statistic,
      'g2P', spearmanr(dE, [x['g2P'] for x in t2]).statistic, 'g0', spearmanr(dE, [x['g0'] for x in t2]).statistic)
for x in t2: print(x)

