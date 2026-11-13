"""Columns of the threshold table (Table S2 of the Supporting Information), written to results/table1.json."""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from kekule_core import (hexagon, naphthalene, benzocyclobutadiene, phenylene, triangular4, bipartite,
                        structure_counts)

MOLS = {'benzene': hexagon(), 'naphthalene': naphthalene(), 'benzocyclobutadiene': benzocyclobutadiene(),
        'biphenylene': phenylene([]), 'linear [3]': phenylene(['p']), 'angular [3]': phenylene(['L']),
        'linear [4]': phenylene(['p', 'p']), 'bent [4]': phenylene(['p', 'L']),
        'helical [4]': phenylene(['L', 'R']), 'zigzag [4]': phenylene(['L', 'L']),
        'triangular [4]': triangular4()}

# exact thresholds from scripts/01_thresholds.py (Rohn criterion); values within 1e-4 of 1 mean "never closes"
_thr = json.load(open(ROOT/'results'/'thresholds.json'))
# rounded to 12 digits so that exact ties (helical = zigzag [4]phenylene, both of class AA) stay ties;
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

# ---------------- Table S2 ----------------
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
