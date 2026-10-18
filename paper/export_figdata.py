"""Compute every quantity that appears in Figs. 1-3 and write it to fig_data/*.csv.
The plotting script (plot_figs.py) reads only these files.
Calculations are the same as in make_figs.py; nothing here changes the data."""
import sys, os, json, csv
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import phenylene, triangular4, bipartite, rho_rohn

os.chdir(Path(__file__).resolve().parent)
os.makedirs('fig_data', exist_ok=True)
def write(name, header, rows):
    with open(f'fig_data/{name}', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header)
        for r in rows: w.writerow([f'{x:.10g}' if isinstance(x, (float, np.floating)) else x for x in r])

MOLS = {'biphenylene': phenylene([]), 'linear3': phenylene(['p']), 'angular3': phenylene(['L']),
        'linear4': phenylene(['p', 'p']), 'bent4': phenylene(['p', 'L']), 'zigzag4': phenylene(['L', 'L']),
        'triangular4': triangular4()}

# ---------------- Fig 1 ----------------
g = MOLS['biphenylene']; T, L, R, bix = bipartite(g)
r, (y, z) = rho_rohn(T)
sign = np.array([-y[i]*z[j] for i, j in bix])           # +1 strengthened, -1 weakened at the threshold
Ti = np.linalg.inv(T)
crit = [k for k, (i, j) in enumerate(bix) if Ti[j, i] < 0][0]
write('fig1_atoms.csv', ['atom', 'x', 'y'], [(a, *g['pos'][a]) for a in range(g['n'])])
write('fig1_bonds.csv', ['atom1', 'atom2', 'cooperative_sign', 'single_bond_route'],
      [(u, v, int(sign[k]), int(k == crit)) for k, (u, v) in enumerate(g['edges'])])
def gap_single(t):
    U = T.copy(); i, j = bix[crit]; U[i, j] = 1 + t; return 2*np.linalg.svd(U, compute_uv=False)[-1]
def gap_coop(e):
    U = T.copy()
    for (i, j), s_ in zip(bix, sign): U[i, j] = 1 + e*s_
    return 2*np.linalg.svd(U, compute_uv=False)[-1]
xs = np.geomspace(0.005, 3.0, 400); es = np.geomspace(0.005, r, 400)
write('fig1_curves.csv', ['route', 'largest_relative_change_percent', 'huckel_gap_abs_beta'],
      [('single', 100*t, gap_single(t)) for t in xs] + [('cooperative', 100*e, gap_coop(e)) for e in es])
write('fig1_summary.csv', ['quantity', 'value'], [('rho_biphenylene', r), ('single_bond_relative_increase', 3.0)])

# ---------------- Fig 2 ----------------
summary = []
for key in ('linear3', 'angular3', 'linear4', 'bent4', 'zigzag4', 'triangular4'):
    g = MOLS[key]; T, L, R, bix = bipartite(g)
    r, (y, z) = rho_rohn(T)
    W = T - r*(T*np.outer(y, z)); U, s, Vt = np.linalg.svd(W)
    w = np.zeros(g['n']); w[L] = U[:, -1]**2; w[R] = Vt[-1]**2; w /= 2
    write(f'fig2_{key}_atoms.csv', ['atom', 'x', 'y', 'closing_orbital_weight'],
          [(a, *g['pos'][a], w[a]) for a in range(g['n'])])
    write(f'fig2_{key}_bonds.csv', ['atom1', 'atom2', 'coupling_change'],
          [(u, v, W[bix[k]] - 1) for k, (u, v) in enumerate(g['edges'])])
    summary.append((key, r, s[-1]))
write('fig2_summary.csv', ['skeleton', 'rho', 'sigma_min_at_pattern'], summary)

# ---------------- Fig 3 ----------------
HART = 27.211386245988
lad, rnd = [], []
for key in ('biphenylene', 'angular3'):
    Rj = json.load(open(ROOT/'dft_validation'/f'{key}_points.json')); ref = Rj['ref']
    for lab, v in Rj.items():
        dE = (v['E']-ref['E'])*HART; dgap = ref['gap']-v['gap']
        if lab == 'ref' or lab.startswith(('P@', 'BLA@')):
            kinds = ['P', 'BLA'] if lab == 'ref' else [lab.split('@')[0]]
            for kind in kinds:
                lad.append((key, kind, v['Em'], dE, v['gap'], dgap,
                            v.get('y0', ''), v.get('S2_BS', '')))
        else:
            rnd.append((key, lab, dE, v['gap'], dgap))
lad.sort(key=lambda t: (t[0], t[1], t[2]))
write('fig3_ladders.csv', ['molecule', 'pattern', 'model_energy_eV', 'b3lyp_deformation_energy_eV', 'b3lyp_gap_eV',
                           'gap_decrease_eV', 'y0', 'S2_broken_symmetry'], lad)
write('fig3_random.csv', ['molecule', 'label', 'b3lyp_deformation_energy_eV', 'b3lyp_gap_eV', 'gap_decrease_eV'], rnd)
print('wrote', sorted(os.listdir('fig_data')))
