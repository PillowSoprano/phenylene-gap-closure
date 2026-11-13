"""Compute the data plotted in Figs. S2 and S3 of the Supporting Information and write it to fig_data/*.csv.
plot_figs.py reads these files and does no further calculation for these figures."""
import sys, os, json, csv
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from kekule_core import phenylene, triangular4, bipartite, rho_rohn

os.chdir(Path(__file__).resolve().parent)
os.makedirs('fig_data', exist_ok=True)
def write(name, header, rows):
    with open(f'fig_data/{name}', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header)
        for r in rows: w.writerow([f'{x:.10g}' if isinstance(x, (float, np.floating)) else x for x in r])

MOLS = {'biphenylene': phenylene([]), 'linear3': phenylene(['p']), 'angular3': phenylene(['L']),
        'linear4': phenylene(['p', 'p']), 'bent4': phenylene(['p', 'L']), 'zigzag4': phenylene(['L', 'L']),
        'triangular4': triangular4()}

# ---------------- Fig. S2: closing patterns and closing orbitals ----------------
summary = []
for key in ('linear3', 'angular3', 'linear4', 'bent4', 'zigzag4', 'triangular4'):
    g = MOLS[key]; T, L, R, bix = bipartite(g)
    r, (y, z) = rho_rohn(T)
    W = T - r*(T*np.outer(y, z)); U, s, Vt = np.linalg.svd(W)
    w = np.zeros(g['n']); w[L] = U[:, -1]**2; w[R] = Vt[-1]**2; w /= 2
    write(f'closing_{key}_atoms.csv', ['atom', 'x', 'y', 'closing_orbital_weight'],
          [(a, *g['pos'][a], w[a]) for a in range(g['n'])])
    write(f'closing_{key}_bonds.csv', ['atom1', 'atom2', 'coupling_change'],
          [(u, v, W[bix[k]] - 1) for k, (u, v) in enumerate(g['edges'])])
    summary.append((key, r, s[-1]))
write('closing_summary.csv', ['skeleton', 'rho', 'sigma_min_at_pattern'], summary)

# ---------------- Fig. S3: B3LYP gap decrease against deformation energy ----------------
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
write('dft_ladders.csv', ['molecule', 'pattern', 'model_energy_eV', 'b3lyp_deformation_energy_eV', 'b3lyp_gap_eV',
                           'gap_decrease_eV', 'y0', 'S2_broken_symmetry'], lad)
write('dft_random.csv', ['molecule', 'label', 'b3lyp_deformation_energy_eV', 'b3lyp_gap_eV', 'gap_decrease_eV'], rnd)
print('wrote', sorted(os.listdir('fig_data')))
