"""Efficiency of the bond-length patterns in <name>_points.json. Usable as a module or script."""
import json, numpy as np

HART = 27.211386245988

def load(name, path='.'):
    R = json.load(open(f'{path}/{name}_points.json'))
    ref = R['ref']
    for lab, r in R.items():
        r['dE'] = (r['E'] - ref['E'])*HART
        r['dgap'] = ref['gap'] - r['gap']
        r['eta'] = r['dgap']/r['dE'] if r['dE'] > 1e-6 else np.nan
    return R

def efficiency(R):
    """eta = gap decrease per eV of B3LYP deformation energy at model energy 1.5 eV: predicted pattern,
    bond-equalizing control, and the median of the eight random patterns (each with its better sign)."""
    etaP = R['P@1.5']['eta']; etaB = R['BLA@1.5']['eta'] if 'BLA@1.5' in R else float('nan')
    best = [max(R[f'R{i}+@1.5']['eta'], R[f'R{i}-@1.5']['eta']) for i in range(8)
            if f'R{i}+@1.5' in R and f'R{i}-@1.5' in R]
    med = float(np.median(best)) if best else np.nan
    return dict(eta_P=etaP, eta_BLA=etaB, eta_rand_best=best, eta_rand_median=med,
                ratio_random=etaP/med if med > 0 else np.inf, ratio_BLA=etaP/etaB)

if __name__ == '__main__':
    for nm in ('biphenylene', 'angular3'):
        print(nm, json.dumps(efficiency(load(nm)), default=float))
