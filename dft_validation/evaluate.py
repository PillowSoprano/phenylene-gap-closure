"""Apply the locked criteria of PREREGISTRATION.md to <name>_points.json. Usable as a module or script."""
import json, numpy as np

HART = 27.211386245988
NAMES = {'biphenylene': 'biphenylene', 'angular3': 'angular[3]phenylene'}

def load(name, path='.'):
    R = json.load(open(f'{path}/{name}_points.json'))
    ref = R['ref']
    for lab, r in R.items():
        r['dE'] = (r['E'] - ref['E'])*HART
        r['dgap'] = ref['gap'] - r['gap']
        r['eta'] = r['dgap']/r['dE'] if r['dE'] > 1e-6 else np.nan
        if 'E_BS' in r:
            r['E_S'] = min(r['E'], r['E_BS'])
            r['dEST'] = (r['E_T'] - r['E_S'])*HART
            r['dBS'] = (r['E_BS'] - r['E'])*HART
    return R

def criterion1(R):
    etaP = R['P@1.5']['eta']; etaB = R['BLA@1.5']['eta'] if 'BLA@1.5' in R else float('nan')
    best = []
    for i in range(8):
        a, b = R.get(f'R{i}+@1.5'), R.get(f'R{i}-@1.5')
        if a is None or b is None: continue
        best.append(max(a['eta'], b['eta']))
    med = float(np.median(best)) if best else np.nan
    return dict(eta_P=etaP, eta_BLA=etaB, eta_rand_best=best, eta_rand_median=med,
                ratio=etaP/med if med > 0 else np.inf, n_rand=len(best), passed=bool(len(best) == 8 and etaP >= 1.5*med))

def ladder(R, kind='P'):
    pts = [R['ref']] + [R[k] for k in sorted((k for k in R if k.startswith(kind + '@')), key=lambda k: R[k]['Em'])]
    return pts

def criterion2(Ra, Rb):
    """Ra = biphenylene, Rb = angular[3]; along each molecule's own P ladder."""
    La, Lb = ladder(Ra), ladder(Rb)
    ea = np.array([p['dE'] for p in La]); ya = np.array([p['y0'] for p in La])
    eb = np.array([p['dE'] for p in Lb]); yb = np.array([p['y0'] for p in Lb])
    assert np.all(np.diff(ea) > 0) and np.all(np.diff(eb) > 0), 'dE not monotone along ladder'
    top = min(ea[-1], eb[-1]); grid = np.arange(0, top + 1e-9, 0.25)
    ia, ib = np.interp(grid, ea, ya), np.interp(grid, eb, yb)
    c2a = bool(np.all(ia > ib))
    def onset(L):
        for p in L:
            if p['S2_BS'] > 0.1: return p['dE']
        return None
    oa, ob = onset(La), onset(Lb)
    if oa is None and ob is None: c2b = 'inconclusive'
    elif oa is not None and ob is not None: c2b = 'consistent' if oa < ob else 'contradicts'
    else: c2b = 'consistent' if oa is not None else 'contradicts'
    return dict(grid=grid.tolist(), y0_biph=ia.tolist(), y0_ang=ib.tolist(), c2a=c2a, onset_biph=oa, onset_ang=ob,
                c2b=c2b, passed=bool(c2a and c2b != 'contradicts'))

def verdict(c1a, c1b, c2):
    if not (c1a['passed'] and c1b['passed']): return '判据 1 不通过：结论停在模型层面，归档。'
    if not c2['passed']: return '判据 1 通过、判据 2 不通过：只写机制，不写反转。'
    return '判据 1、2 都通过：可以写化学故事（机制 + 反转）。'

if __name__ == '__main__':
    Ra, Rb = load('biphenylene'), load('angular3')
    for nm, R in [('biphenylene', Ra), ('angular3', Rb)]:
        print(nm, json.dumps(criterion1(R), default=float))
    print(json.dumps(criterion2(Ra, Rb), default=float))
