import numpy as np, json
from chem import *
out = {}
for name, ch in [('biphenylene', []), ('angular3', ['L'])]:
    g = phenylene(ch); bm = BondModel(g, K=43.0)
    d, E = bm.relax(); U, s, Vt = np.linalg.svd(bm.Tw(d))
    assert abs(s[-1]-s[-2]) > 1e-6
    a = -2*bm.a*U[:, -1][bm.ii]*Vt[-1][bm.jj]
    h = 1e-5; eye = np.eye(bm.m)
    F = np.column_stack([(bm.energy(d+h*eye[k])[1]-bm.energy(d-h*eye[k])[1])/(2*h) for k in range(bm.m)]); F = (F+F.T)/2
    gfd = np.array([(bm.gap(d+h*v)-bm.gap(d-h*v))/(2*h) for v in eye]); assert np.abs(gfd-a).max() < 1e-6
    unit = lambda x: x/np.sqrt(0.5*x@F@x)          # scaled so that model energy 1/2 x^T F x = 1 eV
    P = unit(-np.linalg.solve(F, a))
    BLA = unit(-d)
    rng = np.random.default_rng(2026)
    R = [unit(rng.normal(size=bm.m)) for _ in range(8)]
    # first-order model gap change per sqrt(eV) for reference
    out[name] = dict(P=P.tolist(), BLA=BLA.tolist(), R=[r.tolist() for r in R], edges=g['edges'],
                     model_gap_slope=dict(P=float(a@P), BLA=float(a@BLA), R=[float(a@r) for r in R]))
    print(name, 'model dgap per unit (eV^0.5): P', round(a@P,3), 'BLA', round(a@BLA,3), 'rand', np.round([a@r for r in R],3))
    print('  bond change pattern P at 1 eV (A):', np.round(P,3))
json.dump(out, open('directions.json', 'w'), indent=1)
