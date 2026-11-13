"""Cheapest closed-gap geometry (minimum ground-state energy on det T = 0) for the classes in
l2nearest_points.json, by SLSQP from the nearest closed-gap point and from perturbed copies; compared
with U, the energy at the nearest point (Section 5)."""
import json, sys, numpy as np
from ssh_core import Model
rng = np.random.default_rng(1108)
for i, p in json.load(open(sys.argv[1] if len(sys.argv) > 1 else '../beyond_huckel/l2nearest_points.json')).items():
    m = Model([tuple(e) for e in p['edges']]); ts, rs = m.relax(); es = rs['energy']
    tn = np.array(p['t_near']); assert np.allclose(ts, p['tS'], atol=1e-7)
    best = None
    for st in [tn] + [tn + rng.normal(0, .02, len(tn)) for _ in range(6)]:
        t, r = m.close(ts, st, energy=True)
        if r['success'] and (best is None or r['energy'] < best[1]['energy']): best = (t, r)
    tdist, rd = m.close(ts, tn)
    print(i, 'U (nearest) %.4f  cheapest %.4f  lower by %.1f%%' % (rd['energy'] - es, best[1]['energy'] - es, 100*(1 - (best[1]['energy'] - es)/(rd['energy'] - es))))
