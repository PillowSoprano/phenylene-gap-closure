import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np, json
from scipy.stats import spearmanr
from gapclosure import phenylene, triangular4, benzocyclobutadiene, bipartite, BondModel
t2 = json.load(open(ROOT/'results'/'table2.json')); t1 = {r['name']: r for r in json.load(open(ROOT/'results'/'table1.json'))}
names = [x['name'] for x in t2]; dE = [x['dE'] for x in t2]
def sp(v): r = spearmanr(dE, v); return f'{r.statistic:.3f} (p = {r.pvalue:.1e})'
print('Spearman with closure energy (9 skeletons):')
print('  g0^2/P          ', sp([x['g2P'] for x in t2]))
print('  local estimate  ', sp([x['E_hess'] for x in t2]))
print('  rho             ', sp([t1[n]['rho'] for n in names]))
print('  uniform gap     ', sp([t1[n]['gap'] for n in names]))
print('  relaxed gap g0  ', sp([x['g0'] for x in t2]))
# norms of the two routes for biphenylene
r = 0.188828952032
print(f'\nbiphenylene: L_inf 3.00 vs {r:.3f}; L1 3.00 vs {14*r:.2f}; sum of squares 9.00 vs {14*r*r:.2f}')
t0, al = 2.5, 4.1
print(f'bond-length change for x4 coupling: {3*t0/al:.2f} A ; for 18.9%: {r*t0/al:.3f} A')
# single-bond route energy in the model (from relaxed geometry, one inter-ring bond only)
g = phenylene([]); bm = BondModel(g, K=43.0); d, E0 = bm.relax()
k = [i for i, e in enumerate(g['edges']) if set(e) in ({2, 6}, {7, 3})][0]
from scipy.optimize import brentq
f = lambda s: bm.det(np.where(np.arange(bm.m) == k, d[k] + s, d))[0]
ss = np.linspace(0, -3, 3001); vals = [f(s) for s in ss]; i = np.where(np.diff(np.sign(vals)))[0]
if len(i):
    s0 = brentq(f, ss[i[0]], ss[i[0]+1]); dd = d.copy(); dd[k] += s0
    print(f'model single-bond closure from relaxed geometry: bond shortened by {-s0:.2f} A, energy +{bm.energy(dd)[0]-E0:.1f} eV (vs 4.18 eV cooperative)')
# model gap crossing along x* (local approximation)
G = {x['name']: x for x in t2}
def cross(a, b):
    ga, Ca = G[a]['g0'], G[a]['g0']**2/(2*G[a]['E_hess']); gb, Cb = G[b]['g0'], G[b]['g0']**2/(2*G[b]['E_hess'])
    return ((ga-gb)/(np.sqrt(2*Ca)-np.sqrt(2*Cb)))**2
print(f'model (local) gap crossing biphenylene / angular[3]: {cross("biphenylene", "angular [3]"):.2f} eV')
# isotropic baseline for the efficiency ratio: median over random directions of 1/|cos|, best of two signs is |cos|
rng = np.random.default_rng(0)
for dname, dd_ in (('biphenylene', 14), ('angular [3]', 22)):
    v = rng.normal(size=(200000, dd_)); c = np.abs(v[:, 0])/np.linalg.norm(v, axis=1)
    print(f'isotropic ceiling d={dd_}: median 1/|cos| over 8-sample medians ~ {np.median([1/np.median(c[j*8:(j+1)*8]) for j in range(20000)]):.2f}; 1/median|cos| = {1/np.median(c):.2f}')
# model-predicted ratio for the actual 8 random directions (linear gap response, equal model energy)
D = json.load(open(ROOT/'dft_validation'/'directions.json'))
for key, ch in (('biphenylene', []), ('angular3', ['L'])):
    g = phenylene(ch); bm = BondModel(g, K=43.0); d, E = bm.relax(); U, s, Vt = np.linalg.svd(bm.Tw(d))
    a = -2*bm.a*U[:, -1][bm.ii]*Vt[-1][bm.jj]
    P = np.array(D[key]['P']); R = [np.array(x) for x in D[key]['R']]
    rat = abs(a @ P)/np.median([abs(a @ x) for x in R])
    print(f'{key}: model-predicted eta_P / median random = {rat:.2f}  (B3LYP measured: {"5.48" if key=="biphenylene" else "6.59"})')
