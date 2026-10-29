"""B3LYP/def2-SVP gradient of the Kohn-Sham HOMO-LUMO gap with respect to the C-C bond lengths, by central
finite differences at the B3LYP equilibrium geometry. Each bond-length change is realized with the same
minimum-norm Wilson-B displacement as in run_points.py (hydrogens follow their carbons)."""
import sys, os, json, time, numpy as np
from pyscf import gto, dft, lib
from chem import build, atom_str
lib.num_threads(int(os.environ.get('NT', '2')))
HART = 27.211386245988
name, ch = sys.argv[1], json.loads(sys.argv[2]); h = float(os.environ.get('H', '0.005'))
g, sym, _ = build(ch); nC = g['n']; X0 = np.load(f'{name}_eq.npy')
H_of = {k: int(np.argmin(np.linalg.norm(X0[:nC]-X0[k], axis=1))) for k in range(nC, len(sym))}
def bl(X): return np.array([np.linalg.norm(X[u]-X[v]) for u, v in g['edges']])
r_eq = bl(X0)
def displace(dr, iters=8):
    X = X0.copy(); target = r_eq + dr
    for _ in range(iters):
        C = X[:nC, :2]; r = bl(X); B = np.zeros((len(g['edges']), 2*nC))
        for k, (u, v) in enumerate(g['edges']):
            e = (C[u]-C[v])/r[k]; B[k, 2*u:2*u+2] = e; B[k, 2*v:2*v+2] = -e
        dx = np.linalg.lstsq(B, target-r, rcond=None)[0].reshape(nC, 2); X[:nC, :2] += dx
        for hk, ck in H_of.items(): X[hk, :2] += dx[ck]
    return X, float(np.abs(bl(X)-target).max())
def gap(X, dm0=None):
    mol = gto.M(atom=atom_str(sym, X), basis='def2-svp', verbose=0)
    mf = dft.RKS(mol).density_fit(); mf.xc = 'b3lyp'; mf.conv_tol = 1e-10
    mf.kernel(dm0=dm0); no = mol.nelectron//2
    return (mf.mo_energy[no]-mf.mo_energy[no-1])*HART, mf.make_rdm1(), bool(mf.converged)
outf = f'{name}_gapgrad.json'
res = json.load(open(outf)) if os.path.exists(outf) else {}
t = time.time(); g0, dm0, ok = gap(X0); res['ref'] = dict(gap=g0, converged=ok)
print(f'{name}: reference gap {g0:.5f} eV ({time.time()-t:.0f}s)', flush=True)
m = len(g['edges'])
for k in range(m):
    for sgn in (+1, -1):
        key = f'{k}{"+" if sgn > 0 else "-"}'
        if key in res: continue
        t = time.time(); dr = np.zeros(m); dr[k] = sgn*h
        X, err = displace(dr); gp, _, ok = gap(X, dm0)
        res[key] = dict(gap=gp, maxerr=err, converged=ok); json.dump(res, open(outf, 'w'), indent=1)
        print(f'  bond {k:2d} {g["edges"][k]} {key[-1]}  gap {gp:.6f}  maxerr {err:.1e}  ({time.time()-t:.0f}s)', flush=True)
grad = [(res[f'{k}+']['gap']-res[f'{k}-']['gap'])/(2*h) for k in range(m)]
res['grad_eV_per_A'] = grad; res['h_A'] = h; json.dump(res, open(outf, 'w'), indent=1)
print('done', name, np.round(grad, 3))
