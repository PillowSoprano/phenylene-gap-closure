import sys, time, numpy as np, json
from pyscf import gto, dft, lib
from pyscf.geomopt.berny_solver import optimize
from chem import *
lib.num_threads(2)
name, ch = sys.argv[1], json.loads(sys.argv[2])
g, sym, xyz = build(ch)
mol = gto.M(atom=atom_str(sym, xyz), basis='def2-svp', verbose=0)
mf = dft.RKS(mol).density_fit(); mf.xc = 'b3lyp'
t = time.time()
mol_eq = optimize(mf, maxsteps=60)
X = mol_eq.atom_coords(unit='Angstrom')
np.save(f'{name}_eq.npy', X)
mf2 = dft.RKS(mol_eq).density_fit(); mf2.xc = 'b3lyp'; e = mf2.kernel()
mo_e = mf2.mo_energy; nocc = mol_eq.nelectron//2
print(name, 'optimized in', round(time.time()-t), 's; E =', e, '; gap (eV) =', (mo_e[nocc]-mo_e[nocc-1])*27.211386, flush=True)
P = X[:g['n']]
print('C-C bond lengths:', [round(float(np.linalg.norm(P[u]-P[v])), 4) for u, v in g['edges']])
print('max |z|:', float(abs(X[:, 2]).max()))
