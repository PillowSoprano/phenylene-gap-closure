"""Adiabatic T1 (UB3LYP-D3-free, 6-31G*) for small phenylenes built from dataset coordinates."""
import json, sys, time, numpy as np
from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize
def build(coords, bonds):
    X=np.array(coords)*1.40; n=len(X)
    deg=np.zeros(n,int); nb={i:[] for i in range(n)}
    for a,b in bonds: deg[a]+=1; deg[b]+=1; nb[a].append(b); nb[b].append(a)
    atoms=[('C',(x,y,0.0)) for x,y in X]
    for i in range(n):
        if deg[i]==2:
            v=2*X[i]-X[nb[i][0]]-X[nb[i][1]]; v/=np.linalg.norm(v)
            h=X[i]+1.08*v; atoms.append(('H',(h[0],h[1],0.0)))
    return atoms
def run(atoms, spin, basis='6-31g*'):
    mol=gto.M(atom=atoms,basis=basis,spin=spin,verbose=0)
    mf=dft.UKS(mol).density_fit(); mf.xc='b3lyp'; mf.conv_tol=1e-9
    mol_eq=optimize(mf,maxsteps=100)
    mf2=dft.UKS(mol_eq).density_fit(); mf2.xc='b3lyp'; mf2.conv_tol=1e-10
    e=mf2.kernel()
    s2=mf2.spin_square()[0]
    geo=[(mol_eq.atom_symbol(i),tuple(mol_eq.atom_coord(i,unit='Angstrom'))) for i in range(mol_eq.natm)]
    return e,s2,geo
if __name__=='__main__':
    D={d['id']:d for d in json.load(open('patterns_N2-9.json'))}; C=json.load(open('coords_N2-9.json'))
    out=[]
    for cid in sys.argv[1:]:
        d=D[cid]; bonds=[(a,b) for a,b,_ in d['cp']]; atoms=build(C[cid],bonds)
        t=time.time(); e0,_,geo0=run(atoms,0); e1,s2,_=run(geo0,2)
        r=dict(id=cid,la=d['la_tree'],N=d['N'],rho=d['rho'],ET1_eV=(e1-e0)*27.211386,S2_T1=s2,sec=time.time()-t)
        print(r,flush=True); out.append(r)
        json.dump(out,open('dft_t1_%s.json'%('_'.join(sys.argv[1:])[:40]),'w'))
