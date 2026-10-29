import sys, json; from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import *
mols={'benzene':hexagon(),'naphthalene':naphthalene(),'benzocyclobutadiene':benzocyclobutadiene(),'biphenylene':phenylene([]),'linear[3]':phenylene(['p']),'angular[3]':phenylene(['L']),
 'linear[4]':phenylene(['p','p']),'bent[4]':phenylene(['p','L']),'helical[4]':phenylene(['L','R']),'zigzag[4]':phenylene(['L','L']),'triangular[4]':triangular4()}
rhos=json.load(open(ROOT/'results'/'thresholds.json'))
out={}
for name,g in mols.items():
    T,L,R,bix=bipartite(g); n=len(T); Ti=np.linalg.inv(T); sv=np.linalg.svd(T,compute_uv=False)
    L1=sv[-1]/sv[0]
    L2=1/max(abs(np.linalg.eigvals(np.abs(Ti)@T)))
    U1=min(1.0,1/max(np.diag(T@np.abs(Ti))))
    terms=kekule_terms(T)
    sets=[(s,frozenset((i,p[i]) for i in range(n))) for s,p in terms]
    U2=1.0; best=None
    for s0,M0 in sets:
        if s0>0: continue   # sign normalized so det T > 0; minority = negative
        # f(q) = sum_M s_M q^{|M \ M0|}; f(1)>0, f(0)=-1; largest root in (0,1)
        coef=np.zeros(n+1)
        for s,M in sets: coef[len(M-M0)]+=s
        r=np.roots(coef[::-1]); r=[x.real for x in r if abs(x.imag)<1e-9 and 0<x.real<1]
        q=max(r); e=(1-q)/(1+q)
        if e<U2: U2=e
    nminus=sum(1 for s,_ in sets if s<0)
    rho=rhos[name]['rho']
    out[name]=dict(n=n,minority=nminus,L_gap=L1,L_BS=L2,rho=rho,U_Rohn=U1,U_Kek=U2,Rump=(3+2*2**.5)*n*L2)
    print(f"{name:20s} minus={nminus:2d} gap-bound={L1:.4f} BS={L2:.4f} rho={rho:.4f} Kek={U2:.4f} Rohn={U1:.4f} Rump={(3+2*2**.5)*n*L2:.1f}")
json.dump(out,open(ROOT/'results'/'bounds.json','w'),indent=1)
