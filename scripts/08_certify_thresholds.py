import sys, json, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import phenylene, triangular4, benzocyclobutadiene, bipartite, kekule_terms
mols={'benzocyclobutadiene':benzocyclobutadiene(),'biphenylene':phenylene([]),'linear[3]':phenylene(['p']),'angular[3]':phenylene(['L']),
 'linear[4]':phenylene(['p','p']),'zigzag[4]':phenylene(['L','L']),'helical[4]':phenylene(['L','R']),
 'bent[4]':phenylene(['p','L']),'triangular[4]':triangular4()}
rhos=json.load(open(ROOT/'results'/'thresholds.json'))
def mindet(g,eps):
    T,L,R,bix=bipartite(g); terms=kekule_terms(T); m=len(bix); eix={e:k for k,e in enumerate(bix)}
    TE=[(s,[eix[(i,p[i])] for i in range(len(T))]) for s,p in terms]; h=m//2
    def side(lo,hi):
        k=hi-lo; idx=np.arange(1<<k); S=((idx[:,None]>>np.arange(k))&1)*2-1
        return np.array([np.prod(1+eps*S[:,[e-lo for e in ix if lo<=e<hi]],axis=1) for s,ix in TE]).T
    A=side(0,h)*np.array([s for s,_ in TE]); B=side(h,m); mn=np.inf
    for a0 in range(0,A.shape[0],1<<10): mn=min(mn,(A[a0:a0+(1<<10)]@B.T).min())
    return mn, len(TE), len(T)
d=1e-8
for name,g in mols.items():
    r=rhos[name]['rho']; lo,K,n=mindet(g,r*(1-d)); hi,_,_=mindet(g,r*(1+d))
    bound=K*(1+r)**n*(n+2)*2**-52
    print(f'{name}: rho={r:.10f}  min det at rho(1-1e-8)={lo:.3e}  at rho(1+1e-8)={hi:.3e}  rounding bound {bound:.1e}')
